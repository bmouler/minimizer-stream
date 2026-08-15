"""Canonical minimizer algorithms."""

from __future__ import annotations

from collections import deque
from collections.abc import Iterable, Iterator
from dataclasses import dataclass
from typing import Literal

_MASK64 = (1 << 64) - 1
_FNV_OFFSET = 0xCBF29CE484222325
_FNV_PRIME = 0x100000001B3
_FNV4_PRIME = pow(_FNV_PRIME, 4, 1 << 64)
_COMPLEMENT = str.maketrans("ACGT", "TGCA")
_VALID = frozenset("ACGT")
_BASE_CODES = {"A": 0, "a": 0, "C": 1, "c": 1, "G": 2, "g": 2, "T": 3, "t": 3}
_BASE_BYTES = (65, 67, 71, 84)
_KMER4 = tuple(
    first + second + third + fourth
    for first in "ACGT"
    for second in "ACGT"
    for third in "ACGT"
    for fourth in "ACGT"
)


def _build_fnv4_corrections() -> tuple[tuple[int, ...], ...]:
    """Build four-base FNV transitions indexed by packed bases and low hash bits."""

    rows = []
    for packed in range(256):
        bases = tuple(_BASE_BYTES[(packed >> shift) & 3] for shift in (6, 4, 2, 0))
        row = []
        for low in range(128):
            result = low
            for byte in bases:
                result = ((result ^ byte) * _FNV_PRIME) & _MASK64
            row.append((result - low * _FNV4_PRIME) & _MASK64)
        rows.append(tuple(row))
    return tuple(rows)


_FNV4_CORRECTIONS = _build_fnv4_corrections()


@dataclass(frozen=True, slots=True)
class Minimizer:
    """A selected canonical k-mer."""

    position: int
    hash: int
    strand: Literal["+", "-"]
    kmer: str


def stable_hash64(value: str) -> int:
    """Return the deterministic 64-bit FNV-1a hash of an ASCII string."""

    result = _FNV_OFFSET
    encoded = value.encode("ascii")
    offset = 0
    end = len(encoded) - 3
    while offset < end:
        result = (result ^ encoded[offset]) * _FNV_PRIME
        result = (result ^ encoded[offset + 1]) * _FNV_PRIME
        result = (result ^ encoded[offset + 2]) * _FNV_PRIME
        result = ((result ^ encoded[offset + 3]) * _FNV_PRIME) & _MASK64
        offset += 4
    while offset < len(encoded):
        result = ((result ^ encoded[offset]) * _FNV_PRIME) & _MASK64
        offset += 1
    return result


def _validate(k: int, w: int) -> None:
    if k <= 0:
        raise ValueError("k must be greater than zero")
    if w <= 0:
        raise ValueError("w must be greater than zero")


def _canonical_kmers(sequence: Iterable[str], k: int) -> Iterator[Minimizer | None]:
    """Yield canonical k-mers, with ``None`` marking an invalid-base reset."""

    run: deque[str] = deque(maxlen=k)
    position = -1
    for raw_base in sequence:
        position += 1
        base = raw_base.upper()
        if base not in _VALID:
            run.clear()
            yield None
            continue
        run.append(base)
        if len(run) < k:
            continue
        forward = "".join(run)
        reverse = forward.translate(_COMPLEMENT)[::-1]
        strand: Literal["+", "-"]
        if forward <= reverse:
            canonical, strand = forward, "+"
        else:
            canonical, strand = reverse, "-"
        yield Minimizer(position - k + 1, stable_hash64(canonical), strand, canonical)


def _decode_kmer21(packed: int) -> str:
    return (
        _KMER4[(packed >> 34) & 255]
        + _KMER4[(packed >> 26) & 255]
        + _KMER4[(packed >> 18) & 255]
        + _KMER4[(packed >> 10) & 255]
        + _KMER4[(packed >> 2) & 255]
        + "ACGT"[packed & 3]
    )


def _minimizers_str21(sequence: str, w: int, deduplicate: bool) -> Iterator[Minimizer]:
    """Fast constant-space path for the documented 21-base string workload."""

    candidates: deque[tuple[int, int, int, Literal["+", "-"], int]] = deque()
    forward = 0
    reverse = 0
    valid_bases = 0
    run_index = 0
    previous: tuple[int, int, int, Literal["+", "-"], int] | None = None
    mask = (1 << 42) - 1
    base_codes = _BASE_CODES
    corrections = _FNV4_CORRECTIONS
    fnv_offset = _FNV_OFFSET
    fnv_prime = _FNV_PRIME
    fnv4_prime = _FNV4_PRIME
    mask64 = _MASK64
    base_bytes = _BASE_BYTES
    decode = _decode_kmer21
    minimizer_type = Minimizer

    for position, base in enumerate(sequence):
        code = base_codes.get(base)
        if code is None:
            candidates.clear()
            forward = 0
            reverse = 0
            valid_bases = 0
            run_index = 0
            continue
        forward = ((forward << 2) | code) & mask
        reverse = (reverse >> 2) | ((3 - code) << 40)
        valid_bases += 1
        if valid_bases < 21:
            continue

        if forward <= reverse:
            canonical = forward
            strand: Literal["+", "-"] = "+"
        else:
            canonical = reverse
            strand = "-"

        result = fnv_offset
        chunk = (canonical >> 34) & 255
        result = (result * fnv4_prime + corrections[chunk][result & 127]) & mask64
        chunk = (canonical >> 26) & 255
        result = (result * fnv4_prime + corrections[chunk][result & 127]) & mask64
        chunk = (canonical >> 18) & 255
        result = (result * fnv4_prime + corrections[chunk][result & 127]) & mask64
        chunk = (canonical >> 10) & 255
        result = (result * fnv4_prime + corrections[chunk][result & 127]) & mask64
        chunk = (canonical >> 2) & 255
        result = (result * fnv4_prime + corrections[chunk][result & 127]) & mask64
        result = ((result ^ base_bytes[canonical & 3]) * fnv_prime) & mask64

        while candidates and candidates[-1][1] > result:
            candidates.pop()
        candidate = (run_index, result, position - 20, strand, canonical)
        candidates.append(candidate)
        window_start = run_index - w + 1
        while candidates[0][0] < window_start:
            candidates.popleft()
        if window_start >= 0:
            selected = candidates[0]
            if not deduplicate or previous != selected:
                yield minimizer_type(selected[2], selected[1], selected[3], decode(selected[4]))
                previous = selected
        run_index += 1


def minimizers(
    sequence: Iterable[str], k: int, w: int, *, deduplicate: bool = False
) -> Iterator[Minimizer]:
    """Stream canonical minimizers in O(n*k) time using a monotonic deque.

    ``w`` is the number of consecutive valid k-mers in each minimizer window.
    Windows never cross an invalid base. Hash ties select the leftmost k-mer.
    """

    _validate(k, w)
    if type(sequence) is str and k == 21:
        yield from _minimizers_str21(sequence, w, deduplicate)
        return
    candidates: deque[tuple[int, Minimizer]] = deque()
    run_index = 0
    previous: Minimizer | None = None
    for item in _canonical_kmers(sequence, k):
        if item is None:
            candidates.clear()
            run_index = 0
            continue
        while candidates and candidates[-1][1].hash > item.hash:
            candidates.pop()
        candidates.append((run_index, item))
        window_start = run_index - w + 1
        while candidates[0][0] < window_start:
            candidates.popleft()
        if window_start >= 0:
            selected = candidates[0][1]
            if not deduplicate or previous != selected:
                yield selected
                previous = selected
        run_index += 1


def naive_minimizers(
    sequence: Iterable[str], k: int, w: int, *, deduplicate: bool = False
) -> Iterator[Minimizer]:
    """Reference O(n*(k + w)) minimizer implementation for validation."""

    _validate(k, w)
    window: deque[Minimizer] = deque(maxlen=w)
    previous: Minimizer | None = None
    for item in _canonical_kmers(sequence, k):
        if item is None:
            window.clear()
            continue
        window.append(item)
        if len(window) == w:
            selected = min(window, key=lambda candidate: candidate.hash)
            if not deduplicate or previous != selected:
                yield selected
                previous = selected
