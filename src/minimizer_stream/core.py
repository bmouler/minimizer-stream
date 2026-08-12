"""Canonical minimizer algorithms."""

from __future__ import annotations

from collections import deque
from collections.abc import Iterable, Iterator
from dataclasses import dataclass

_MASK64 = (1 << 64) - 1
_COMPLEMENT = str.maketrans("ACGT", "TGCA")
_VALID = frozenset("ACGT")


@dataclass(frozen=True, slots=True)
class Minimizer:
    """A selected canonical k-mer."""

    position: int
    hash: int
    strand: str
    kmer: str


def stable_hash64(value: str) -> int:
    """Return the deterministic 64-bit FNV-1a hash of an ASCII string."""

    result = 0xCBF29CE484222325
    for byte in value.encode("ascii"):
        result ^= byte
        result = (result * 0x100000001B3) & _MASK64
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
        if forward <= reverse:
            canonical, strand = forward, "+"
        else:
            canonical, strand = reverse, "-"
        yield Minimizer(position - k + 1, stable_hash64(canonical), strand, canonical)


def minimizers(
    sequence: Iterable[str], k: int, w: int, *, deduplicate: bool = False
) -> Iterator[Minimizer]:
    """Stream canonical minimizers in O(n*k) time using a monotonic deque.

    ``w`` is the number of consecutive valid k-mers in each minimizer window.
    Windows never cross an invalid base. Hash ties select the leftmost k-mer.
    """

    _validate(k, w)
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
