from __future__ import annotations

import random

import pytest

from minimizer_stream import Minimizer, minimizers, naive_minimizers, stable_hash64


def test_stable_hash_known_vector() -> None:
    assert stable_hash64("") == 0xCBF29CE484222325
    assert stable_hash64("ACGT") == 0x9A90178BA8FEDA4E


def test_stable_hash_block_and_tail_vector() -> None:
    assert stable_hash64("ACGTA") == 0xA1B21150290C877D


def test_forward_and_reverse_complements_share_canonical_form() -> None:
    forward = list(minimizers("AAAC", 4, 1))
    reverse = list(minimizers("GTTT", 4, 1))
    assert forward == [Minimizer(0, stable_hash64("AAAC"), "+", "AAAC")]
    assert reverse == [Minimizer(0, stable_hash64("AAAC"), "-", "AAAC")]


def test_palindrome_uses_forward_strand() -> None:
    assert next(minimizers("ACGT", 4, 1)).strand == "+"


def test_leftmost_hash_tie_and_deduplication() -> None:
    selected = list(minimizers("AAAAA", 1, 3))
    assert [item.position for item in selected] == [0, 1, 2]
    assert [item.position for item in minimizers("AAAAA", 1, 3, deduplicate=True)] == [0, 1, 2]
    repeated = list(minimizers("ACGTA", 2, 2))
    deduplicated = list(minimizers("ACGTA", 2, 2, deduplicate=True))
    assert len(repeated) == 3
    assert deduplicated == [repeated[0], repeated[2]]


def test_invalid_base_resets_windows_and_positions_remain_absolute() -> None:
    result = list(minimizers("AAANAAAA", 2, 2))
    assert [item.position for item in result] == [0, 4, 5]
    assert list(minimizers("AANAA", 2, 2)) == []


def test_short_sequence_and_lowercase() -> None:
    assert list(minimizers("AC", 3, 1)) == []
    assert list(minimizers("acgt", 4, 1)) == list(minimizers("ACGT", 4, 1))


def test_k21_string_path_preserves_canonical_strands_and_resets() -> None:
    forward = list(minimizers("A" * 21, 21, 1))
    reverse = list(minimizers("T" * 21, 21, 1))
    assert forward == [Minimizer(0, stable_hash64("A" * 21), "+", "A" * 21)]
    assert reverse == [Minimizer(0, stable_hash64("A" * 21), "-", "A" * 21)]
    assert list(minimizers("acgt" * 6, 21, 2)) == list(naive_minimizers("acgt" * 6, 21, 2))
    assert list(minimizers("ACGT" * 5, 21, 1)) == []

    for invalid in ("N", "é"):
        sequence = "A" * 23 + invalid + "T" * 23
        result = list(minimizers(sequence, 21, 3))
        assert result == list(naive_minimizers(sequence, 21, 3))
        assert [item.position for item in result] == [0, 24]


def test_k21_string_path_preserves_ties_expiry_and_deduplication() -> None:
    ties = "A" * 25
    tied = list(minimizers(ties, 21, 3))
    assert tied == list(naive_minimizers(ties, 21, 3))
    assert [item.position for item in tied] == [0, 1, 2]

    sequence = "ACGGCAAGTGAGCAAGGTCGAAGC"
    repeated = list(minimizers(sequence, 21, 3))
    deduplicated = list(minimizers(sequence, 21, 3, deduplicate=True))
    assert repeated == list(naive_minimizers(sequence, 21, 3))
    assert deduplicated == list(naive_minimizers(sequence, 21, 3, deduplicate=True))
    assert len(repeated) == 2
    assert deduplicated == [repeated[0]]


def test_k21_non_string_iterable_retains_general_streaming_path() -> None:
    sequence = "GATTACA" * 5
    assert list(minimizers(iter(sequence), 21, 4)) == list(naive_minimizers(iter(sequence), 21, 4))


@pytest.mark.parametrize("function", [minimizers, naive_minimizers])
@pytest.mark.parametrize(
    ("k", "w", "message"),
    [
        (0, 1, "k must be greater than zero"),
        (1, 0, "w must be greater than zero"),
    ],
)
def test_invalid_parameters(function: object, k: int, w: int, message: str) -> None:
    with pytest.raises(ValueError, match=f"^{message}$"):
        list(function("ACGT", k, w))  # type: ignore[operator]


def test_naive_matches_streaming_across_deterministic_random_inputs() -> None:
    generator = random.Random(2026)
    for length in range(40):
        sequence = "".join(generator.choices("ACGTN", k=length))
        for k in range(1, 6):
            for w in range(1, 6):
                assert list(minimizers(sequence, k, w)) == list(naive_minimizers(sequence, k, w))
                assert list(minimizers(sequence, k, w, deduplicate=True)) == list(
                    naive_minimizers(sequence, k, w, deduplicate=True)
                )
