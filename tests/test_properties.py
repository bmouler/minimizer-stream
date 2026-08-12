from hypothesis import given
from hypothesis import strategies as st

from minimizer_stream import minimizers, naive_minimizers, stable_hash64

_SEQUENCES = st.text(alphabet="ACGTN", max_size=60)
_K = st.integers(min_value=1, max_value=6)
_W = st.integers(min_value=1, max_value=6)
_DEDUPLICATE = st.booleans()
_COMPLEMENT = str.maketrans("ACGT", "TGCA")


def _canonical(kmer: str) -> str:
    reverse_complement = kmer.translate(_COMPLEMENT)[::-1]
    return min(kmer, reverse_complement)


@given(seq=_SEQUENCES, k=_K, w=_W, deduplicate=_DEDUPLICATE)
def test_streaming_minimizers_equal_naive_oracle(
    seq: str, k: int, w: int, deduplicate: bool
) -> None:
    assert list(minimizers(seq, k, w, deduplicate=deduplicate)) == list(
        naive_minimizers(seq, k, w, deduplicate=deduplicate)
    )


@given(value=_SEQUENCES)
def test_stable_hash_is_deterministic_and_pinned(value: str) -> None:
    assert stable_hash64(value) == stable_hash64(value)
    assert stable_hash64("ACGT") == 0x9A90178BA8FEDA4E


@given(seq=_SEQUENCES, k=_K, w=_W, deduplicate=_DEDUPLICATE)
def test_emitted_minimizers_reverify(seq: str, k: int, w: int, deduplicate: bool) -> None:
    for item in minimizers(seq, k, w, deduplicate=deduplicate):
        assert item.kmer == _canonical(seq[item.position : item.position + k])
        assert item.hash == stable_hash64(item.kmer)
