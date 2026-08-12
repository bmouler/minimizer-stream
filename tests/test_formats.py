from __future__ import annotations

import pytest

from minimizer_stream import SequenceRecord, parse_records


def test_multiline_fasta_is_parsed_record_by_record() -> None:
    consumed: list[str] = []

    def lines() -> object:
        for line in [">first\n", "AC\n", "GT\n", ">second description\n", "NN\n"]:
            consumed.append(line)
            yield line

    records = parse_records(lines())  # type: ignore[arg-type]
    assert next(records) == SequenceRecord("first", "ACGT")
    assert consumed == [">first\n", "AC\n", "GT\n", ">second description\n"]
    assert list(records) == [SequenceRecord("second description", "NN")]


def test_fastq_with_optional_repeated_name() -> None:
    text = "@one\nACGT\n+one\n!!!!\n@two\nTGCA\n+\n####\n"
    assert list(parse_records(text.splitlines(keepends=True))) == [
        SequenceRecord("one", "ACGT"),
        SequenceRecord("two", "TGCA"),
    ]


@pytest.mark.parametrize(
    ("text", "message"),
    [
        ("", "empty"),
        ("ACGT\n", "expected a FASTA or FASTQ header"),
        (">\nAC\n", "record name is empty"),
        (">one\n", "final FASTA record has no sequence"),
        (">one\n>two\nAC\n", "previous FASTA record has no sequence"),
        (">one\n@bad\n", "unexpected FASTQ marker"),
        (">one\n+bad\n", "unexpected FASTQ marker"),
        (">one\n\n", "blank lines"),
        ("@one\nAC\n+\n", "incomplete FASTQ record"),
        ("@one\n@bad\n+\n!!\n", "invalid FASTQ sequence"),
        ("@one\nAC\nwrong\n!!\n", r"expected '\+' separator"),
        ("@one\nAC\n+two\n!!\n", "names do not match"),
        ("@one\nAC\n+\n!\n", "lengths differ"),
        ("@one\nAC\n+\n!!\n>two\n", "expected a FASTQ header"),
        ("@\nAC\n+\n!!\n", "record name is empty"),
        (">one\textra\nAC\n", "control character"),
        (">=formula\nAC\n", "unsafe formula prefix"),
    ],
)
def test_malformed_input_is_rejected(text: str, message: str) -> None:
    with pytest.raises(ValueError, match=message):
        list(parse_records(text.splitlines(keepends=True)))
