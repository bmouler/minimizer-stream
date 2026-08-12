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


def test_parser_preserves_spaces_and_reports_physical_line_numbers() -> None:
    assert list(parse_records([">one\n", " AC \n"])) == [SequenceRecord("one", "AC")]
    with pytest.raises(ValueError, match=r"^line 2: blank lines are not allowed$"):
        list(parse_records([">one\n", "\n"]))
    with pytest.raises(ValueError, match=r"^line 1: record name is empty$"):
        list(parse_records([">\n", "AC\n"]))


def test_multiline_fasta_joins_without_separators() -> None:
    assert list(parse_records([">one\n", "AC\n", "GT\n"])) == [SequenceRecord("one", "ACGT")]


def test_del_control_character_is_rejected() -> None:
    with pytest.raises(ValueError, match=r"^line 1: record name contains a control character$"):
        list(parse_records([">\x7f\n", "AC\n"]))


@pytest.mark.parametrize("marker", ["+", ">"])
def test_fastq_sequence_rejects_all_reserved_markers(marker: str) -> None:
    with pytest.raises(ValueError, match=r"^line 2: invalid FASTQ sequence$"):
        list(parse_records(["@one\n", f"{marker}A\n", "+\n", "!!\n"]))


def test_fastq_sequence_preserves_trailing_spaces() -> None:
    assert list(parse_records(["@one\n", "AC \n", "+\n", "!!!\n"])) == [
        SequenceRecord("one", "AC ")
    ]


def test_sequence_preserves_terminal_x_before_newline() -> None:
    assert list(parse_records([">one\n", "ACX\n"])) == [SequenceRecord("one", "ACX")]


def test_errors_report_the_originating_header_line() -> None:
    with pytest.raises(ValueError, match=r"^line 3: record name is empty$"):
        list(parse_records([">one\n", "AC\n", ">\n", "GT\n"]))
    with pytest.raises(ValueError, match=r"^line 1: record name is empty$"):
        list(parse_records(["@\n", "AC\n", "+\n", "!!\n"]))


def test_exact_terminal_record_errors() -> None:
    with pytest.raises(ValueError, match=r"^input is empty$"):
        list(parse_records([]))
    with pytest.raises(ValueError, match=r"^final FASTA record has no sequence$"):
        list(parse_records([">one\n"]))


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
