"""Streaming FASTA and FASTQ parsing."""

from __future__ import annotations

from collections.abc import Iterable, Iterator
from dataclasses import dataclass


@dataclass(frozen=True, slots=True)
class SequenceRecord:
    """One named biological sequence."""

    name: str
    sequence: str


_SPREADSHEET_FORMULA_PREFIXES = frozenset("=+-@")


def _nonempty_lines(lines: Iterable[str]) -> Iterator[tuple[int, str]]:
    for number, raw_line in enumerate(lines, 1):
        line = raw_line.rstrip("\r\n")
        if not line:
            raise ValueError(f"line {number}: blank lines are not allowed")
        yield number, line


def parse_records(lines: Iterable[str]) -> Iterator[SequenceRecord]:
    """Parse FASTA or four-line FASTQ incrementally, retaining one record at a time."""

    source = iter(_nonempty_lines(lines))
    try:
        first_number, first = next(source)
    except StopIteration:
        raise ValueError("input is empty") from None
    if first.startswith(">"):
        yield from _parse_fasta(source, first_number, first)
    elif first.startswith("@"):
        yield from _parse_fastq(source, first_number, first)
    else:
        raise ValueError(f"line {first_number}: expected a FASTA or FASTQ header")


def _name(number: int, header: str) -> str:
    name = header[1:].strip()
    if not name:
        raise ValueError(f"line {number}: record name is empty")
    if any(ord(character) < 32 or ord(character) == 127 for character in name):
        raise ValueError(f"line {number}: record name contains a control character")
    if name[0] in _SPREADSHEET_FORMULA_PREFIXES:
        raise ValueError(f"line {number}: record name has an unsafe formula prefix")
    return name


def _parse_fasta(
    source: Iterator[tuple[int, str]], first_number: int, first: str
) -> Iterator[SequenceRecord]:
    name = _name(first_number, first)
    pieces: list[str] = []
    for number, line in source:
        if line.startswith(">"):
            if not pieces:
                raise ValueError(f"line {number}: previous FASTA record has no sequence")
            yield SequenceRecord(name, "".join(pieces))
            name = _name(number, line)
            pieces = []
        elif line.startswith(("@", "+")):
            raise ValueError(f"line {number}: unexpected FASTQ marker in FASTA sequence")
        else:
            pieces.append(line.strip())
    if not pieces:
        raise ValueError("final FASTA record has no sequence")
    yield SequenceRecord(name, "".join(pieces))


def _parse_fastq(
    source: Iterator[tuple[int, str]], first_number: int, first: str
) -> Iterator[SequenceRecord]:
    header_number, header = first_number, first
    while True:
        name = _name(header_number, header)
        try:
            sequence_number, sequence = next(source)
            plus_number, plus = next(source)
            quality_number, quality = next(source)
        except StopIteration:
            raise ValueError(f"line {header_number}: incomplete FASTQ record") from None
        if sequence.startswith(("@", "+", ">")):
            raise ValueError(f"line {sequence_number}: invalid FASTQ sequence")
        if not plus.startswith("+"):
            raise ValueError(f"line {plus_number}: expected '+' separator")
        repeated_name = plus[1:].strip()
        if repeated_name and repeated_name != name:
            raise ValueError(f"line {plus_number}: FASTQ names do not match")
        if len(sequence) != len(quality):
            raise ValueError(f"line {quality_number}: sequence and quality lengths differ")
        yield SequenceRecord(name, sequence)
        try:
            header_number, header = next(source)
        except StopIteration:
            return
        if not header.startswith("@"):
            raise ValueError(f"line {header_number}: expected a FASTQ header")
