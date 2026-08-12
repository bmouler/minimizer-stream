"""Command-line interface for minimizer-stream."""

from __future__ import annotations

import argparse
import json
import sys
from collections.abc import Callable, Iterable, Iterator
from contextlib import contextmanager
from pathlib import Path
from typing import Literal, Never, TextIO, TypeAlias, cast

from .core import Minimizer, minimizers
from .formats import parse_records

Formatter: TypeAlias = Callable[[str, Minimizer], str]
OutputFormat: TypeAlias = Literal["tsv", "jsonl"]


class _Arguments(argparse.Namespace):
    input: str
    k: int
    w: int
    format: OutputFormat
    deduplicate: bool


class _Parser(argparse.ArgumentParser):
    def error(self, message: str) -> Never:
        raise ValueError(message)


def _argument_parser() -> _Parser:
    parser = _Parser(prog="minimizer-stream", description=__doc__)
    parser.add_argument("input", help="FASTA/FASTQ path, or - for standard input")
    parser.add_argument("-k", type=int, required=True, help="k-mer length")
    parser.add_argument("-w", type=int, required=True, help="window size in k-mers")
    parser.add_argument("--format", choices=("tsv", "jsonl"), default="tsv")
    parser.add_argument(
        "--deduplicate",
        action="store_true",
        help="omit consecutive selections of the same minimizer",
    )
    return parser


def _tsv(record: str, item: Minimizer) -> str:
    return f"{record}\t{item.position}\t{item.hash}\t{item.strand}\t{item.kmer}"


def _jsonl(record: str, item: Minimizer) -> str:
    return json.dumps(
        {
            "record": record,
            "position": item.position,
            "hash": item.hash,
            "strand": item.strand,
            "kmer": item.kmer,
        },
        separators=(",", ":"),
    )


@contextmanager
def _input_stream(path: str) -> Iterator[TextIO]:
    if path == "-":
        yield sys.stdin
        return
    with Path(path).open(encoding="ascii") as stream:
        yield stream


def _emit(lines: Iterable[str], output: TextIO) -> None:
    for line in lines:
        output.write(line)
        output.write("\n")


def main(argv: list[str] | None = None) -> int:
    """Run the command and return a process exit status."""

    try:
        args = cast(_Arguments, _argument_parser().parse_args(argv))
        if args.k <= 0:
            raise ValueError("k must be greater than zero")
        if args.w <= 0:
            raise ValueError("w must be greater than zero")
        formatters: dict[OutputFormat, Formatter] = {"tsv": _tsv, "jsonl": _jsonl}
        formatter = formatters[args.format]
        with _input_stream(args.input) as stream:
            _emit(
                (
                    formatter(record.name, item)
                    for record in parse_records(stream)
                    for item in minimizers(
                        record.sequence,
                        args.k,
                        args.w,
                        deduplicate=args.deduplicate,
                    )
                ),
                sys.stdout,
            )
    except (OSError, UnicodeError, ValueError) as error:
        print(f"minimizer-stream: error: {error}", file=sys.stderr)
        return 2
    return 0
