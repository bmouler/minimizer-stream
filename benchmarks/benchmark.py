"""Deterministic end-to-end parser and minimizer benchmark."""

from __future__ import annotations

import argparse
import hashlib
import json
import random
import statistics
import time
from collections.abc import Callable, Iterable

from minimizer_stream import Minimizer, minimizers, naive_minimizers, parse_records

Output = list[tuple[str, Minimizer]]
DEFAULT_WORKLOAD = (8, 8_000, 21, 100)
DEFAULT_EXPECTED_CHECKSUM = "faaba963089826ffa15afd9f9d9e0d63e3dd5d6de96ea6a3632743fb000fb008"


def _fixtures(records: int, length: int) -> tuple[tuple[str, ...], tuple[str, ...]]:
    generator = random.Random(20260815)
    fasta: list[str] = []
    fastq: list[str] = []
    for index in range(records):
        sequence = "".join(
            generator.choices(
                "ACGTNacgt",
                weights=(100, 100, 100, 100, 1, 20, 20, 20, 20),
                k=length,
            )
        )
        fasta.append(f">fasta-{index} representative record\n")
        for offset in range(0, length, 80):
            fasta.append(sequence[offset : offset + 80] + "\n")
        fastq.extend(
            (
                f"@fastq-{index} representative record\n",
                sequence + "\n",
                "+\n",
                "I" * length + "\n",
            )
        )
    return tuple(fasta), tuple(fastq)


def _workload(
    fixtures: tuple[tuple[str, ...], tuple[str, ...]],
    function: Callable[..., Iterable[Minimizer]],
    k: int,
    w: int,
) -> Output:
    output: Output = []
    append = output.append
    for lines in fixtures:
        for record in parse_records(iter(lines)):
            name = record.name
            for item in function(record.sequence, k, w, deduplicate=True):
                append((name, item))
    return output


def _checksum(output: Output) -> str:
    digest = hashlib.sha256()
    for name, item in output:
        digest.update(name.encode("ascii"))
        digest.update(b"\0")
        digest.update(item.position.to_bytes(8, "little", signed=True))
        digest.update(item.hash.to_bytes(8, "little"))
        digest.update(item.strand.encode("ascii"))
        digest.update(item.kmer.encode("ascii"))
        digest.update(b"\n")
    return digest.hexdigest()


def _measure(
    fixtures: tuple[tuple[str, ...], tuple[str, ...]],
    k: int,
    w: int,
    warmups: int,
    samples: int,
    expected: Output,
) -> list[float]:
    for _ in range(warmups):
        if _workload(fixtures, minimizers, k, w) != expected:
            raise RuntimeError("minimizer result changed during warmup")
    timings = []
    for _ in range(samples):
        started = time.perf_counter()
        measured = _workload(fixtures, minimizers, k, w)
        timings.append(time.perf_counter() - started)
        if measured != expected:
            raise RuntimeError("minimizer result changed between timed samples")
    return timings


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--records", type=int, default=8)
    parser.add_argument("--length", type=int, default=8_000)
    parser.add_argument("--k", type=int, default=21)
    parser.add_argument("--w", type=int, default=100)
    parser.add_argument("--warmups", type=int, default=2)
    parser.add_argument("--samples", type=int, default=11)
    parser.add_argument("--json", action="store_true")
    args = parser.parse_args()
    fixtures = _fixtures(args.records, args.length)
    result = _workload(fixtures, minimizers, args.k, args.w)
    oracle = _workload(fixtures, naive_minimizers, args.k, args.w)
    if result != oracle:
        raise RuntimeError("streaming and reference implementations disagree")
    checksum = _checksum(result)
    if (
        args.records,
        args.length,
        args.k,
        args.w,
    ) == DEFAULT_WORKLOAD and checksum != DEFAULT_EXPECTED_CHECKSUM:
        raise RuntimeError(
            f"expected default result digest {DEFAULT_EXPECTED_CHECKSUM}, got {checksum}"
        )
    timings = _measure(fixtures, args.k, args.w, args.warmups, args.samples, result)
    report = {
        "records": args.records * 2,
        "bases": args.records * args.length * 2,
        "outputs": len(result),
        "identical": True,
        "checksum": checksum,
        "samples": args.samples,
        "median_seconds": statistics.median(timings),
        "min_seconds": min(timings),
        "max_seconds": max(timings),
    }
    if args.json:
        print(json.dumps(report, sort_keys=True))
    else:
        print(
            f"records={report['records']} bases={report['bases']} "
            f"outputs={report['outputs']} identical=true checksum={report['checksum']}"
        )
        print(
            f"samples={report['samples']} median_seconds={report['median_seconds']:.6f} "
            f"min_seconds={report['min_seconds']:.6f} "
            f"max_seconds={report['max_seconds']:.6f}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
