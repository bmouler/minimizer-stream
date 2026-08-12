"""Deterministic comparison of deque and naive minimizer implementations."""

from __future__ import annotations

import argparse
import random
import statistics
import time

from minimizer_stream import minimizers, naive_minimizers


def _measure(function: object, sequence: str, k: int, w: int, repeats: int) -> float:
    samples = []
    for _ in range(repeats):
        started = time.perf_counter()
        list(function(sequence, k, w))  # type: ignore[operator]
        samples.append(time.perf_counter() - started)
    return statistics.median(samples)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--length", type=int, default=250_000)
    parser.add_argument("--k", type=int, default=21)
    parser.add_argument("--w", type=int, default=50)
    parser.add_argument("--repeats", type=int, default=5)
    args = parser.parse_args()
    sequence = "".join(random.Random(2026).choices("ACGT", k=args.length))
    fast_result = list(minimizers(sequence, args.k, args.w))
    slow_result = list(naive_minimizers(sequence, args.k, args.w))
    if fast_result != slow_result:
        raise RuntimeError("implementations disagree")
    fast = _measure(minimizers, sequence, args.k, args.w, args.repeats)
    slow = _measure(naive_minimizers, sequence, args.k, args.w, args.repeats)
    print(f"outputs={len(fast_result)} identical=true")
    print(f"deque_seconds={fast:.6f} naive_seconds={slow:.6f} speedup={slow / fast:.2f}x")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
