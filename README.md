# minimizer-stream

[![CI](https://github.com/bmouler/minimizer-stream/actions/workflows/ci.yml/badge.svg)](https://github.com/bmouler/minimizer-stream/actions/workflows/ci.yml) [![branch coverage](https://img.shields.io/badge/branch%20coverage-100%25-brightgreen)](https://github.com/bmouler/minimizer-stream/actions/workflows/ci.yml) [![Python 3.11+](https://img.shields.io/badge/python-3.11%2B-blue)](https://www.python.org/) [![MIT](https://img.shields.io/badge/license-MIT-blue)](LICENSE)

A dependency-free Python library and command-line program that streams canonical minimizer sketches from FASTA and four-line FASTQ. Output includes record, zero-based position, stable 64-bit FNV-1a hash, canonical strand, and canonical k-mer.

## Install

```bash
python -m pip install .
```

For development and verification:

```bash
python -m pip install -e '.[test]'
ruff check .
pytest --cov=minimizer_stream --cov-branch --cov-fail-under=100
```

## Quickstart

```bash
minimizer-stream reads.fasta -k 21 -w 10 > sketch.tsv
minimizer-stream reads.fastq -k 21 -w 10 --format jsonl --deduplicate
cat reads.fasta | minimizer-stream - -k 15 -w 5
```

TSV has no header and its columns are `record`, `position`, `hash`, `strand`, `kmer`. JSONL uses the same field names. The Python API is lazy:

```python
from minimizer_stream import minimizers, naive_minimizers

for item in minimizers("ACGTTGCA", k=3, w=2, deduplicate=True):
    print(item.position, item.hash, item.strand, item.kmer)

assert list(minimizers("ACGTTGCA", 3, 2)) == list(naive_minimizers("ACGTTGCA", 3, 2))
```

## Algorithm

Each k-mer is compared lexicographically with its reverse complement; the smaller representation is canonical, with palindromes assigned to the forward strand. Canonical strings are hashed with deterministic 64-bit FNV-1a, never Python's process-randomized `hash`. A monotonic deque retains possible minima for each window, removing expired entries from the front and strictly larger hashes from the back. Keeping equal hashes preserves deterministic leftmost tie-breaking.

Constructing, reverse-complementing, and hashing each canonical k-mer costs $O(k)$ in this pure-Python implementation. The monotonic-deque path is therefore $O(nk)$ total with $O(w+k)$ working memory, versus $O(n(k+w))$ for `naive_minimizers`, which additionally scans each minimizer window. Non-ACGT characters reset both k-mer construction and the minimizer window, so no output spans an invalid base.

## Reproducible capability evidence

The deterministic benchmark seeds its input and first requires both implementations to produce exactly identical objects:

```bash
python benchmarks/benchmark.py --length 50000 --k 21 --w 50 --repeats 3
```

On an Apple M3 Max with Python 3.11.12, that command produced 49,931 identical selections and median times of 0.130359 seconds for the deque implementation and 0.195500 seconds for the oracle, a 1.50x speedup. Timing depends on hardware and Python version; the durable claim is removing the extra $O(nw)$ window scan through amortized $O(n)$ deque maintenance, while shared k-mer construction remains $O(nk)$. The exhaustive deterministic randomized test also compares both APIs across lengths 0–39 and all combinations of `k` and `w` from 1–5.

## Verification

CI runs Ruff and the full test suite on Python 3.11 and 3.12. Its exact coverage command is:

```bash
pytest --cov=minimizer_stream --cov-branch --cov-fail-under=100
```

Tests include reverse complements, palindromes, ties, invalid-base resets, short and lowercase sequences, FASTA, FASTQ, malformed records, invalid parameters, subprocess execution from a clean temporary directory, TSV, JSONL, standard input, and deduplication.

## Limitations

- FASTQ must use the conventional four-line form; wrapped sequence or quality lines are rejected.
- FASTA permits wrapped sequence lines but rejects blank lines and empty records.
- Record names containing control characters or beginning with `=`, `+`, `-`, or `@` are rejected so default TSV output keeps five columns and is safe from spreadsheet formula interpretation.
- Input is decoded as ASCII; compressed files must be decompressed before use.
- Record sequences are retained one at a time. The parser does not load the full file, but an individual record must fit in memory.
- `w` counts consecutive valid k-mers, not bases. A valid window therefore spans `k + w - 1` bases.
- Hash collisions are possible for any 64-bit sketch; the canonical k-mer is emitted so callers can distinguish them when needed.
