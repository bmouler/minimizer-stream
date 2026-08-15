# Changelog

## [Unreleased]

- Added a packed rolling two-bit fast path for exact-string `k=21` minimizers and accelerated bit-identical FNV-1a hashing.
- Added a deterministic end-to-end FASTA/FASTQ parse and minimizer benchmark with exact oracle equality.


## [1.0.0] - 2026-08-12

First stable release.

- Streaming canonical minimizer sketches for FASTA and FASTQ.
- Added a deterministic Hypothesis property-based suite alongside differential-oracle tests.
- Documented a 98.98% mutation score: 388/392 mutants killed, with four behavior-equivalent survivors.
- Adopted strict mypy typing.
- Expanded CI to Linux and macOS on Python 3.11–3.13.
