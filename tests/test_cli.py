from __future__ import annotations

import json
import runpy
import subprocess
import sys
from pathlib import Path

import pytest

from minimizer_stream import minimizers
from minimizer_stream.cli import main


def test_module_cli_tsv_from_clean_directory(tmp_path: Path) -> None:
    input_path = tmp_path / "reads.fa"
    input_path.write_text(">alpha\nACGTA\n>beta\nAAANAAAA\n", encoding="ascii")
    command = Path(sys.executable).with_name("minimizer-stream")
    result = subprocess.run(
        [str(command), str(input_path), "-k", "2", "-w", "2"],
        cwd=tmp_path,
        capture_output=True,
        text=True,
        check=False,
    )
    expected = [
        f"{name}\t{item.position}\t{item.hash}\t{item.strand}\t{item.kmer}"
        for name, sequence in [("alpha", "ACGTA"), ("beta", "AAANAAAA")]
        for item in minimizers(sequence, 2, 2)
    ]
    assert result.returncode == 0
    assert result.stdout.splitlines() == expected
    assert result.stderr == ""


def test_module_cli_jsonl_fastq_stdin_and_dedup(tmp_path: Path) -> None:
    command = Path(sys.executable).with_name("minimizer-stream")
    result = subprocess.run(
        [
            str(command),
            "-",
            "-k",
            "2",
            "-w",
            "2",
            "--format",
            "jsonl",
            "--deduplicate",
        ],
        cwd=tmp_path,
        input="@read\nACGTA\n+\n!!!!!\n",
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 0
    payloads = [json.loads(line) for line in result.stdout.splitlines()]
    assert payloads == [
        {
            "record": "read",
            "position": item.position,
            "hash": item.hash,
            "strand": item.strand,
            "kmer": item.kmer,
        }
        for item in minimizers("ACGTA", 2, 2, deduplicate=True)
    ]


def test_main_success_paths(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    input_path = tmp_path / "read.fa"
    input_path.write_text(">read\nACGTA\n", encoding="ascii")
    assert main([str(input_path), "-k", "2", "-w", "2"]) == 0
    assert capsys.readouterr().out.count("\n") == 3
    assert main([str(input_path), "-k", "2", "-w", "2", "--format", "jsonl"]) == 0
    jsonl_output = capsys.readouterr().out.splitlines()
    assert json.loads(jsonl_output[0])["record"] == "read"
    assert jsonl_output[0] == (
        '{"record":"read","position":0,"hash":650881229941436253,"strand":"+","kmer":"AC"}'
    )


def test_main_passes_deduplication_to_stream(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    input_path = tmp_path / "read.fa"
    input_path.write_text(">read\nACGTA\n", encoding="ascii")
    assert main([str(input_path), "-k", "2", "-w", "2", "--format", "jsonl", "--deduplicate"]) == 0
    assert len(capsys.readouterr().out.splitlines()) == 2


def test_help_documents_the_cli_contract(capsys: pytest.CaptureFixture[str]) -> None:
    with pytest.raises(SystemExit, match="0"):
        main(["--help"])
    assert capsys.readouterr().out == (
        "usage: minimizer-stream [-h] -k K -w W [--format {tsv,jsonl}] [--deduplicate]\n"
        "                        input\n"
        "\n"
        "Command-line interface for minimizer-stream.\n"
        "\n"
        "positional arguments:\n"
        "  input                 FASTA/FASTQ path, or - for standard input\n"
        "\n"
        "options:\n"
        "  -h, --help            show this help message and exit\n"
        "  -k K                  k-mer length\n"
        "  -w W                  window size in k-mers\n"
        "  --format {tsv,jsonl}\n"
        "  --deduplicate         omit consecutive selections of the same minimizer\n"
    )


def test_module_entrypoint_in_process(
    monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    monkeypatch.setattr(sys, "argv", ["minimizer-stream", "-", "-k", "2", "-w", "2"])
    monkeypatch.setattr(sys, "stdin", iter([">read\n", "ACGT\n"]))
    with pytest.raises(SystemExit, match="0"):
        runpy.run_module("minimizer_stream.__main__", run_name="__main__")
    assert capsys.readouterr().out


@pytest.mark.parametrize(
    ("arguments", "message"),
    [
        (["-", "-k", "0", "-w", "1"], "k must be greater than zero"),
        (["-", "-k", "1", "-w", "0"], "w must be greater than zero"),
        (["-", "-k", "1"], "required"),
        (["-", "-w", "1"], "required"),
        (["-", "-k", "1", "-w", "1", "--format", "xml"], "invalid choice"),
    ],
)
def test_main_reports_argument_errors(
    arguments: list[str], message: str, capsys: pytest.CaptureFixture[str]
) -> None:
    assert main(arguments) == 2
    error = capsys.readouterr().err
    if "must be greater" in message:
        assert error == f"minimizer-stream: error: {message}\n"
    else:
        assert message in error


def test_main_reports_malformed_and_missing_files(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    monkeypatch.setattr(sys, "stdin", iter(["bad\n"]))
    assert main(["-", "-k", "1", "-w", "1"]) == 2
    assert "expected a FASTA or FASTQ header" in capsys.readouterr().err
    assert main([str(tmp_path / "missing.fa"), "-k", "1", "-w", "1"]) == 2
    assert "No such file" in capsys.readouterr().err


def test_main_reports_non_ascii_input(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    input_path = tmp_path / "non-ascii.fa"
    input_path.write_bytes(b">read\nAC\xffT\n")
    assert main([str(input_path), "-k", "1", "-w", "1"]) == 2
    assert "ascii" in capsys.readouterr().err
