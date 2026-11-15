#!/usr/bin/env python3
"""Assert stable and semantic outputs from the tracked workflow suites."""

from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from pathlib import Path

import yaml
from jsonschema import Draft202012Validator

ROOT = Path(__file__).resolve().parents[1]
INFERENCE_OUTPUT_SCHEMA = yaml.safe_load(
    (ROOT / "workflow" / "schemas" / "auto_preprocess_output.schema.yaml").read_text(
        encoding="utf-8"
    )
)
Draft202012Validator.check_schema(INFERENCE_OUTPUT_SCHEMA)
INFERENCE_OUTPUT_VALIDATOR = Draft202012Validator(INFERENCE_OUTPUT_SCHEMA)
SHORT_PE_INFERENCE_PAYLOAD = {
    "recommendations": {
        "read1": {
            "headcrop_bases": 0,
            "tailcrop_base_position": 250,
            "maximum_length_after_headcrop_bases": 250,
        },
        "read2": {
            "headcrop_bases": 0,
            "tailcrop_base_position": 198,
            "maximum_length_after_headcrop_bases": 198,
        },
    },
    "tool_args": {
        "fastp": {
            "--max_len1": 250,
            "--max_len2": 198,
            "--average_qual": 30,
            "--length_required": 30,
        }
    },
    "warnings": [],
}
SHORT_PE_INFERENCE_LOG = "\n".join(
    [
        "Auto-preprocess inference started",
        "Auto-preprocess inference complete",
        (
            "Recommendations: read1(headcrop_bases=0, tailcrop_base_position=250, "
            "maximum_length_after_headcrop_bases=250); read2(headcrop_bases=0, "
            "tailcrop_base_position=198, maximum_length_after_headcrop_bases=198)"
        ),
        (
            "Tool arguments: fastp --max_len1 250 --max_len2 198 "
            "--average_qual 30 --length_required 30"
        ),
        "Warnings: none",
        "",
    ]
)
HYBRID_SHORT_INFERENCE_PAYLOAD = {
    "recommendations": {
        "read1": {
            "headcrop_bases": 0,
            "tailcrop_base_position": 30,
            "maximum_length_after_headcrop_bases": 30,
        },
        "read2": {
            "headcrop_bases": 0,
            "tailcrop_base_position": 30,
            "maximum_length_after_headcrop_bases": 30,
        },
    },
    "tool_args": {
        "fastp": {
            "--max_len1": 30,
            "--max_len2": 30,
            "--average_qual": 30,
            "--length_required": 30,
        }
    },
    "warnings": [
        {
            "code": "tailcrop_safeguard_adjusted",
            "message": (
                "Retained-length safeguard adjusted a raw tailcrop cutoff for "
                "read1 and read2."
            ),
        },
    ],
}
HYBRID_SHORT_INFERENCE_LOG = "\n".join(
    [
        "Auto-preprocess inference started",
        "Auto-preprocess inference complete",
        (
            "Recommendations: read1(headcrop_bases=0, tailcrop_base_position=30, "
            "maximum_length_after_headcrop_bases=30); read2(headcrop_bases=0, "
            "tailcrop_base_position=30, maximum_length_after_headcrop_bases=30)"
        ),
        (
            "Tool arguments: fastp --max_len1 30 --max_len2 30 "
            "--average_qual 30 --length_required 30"
        ),
        (
            "Warnings: tailcrop_safeguard_adjusted: Retained-length safeguard "
            "adjusted a raw tailcrop cutoff for read1 and read2."
        ),
        "",
    ]
)


@dataclass(frozen=True)
class InferenceExpectation:
    readset: str
    payload: dict
    log: str
    absent_readsets: tuple[str, ...] = ()


@dataclass(frozen=True)
class FastqContentExpectation:
    output_relative: str
    identifiers: tuple[str, ...]
    first_sequence_prefix: str = ""


@dataclass(frozen=True)
class ArtifactInventory:
    directories: frozenset[str]
    files: frozenset[str]
    forbidden: frozenset[str]


@dataclass(frozen=True)
class Case:
    sample: str
    readsets: tuple[str, ...]
    read_comparisons: tuple[tuple[str, str], ...]
    short_readsets: tuple[str, ...] = ()
    long_readsets: tuple[str, ...] = ()
    fastq_contents: tuple[FastqContentExpectation, ...] = ()
    porechop_routed_readsets: tuple[str, ...] = ()
    porechop_bypassed_readsets: tuple[str, ...] = ()
    inference: InferenceExpectation | None = None


CASES = {
    "short-se": Case(
        sample="test_short_se",
        readsets=("test_short_se_short",),
        read_comparisons=(
            (
                "reads/test_short_se_short.fastq",
                "expected_short_se.fastq",
            ),
        ),
        short_readsets=("test_short_se_short",),
    ),
    "short-pe-disabled": Case(
        sample="test_short_pe",
        readsets=("test_short_pe_short",),
        read_comparisons=(
            (
                "reads/test_short_pe_short_R1.fastq",
                "expected_short_pe_R1.fastq",
            ),
            (
                "reads/test_short_pe_short_R2.fastq",
                "expected_short_pe_R2.fastq",
            ),
        ),
        short_readsets=("test_short_pe_short",),
    ),
    "short-pe-infer": Case(
        sample="test_infer_pe",
        readsets=("test_infer_pe_short",),
        read_comparisons=(
            (
                "reads/test_infer_pe_short_R1.fastq",
                "expected_short_infer_pe_R1.fastq",
            ),
            (
                "reads/test_infer_pe_short_R2.fastq",
                "expected_short_infer_pe_R2.fastq",
            ),
        ),
        short_readsets=("test_infer_pe_short",),
        inference=InferenceExpectation(
            readset="test_infer_pe_short",
            payload=SHORT_PE_INFERENCE_PAYLOAD,
            log=SHORT_PE_INFERENCE_LOG,
        ),
    ),
    "long": Case(
        sample="test_long",
        readsets=("test_long_long",),
        read_comparisons=(
            (
                "reads/test_long_long.fastq",
                "expected_long.fastq",
            ),
        ),
        long_readsets=("test_long_long",),
        porechop_routed_readsets=("test_long_long",),
    ),
    "nano-raw": Case(
        sample="test_nano_raw",
        readsets=("test_nano_raw_long",),
        read_comparisons=(
            (
                "reads/test_nano_raw_long.fastq",
                "expected_long.fastq",
            ),
        ),
        long_readsets=("test_nano_raw_long",),
        porechop_routed_readsets=("test_nano_raw_long",),
    ),
    "pacbio-raw": Case(
        sample="test_pacbio_raw",
        readsets=("test_pacbio_raw_long",),
        read_comparisons=(
            (
                "reads/test_pacbio_raw_long.fastq",
                "expected_pacbio_raw.fastq",
            ),
        ),
        long_readsets=("test_pacbio_raw_long",),
        fastq_contents=(
            FastqContentExpectation(
                output_relative="reads/test_pacbio_raw_long.fastq",
                identifiers=("read_keep",),
                first_sequence_prefix="AATGTACTTCGTTCAGTTACGTATTGCT",
            ),
        ),
        porechop_bypassed_readsets=("test_pacbio_raw_long",),
    ),
    "hybrid": Case(
        sample="test_hybrid",
        readsets=("test_hybrid_short", "test_hybrid_long"),
        read_comparisons=(
            (
                "reads/test_hybrid_short_R1.fastq",
                "expected_hybrid_infer_short_R1.fastq",
            ),
            (
                "reads/test_hybrid_short_R2.fastq",
                "expected_hybrid_infer_short_R2.fastq",
            ),
            (
                "reads/test_hybrid_long.fastq",
                "expected_hybrid_infer_long.fastq",
            ),
        ),
        short_readsets=("test_hybrid_short",),
        long_readsets=("test_hybrid_long",),
        porechop_routed_readsets=("test_hybrid_long",),
        inference=InferenceExpectation(
            readset="test_hybrid_short",
            payload=HYBRID_SHORT_INFERENCE_PAYLOAD,
            log=HYBRID_SHORT_INFERENCE_LOG,
            absent_readsets=("test_hybrid_long",),
        ),
    ),
}

SUITES = {
    "inference": ("short-pe-infer", "hybrid"),
    "no-inference": (
        "short-se",
        "short-pe-disabled",
        "long",
        "nano-raw",
        "pacbio-raw",
    ),
}


def build_artifact_inventory(case: Case) -> ArtifactInventory:
    """Build one explicit test-owned inventory for a workflow case."""
    files = {output_relative for output_relative, _ in case.read_comparisons}
    files.update(
        f"qc/{stage}/{readset}.sequali.{extension}"
        for readset in case.readsets
        for stage in ("raw", "processed")
        for extension in ("html", "json")
    )
    files.update(
        f"logs/sequali/{stage}/{readset}.sequali.log"
        for readset in case.readsets
        for stage in ("raw", "processed")
    )
    files.update(
        f"qc/processed/{readset}.fastp.{extension}"
        for readset in case.short_readsets
        for extension in ("html", "json")
    )
    files.update(f"logs/fastp/{readset}.fastp.log" for readset in case.short_readsets)
    files.update(
        f"logs/chopper/{readset}.chopper.log" for readset in case.long_readsets
    )
    files.update(
        f"logs/porechop_abi/{readset}.porechop_abi.log"
        for readset in case.porechop_routed_readsets
    )

    forbidden = {
        f"reads/{readset}.porechop.fastq"
        for readset in (
            *case.porechop_routed_readsets,
            *case.porechop_bypassed_readsets,
        )
    }
    forbidden.update(
        f"logs/porechop_abi/{readset}.porechop_abi.log"
        for readset in case.porechop_bypassed_readsets
    )
    forbidden.update(
        f"qc/{stage}/{readset}.sequali.log"
        for readset in case.readsets
        for stage in ("raw", "processed")
    )
    forbidden.update(
        f"qc/processed/{readset}.{tool}.log"
        for readset, tool in (
            *((readset, "fastp") for readset in case.short_readsets),
            *((readset, "chopper") for readset in case.long_readsets),
            *(
                (readset, "porechop_abi")
                for readset in (
                    *case.porechop_routed_readsets,
                    *case.porechop_bypassed_readsets,
                )
            ),
        )
    )

    if case.inference is None:
        absent_inference_readsets = case.readsets
    else:
        files.update(
            {
                f"reads/{case.inference.readset}.auto_params.json",
                f"logs/auto_preprocess/{case.inference.readset}.auto_params.log",
            }
        )
        absent_inference_readsets = case.inference.absent_readsets
    forbidden.update(
        f"reads/{readset}.auto_params.json" for readset in absent_inference_readsets
    )
    forbidden.update(
        f"logs/auto_preprocess/{readset}.auto_params.log"
        for readset in absent_inference_readsets
    )
    forbidden.update(f"reads/{readset}.auto_params.log" for readset in case.readsets)

    directories = {"reads", "qc", "qc/raw", "qc/processed", "logs"}
    directories.update({"logs/sequali", "logs/sequali/raw", "logs/sequali/processed"})
    if case.short_readsets:
        directories.add("logs/fastp")
    if case.long_readsets:
        directories.add("logs/chopper")
    if case.porechop_routed_readsets:
        directories.add("logs/porechop_abi")
    if case.inference is not None:
        directories.add("logs/auto_preprocess")

    return ArtifactInventory(
        directories=frozenset(directories),
        files=frozenset(files),
        forbidden=frozenset(forbidden),
    )


CASE_INVENTORIES = {
    case_id: build_artifact_inventory(case) for case_id, case in CASES.items()
}


def require_nonempty(path: Path) -> None:
    if not path.is_file() or path.stat().st_size == 0:
        raise AssertionError(f"Required nonempty file is missing or empty: {path}")


def assert_reads(case: Case, workdir: Path, result_root: Path) -> None:
    expected_root = workdir / "expected"
    for output_relative, expected_name in case.read_comparisons:
        output = result_root / output_relative
        expected = expected_root / expected_name
        require_nonempty(expected)
        if output.read_bytes() != expected.read_bytes():
            raise AssertionError(
                f"Processed reads differ from tracked fixture: {output} != {expected}"
            )


def assert_fastq_contents(case: Case, result_root: Path) -> None:
    for expectation in case.fastq_contents:
        output = result_root / expectation.output_relative
        lines = output.read_text(encoding="utf-8").splitlines()
        if not lines or len(lines) % 4:
            raise AssertionError(f"Malformed FASTQ record structure: {output}")

        records = []
        for offset in range(0, len(lines), 4):
            header, sequence, separator, quality = lines[offset : offset + 4]
            if not header.startswith("@") or not separator.startswith("+"):
                raise AssertionError(f"Malformed FASTQ record structure: {output}")
            if len(sequence) != len(quality):
                raise AssertionError(
                    f"FASTQ sequence/quality length mismatch: {output}"
                )
            records.append((header[1:].split(maxsplit=1)[0], sequence))

        identifiers = tuple(identifier for identifier, _ in records)
        if identifiers != expectation.identifiers:
            raise AssertionError(
                f"Unexpected FASTQ identifiers: {identifiers!r} "
                f"!= {expectation.identifiers!r}"
            )
        if expectation.first_sequence_prefix and not records[0][1].startswith(
            expectation.first_sequence_prefix
        ):
            raise AssertionError(
                f"First FASTQ sequence lacks required prefix "
                f"{expectation.first_sequence_prefix!r}: {output}"
            )


def assert_inference(case: Case, result_root: Path) -> None:
    if case.inference is None:
        return

    expectation = case.inference
    json_path = result_root / "reads" / f"{expectation.readset}.auto_params.json"
    log_path = (
        result_root / "logs/auto_preprocess" / f"{expectation.readset}.auto_params.log"
    )
    payload = json.loads(json_path.read_text(encoding="utf-8"))
    INFERENCE_OUTPUT_VALIDATOR.validate(payload)
    if payload != expectation.payload:
        raise AssertionError(f"Unexpected inference payload: {payload!r}")
    if log_path.read_text(encoding="utf-8") != expectation.log:
        raise AssertionError(f"Unexpected inference log content: {log_path}")


def _tree_entry_types(result_root: Path) -> dict[str, str]:
    entries = {}
    pending = [result_root]
    while pending:
        directory = pending.pop()
        for path in directory.iterdir():
            relative = path.relative_to(result_root).as_posix()
            if path.is_symlink():
                entries[relative] = "symlink"
            elif path.is_dir():
                entries[relative] = "directory"
                pending.append(path)
            elif path.is_file():
                entries[relative] = "file"
            else:
                entries[relative] = "other"
    return entries


def assert_artifact_inventory(
    case: Case, inventory: ArtifactInventory, result_root: Path
) -> None:
    if result_root.is_symlink() or not result_root.is_dir():
        raise AssertionError(
            f"Required sample result directory is missing: {result_root}"
        )

    for relative in inventory.forbidden:
        path = result_root / relative
        if path.exists() or path.is_symlink():
            raise AssertionError(f"Forbidden workflow artifact exists: {path}")

    expected_types = {
        **{relative: "directory" for relative in inventory.directories},
        **{relative: "file" for relative in inventory.files},
    }
    actual_types = _tree_entry_types(result_root)
    if set(actual_types) != set(expected_types):
        missing = sorted(set(expected_types) - set(actual_types))
        extra = sorted(set(actual_types) - set(expected_types))
        raise AssertionError(
            f"Unexpected result tree for {case.sample}: missing={missing!r}, "
            f"extra={extra!r}"
        )

    wrong_types = {
        relative: (expected_types[relative], actual_types[relative])
        for relative in expected_types
        if actual_types[relative] != expected_types[relative]
    }
    if wrong_types:
        raise AssertionError(
            f"Unexpected artifact types for {case.sample}: {wrong_types!r}"
        )

    empty_files = sorted(
        relative
        for relative in inventory.files
        if (result_root / relative).stat().st_size == 0
    )
    if empty_files:
        raise AssertionError(
            f"Expected nonempty artifacts are empty for {case.sample}: {empty_files!r}"
        )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--suite", required=True, choices=sorted(SUITES))
    parser.add_argument("--workdir", type=Path, default=Path(".test"))
    args = parser.parse_args()

    workdir = args.workdir.resolve()
    for case_id in SUITES[args.suite]:
        case = CASES[case_id]
        result_root = workdir / "results" / case.sample
        assert_artifact_inventory(case, CASE_INVENTORIES[case_id], result_root)
        assert_reads(case, workdir, result_root)
        assert_fastq_contents(case, result_root)
        assert_inference(case, result_root)
    print(f"Workflow outputs validated for {args.suite} suite")


if __name__ == "__main__":
    main()
