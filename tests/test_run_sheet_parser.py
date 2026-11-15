import importlib.util
from pathlib import Path

import pytest
from jsonschema import ValidationError

ROOT = Path(__file__).resolve().parents[1]
RUN_SHEET_SOURCE = ROOT / "workflow" / "scripts" / "common" / "run_sheet.py"
RUN_SHEET_SPEC = importlib.util.spec_from_file_location(
    "test_run_sheet_helpers", RUN_SHEET_SOURCE
)
assert RUN_SHEET_SPEC is not None and RUN_SHEET_SPEC.loader is not None
RUN_SHEET_HELPERS = importlib.util.module_from_spec(RUN_SHEET_SPEC)
RUN_SHEET_SPEC.loader.exec_module(RUN_SHEET_HELPERS)
RunSheetError = RUN_SHEET_HELPERS.RunSheetError
build_samples_row_validator = RUN_SHEET_HELPERS.build_samples_row_validator
load_run_sheet = RUN_SHEET_HELPERS.load_run_sheet
SAMPLES_VALIDATOR = build_samples_row_validator(
    ROOT / "workflow" / "schemas" / "samples.schema.yaml"
)
COMPLETE_HEADER = ["group", "sample", "data_type", "read1"]
COMPLETE_ROW = ["group1", "sample1", "short_se", "inputs/reads.fastq"]


def _write_sheet(tmp_path, header=COMPLETE_HEADER, rows=(COMPLETE_ROW,)):
    path = tmp_path / "run.tsv"
    path.write_text(
        "\t".join(header) + "\n" + "".join("\t".join(row) + "\n" for row in rows),
        encoding="utf-8",
    )
    return path


def _load(path, validator=SAMPLES_VALIDATOR):
    validate_row = None if validator is None else validator.validate
    return load_run_sheet(path, validate_row)


@pytest.mark.parametrize("missing_header", ["group", "sample", "data_type", "read1"])
def test_each_required_header_is_mandatory(tmp_path, missing_header):
    index = COMPLETE_HEADER.index(missing_header)
    header = COMPLETE_HEADER[:index] + COMPLETE_HEADER[index + 1 :]
    row = COMPLETE_ROW[:index] + COMPLETE_ROW[index + 1 :]

    with pytest.raises(
        RunSheetError,
        match=rf"Run sheet is missing required header\(s\): {missing_header}",
    ):
        _load(_write_sheet(tmp_path, header, [row]))


@pytest.mark.parametrize(
    ("header", "message"),
    [
        (["group", "", "data_type", "read1"], "non-empty TSV header"),
        (
            ["group", "sample", "sample", "data_type", "read1"],
            r"duplicate header\(s\): sample",
        ),
    ],
)
def test_invalid_headers_are_rejected(tmp_path, header, message):
    row = ["value"] * len(header)

    with pytest.raises(RunSheetError, match=message):
        _load(_write_sheet(tmp_path, header, [row]))


def test_header_and_values_are_trimmed_before_validation_and_normalization(
    tmp_path, monkeypatch
):
    monkeypatch.setenv("READ_ROOT", "/expanded/reads")
    monkeypatch.setenv("CONTAM_ROOT", "/expanded/contaminants")
    monkeypatch.setenv("HOME", "/home/tester")
    path = _write_sheet(
        tmp_path,
        [" group ", " sample ", " data_type ", " read1 ", " read2 ", " contam "],
        [
            [
                " group1 ",
                " paired ",
                " short_pe ",
                " relative/R1.fastq ",
                " $READ_ROOT/R2.fastq ",
                " ",
            ]
        ],
    )

    assert _load(path) == [
        {
            "group": "group1",
            "sample": "paired",
            "data_type": "short_pe",
            "read1": "relative/R1.fastq",
            "read2": "/expanded/reads/R2.fastq",
            "contam": "",
        }
    ]
    assert (
        _load(
            _write_sheet(
                tmp_path,
                COMPLETE_HEADER,
                [["group1", "sample1", "short_se", "~/reads.fastq"]],
            )
        )[0]["read1"]
        == "/home/tester/reads.fastq"
    )
    assert (
        _load(
            _write_sheet(
                tmp_path,
                [*COMPLETE_HEADER, "contam"],
                [
                    [
                        "group1",
                        "long",
                        "nano-raw",
                        "relative/long.fastq",
                        "$CONTAM_ROOT/contam.fa",
                    ]
                ],
            )
        )[0]["contam"]
        == "/expanded/contaminants/contam.fa"
    )


def test_missing_run_sheet_is_rejected(tmp_path):
    path = tmp_path / "missing.tsv"

    with pytest.raises(RunSheetError, match=f"Run sheet does not exist: {path}"):
        _load(path)


def test_empty_run_sheet_is_rejected(tmp_path):
    with pytest.raises(RunSheetError, match="must contain at least one data row"):
        _load(_write_sheet(tmp_path, rows=[]))


@pytest.mark.parametrize(
    ("row", "message"),
    [
        (COMPLETE_ROW + ["extra"], "has more fields than the header"),
        (COMPLETE_ROW[:-1], "has fewer fields than the header"),
    ],
)
def test_row_field_count_must_match_header(tmp_path, row, message):
    with pytest.raises(RunSheetError, match=message):
        _load(_write_sheet(tmp_path, rows=[row]))


def test_malformed_tsv_is_rejected(tmp_path):
    path = tmp_path / "run.tsv"
    path.write_text(
        "group\tsample\tdata_type\tread1\n"
        'group1\tsample1\tshort_se\t"unterminated\n',
        encoding="utf-8",
    )

    with pytest.raises(RunSheetError, match="Malformed TSV run sheet"):
        _load(path)


def test_samples_schema_rejection_is_propagated_directly(tmp_path):
    with pytest.raises(ValidationError):
        _load(
            _write_sheet(
                tmp_path,
                rows=[["", "sample1", "short_se", "inputs/reads.fastq"]],
            )
        )


def test_explicit_names_are_reserved_and_counters_are_per_group_and_type(tmp_path):
    header = ["group", "sample", "data_type", "read1", "read2"]
    rows = [
        ["g", "", "short_se", "inputs/one.fastq", ""],
        ["g", "", "short_pe", "inputs/two_R1.fastq", "inputs/two_R2.fastq"],
        ["h", "", "short_se", "inputs/three.fastq", ""],
        ["g", "short_se_1", "short_se", "inputs/four.fastq", ""],
        ["g", "", "short_se", "inputs/five.fastq", ""],
        ["g", "short_pe_1", "short_pe", "inputs/six_R1.fastq", "inputs/six_R2.fastq"],
    ]

    parsed = _load(_write_sheet(tmp_path, header, rows))

    assert [(row["group"], row["sample"]) for row in parsed] == [
        ("g", "short_se_2"),
        ("g", "short_pe_2"),
        ("h", "short_se_1"),
        ("g", "short_se_1"),
        ("g", "short_se_3"),
        ("g", "short_pe_1"),
    ]
    assert [row["read1"] for row in parsed] == [row[3] for row in rows]


def test_duplicate_normalized_group_sample_pair_is_rejected(tmp_path):
    rows = [
        ["group1", "duplicate", "short_se", "inputs/one.fastq"],
        ["group1", "duplicate", "long", "inputs/two.fastq"],
    ]

    with pytest.raises(
        RunSheetError,
        match=r"duplicate normalized \(group, sample\) pair: \('group1', 'duplicate'\)",
    ):
        _load(_write_sheet(tmp_path, rows=rows))


def test_duplicate_sample_text_in_different_groups_is_allowed(tmp_path):
    rows = [
        ["group1", "duplicate", "short_se", "inputs/one.fastq"],
        ["group2", "duplicate", "short_se", "inputs/two.fastq"],
    ]

    assert [row["group"] for row in _load(_write_sheet(tmp_path, rows=rows))] == [
        "group1",
        "group2",
    ]


def test_normalized_basename_over_220_characters_is_rejected(tmp_path):
    group = "g" * 110
    sample = "s" * 110

    with pytest.raises(
        RunSheetError,
        match="Normalized row basename exceeds 220 characters.*221 characters",
    ):
        _load(
            _write_sheet(
                tmp_path,
                rows=[[group, sample, "short_se", "inputs/reads.fastq"]],
            )
        )


def test_220_character_normalized_basename_is_allowed(tmp_path):
    group = "g" * 110
    sample = "s" * 109

    assert (
        _load(
            _write_sheet(
                tmp_path,
                rows=[[group, sample, "short_se", "inputs/reads.fastq"]],
            )
        )[0]["sample"]
        == sample
    )


def test_one_caller_owned_validator_is_reused_for_every_row(tmp_path):
    validated = []

    def validate_row(row):
        validated.append(row)
        SAMPLES_VALIDATOR.validate(row)

    rows = [
        ["group1", f"sample{index}", "short_se", f"inputs/{index}.fastq"]
        for index in range(1, 5)
    ]

    parsed = load_run_sheet(_write_sheet(tmp_path, rows=rows), validate_row)

    assert len(validated) == len(rows)
    assert [id(row) for row in validated] != [id(row) for row in parsed]
    assert [row["sample"] for row in validated] == [row["sample"] for row in parsed]
