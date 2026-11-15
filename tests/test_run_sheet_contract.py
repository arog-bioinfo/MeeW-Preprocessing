import json
import random
import string
import subprocess
from pathlib import Path

import pytest
import yaml

ROOT = Path(__file__).resolve().parents[1]
SNAKEFILE = ROOT / "workflow" / "Snakefile"


def _write_case(tmp_path, header, rows, *, auto_preprocess_enabled=False):
    workdir = tmp_path / "workdir"
    workdir.mkdir()
    run_sheet = workdir / "run.tsv"
    run_sheet.write_text(
        "\t".join(header) + "\n" + "".join("\t".join(row) + "\n" for row in rows),
        encoding="utf-8",
    )

    for row in rows:
        for column in ("read1", "read2", "contam"):
            if column not in header:
                continue
            value = row[header.index(column)]
            if value:
                path = workdir / value
                path.parent.mkdir(parents=True, exist_ok=True)
                path.touch()

    auto_preprocess = {"enabled": False}
    if auto_preprocess_enabled:
        auto_preprocess = yaml.safe_load(
            (ROOT / "config" / "config.yaml").read_text(encoding="utf-8")
        )["auto_preprocess"]

    config = {
        "run_sheet": str(run_sheet),
        "auto_preprocess": auto_preprocess,
        "sequali": {"extra": ""},
        "fastp": {"adapters": "", "extra": ""},
        "chopper": {"extra": ""},
        "porechop_abi": {
            "extra": "",
            "ab_initio_config": "",
            "custom_adapters": "",
        },
    }
    config_path = workdir / "config.yaml"
    config_path.write_text(yaml.safe_dump(config, sort_keys=False), encoding="utf-8")
    return workdir, config_path


def _snakemake(
    tmp_path,
    header,
    rows,
    *extra,
    snakefile=SNAKEFILE,
    auto_preprocess_enabled=False,
):
    workdir, config_path = _write_case(
        tmp_path,
        header,
        rows,
        auto_preprocess_enabled=auto_preprocess_enabled,
    )
    return subprocess.run(
        [
            "snakemake",
            "--snakefile",
            str(snakefile),
            "--directory",
            str(workdir),
            "--configfile",
            str(config_path),
            "--replace-workflow-config",
            "--nolock",
            "--cores",
            "1",
            "--dry-run",
            "--strict-dag-evaluation",
            "functions",
            *extra,
        ],
        cwd=ROOT,
        check=False,
        capture_output=True,
        text=True,
    )


def _output(result):
    return result.stdout + result.stderr


def _write_module_harness(tmp_path, config_path, *, results=None, logs=None):
    pathvars = []
    if results is not None:
        pathvars.append(f"        results={results!r},")
    if logs is not None:
        pathvars.append(f"        logs={logs!r},")
    pathvars_block = ""
    if pathvars:
        pathvars_block = "    pathvars:\n" + "\n".join(pathvars) + "\n"

    harness = tmp_path / "module.smk"
    harness.write_text(
        f"""import json

configfile: {str(config_path)!r}

module preprocessing:
{pathvars_block}    snakefile:
        {str(SNAKEFILE)!r}
    config:
        config

use rule * from preprocessing as imported_*

print(
    "MODULE_FINAL_TARGETS="
    + json.dumps(preprocessing.FINAL_TARGETS(None))
)

MODULE_LOG_PATHS = []
for row in preprocessing.RUN_SHEET_ROWS:
    row_wildcards = {{"group": row["group"], "sample": row["sample"]}}
    for qc_stage in ["raw", "processed"]:
        MODULE_LOG_PATHS.append(
            str(
                rules.imported_sequali_qc.log[0].apply_wildcards(
                    {{**row_wildcards, "qc_stage": qc_stage}}
                )
            )
        )
    if row["data_type"] == "short_se":
        MODULE_LOG_PATHS.append(
            str(
                rules.imported_process_short_se_fastp.log[0].apply_wildcards(
                    row_wildcards
                )
            )
        )
    elif row["data_type"] == "short_pe":
        MODULE_LOG_PATHS.append(
            str(
                rules.imported_process_short_pe_fastp.log[0].apply_wildcards(
                    row_wildcards
                )
            )
        )
    else:
        MODULE_LOG_PATHS.append(
            str(
                rules.imported_process_long_chopper.log[0].apply_wildcards(
                    row_wildcards
                )
            )
        )
        if row["data_type"] in preprocessing.PORECHOP_ABI_LONG_READ_TYPES:
            MODULE_LOG_PATHS.append(
                str(
                    rules.imported_process_long_porechop_abi.log[0].apply_wildcards(
                        row_wildcards
                    )
                )
            )
    if (
        preprocessing.AUTO_PREPROCESS_ENABLED
        and row["data_type"] in preprocessing.SHORT_READ_TYPES
    ):
        MODULE_LOG_PATHS.append(
            str(
                rules.imported_infer_auto_preprocess.log[0].apply_wildcards(
                    row_wildcards
                )
            )
        )
print("MODULE_LOG_PATHS=" + json.dumps(MODULE_LOG_PATHS))
""",
        encoding="utf-8",
    )
    return harness


def _expected_final_targets(rows, results):
    targets = []
    for group, sample, _, *_ in rows:
        basename = f"{group}_{sample}"
        root = results.format(group=group)
        targets.extend(
            [
                f"{root}/qc/raw/{basename}.sequali.html",
                f"{root}/qc/raw/{basename}.sequali.json",
            ]
        )

    for group, sample, data_type, *_ in rows:
        basename = f"{group}_{sample}"
        root = results.format(group=group)
        if data_type == "short_pe":
            targets.extend(
                [
                    f"{root}/reads/{basename}_R1.fastq",
                    f"{root}/reads/{basename}_R2.fastq",
                ]
            )
        else:
            targets.append(f"{root}/reads/{basename}.fastq")

    for group, sample, _, *_ in rows:
        basename = f"{group}_{sample}"
        root = results.format(group=group)
        targets.extend(
            [
                f"{root}/qc/processed/{basename}.sequali.html",
                f"{root}/qc/processed/{basename}.sequali.json",
            ]
        )
    return targets


def _expected_logs(rows, logs, *, auto_preprocess_enabled):
    expected = []
    for group, sample, data_type, *_ in rows:
        basename = f"{group}_{sample}"
        root = logs.format(group=group)
        expected.extend(
            [
                f"{root}/sequali/raw/{basename}.sequali.log",
                f"{root}/sequali/processed/{basename}.sequali.log",
            ]
        )
        if data_type in {"short_se", "short_pe"}:
            expected.append(f"{root}/fastp/{basename}.fastp.log")
            if auto_preprocess_enabled:
                expected.append(f"{root}/auto_preprocess/{basename}.auto_params.log")
        else:
            expected.append(f"{root}/chopper/{basename}.chopper.log")
            if not data_type.startswith("pacbio-"):
                expected.append(f"{root}/porechop_abi/{basename}.porechop_abi.log")
    return expected


@pytest.mark.parametrize(
    (
        "results_override",
        "logs_override",
        "expected_results",
        "expected_logs_root",
        "auto_preprocess_enabled",
    ),
    [
        (None, None, "results/{group}", "results/{group}/logs", False),
        (
            "results/{group,[A-Za-z0-9][A-Za-z0-9._-]*}/preprocessing",
            None,
            "results/{group}/preprocessing",
            "results/{group}/preprocessing/logs",
            True,
        ),
        (
            "relocated-results/{group}/preprocessing",
            "independent-logs/{group}/preprocessing",
            "relocated-results/{group}/preprocessing",
            "independent-logs/{group}/preprocessing",
            True,
        ),
    ],
    ids=["standalone-defaults", "results-with-nested-logs", "independent-logs"],
)
def test_rule_derived_final_targets_and_log_pathvar_overrides(
    tmp_path,
    results_override,
    logs_override,
    expected_results,
    expected_logs_root,
    auto_preprocess_enabled,
):
    header = ["group", "sample", "data_type", "read1", "read2"]
    rows = [
        ["alpha", "se", "short_se", "inputs/se.fastq", ""],
        ["alpha", "pe", "short_pe", "inputs/pe_R1.fastq", "inputs/pe_R2.fastq"],
        ["alpha", "nano", "nano-raw", "inputs/nano.fastq", ""],
        ["beta", "pacbio", "pacbio-hifi", "inputs/pacbio.fastq", ""],
    ]
    workdir, config_path = _write_case(
        tmp_path,
        header,
        rows,
        auto_preprocess_enabled=auto_preprocess_enabled,
    )
    harness = _write_module_harness(
        tmp_path,
        config_path,
        results=results_override,
        logs=logs_override,
    )

    result = subprocess.run(
        [
            "snakemake",
            "--snakefile",
            str(harness),
            "--directory",
            str(workdir),
            "--nolock",
            "--cores",
            "1",
            "--dry-run",
            "--printshellcmds",
            "imported_all",
        ],
        cwd=ROOT,
        check=False,
        capture_output=True,
        text=True,
    )

    assert result.returncode == 0, _output(result)
    rendered = _output(result)
    reported = [
        json.loads(line.removeprefix("MODULE_FINAL_TARGETS="))
        for line in result.stdout.splitlines()
        if line.startswith("MODULE_FINAL_TARGETS=")
    ]
    assert reported == [_expected_final_targets(rows, expected_results)]
    reported_logs = [
        json.loads(line.removeprefix("MODULE_LOG_PATHS="))
        for line in result.stdout.splitlines()
        if line.startswith("MODULE_LOG_PATHS=")
    ]
    expected_logs = _expected_logs(
        rows,
        expected_logs_root,
        auto_preprocess_enabled=auto_preprocess_enabled,
    )
    assert reported_logs == [expected_logs]

    final_targets = reported[0]
    if results_override is None:
        standalone = subprocess.run(
            [
                "snakemake",
                "--snakefile",
                str(SNAKEFILE),
                "--directory",
                str(workdir),
                "--configfile",
                str(config_path),
                "--replace-workflow-config",
                "--nolock",
                "--cores",
                "1",
                "--dry-run",
                "--strict-dag-evaluation",
                "functions",
            ],
            cwd=ROOT,
            check=False,
            capture_output=True,
            text=True,
        )
        assert standalone.returncode == 0, _output(standalone)
        standalone_lines = standalone.stdout.splitlines()
        all_rule_index = standalone_lines.index("rule all:")
        standalone_targets = standalone_lines[all_rule_index + 1].strip()
        assert standalone_targets.startswith("input: ")
        assert standalone_targets.removeprefix("input: ").split(", ") == final_targets

    assert not any(".fastp." in path for path in final_targets)
    assert not any("auto_params" in path for path in final_targets)
    assert not any(".log" in path for path in final_targets)
    assert not any(".porechop.fastq" in path for path in final_targets)
    if logs_override is not None:
        assert not any(
            path.startswith(f"{expected_results.format(group='alpha')}/logs/")
            for path in reported_logs[0]
        )

    observed_paths = final_targets + reported_logs[0]
    for group, sample, data_type, *_ in rows:
        basename = f"{group}_{sample}"
        old_root = expected_results.format(group=group)
        assert f"{old_root}/qc/raw/{basename}.sequali.log" not in observed_paths
        assert f"{old_root}/qc/processed/{basename}.sequali.log" not in observed_paths
        assert f"{old_root}/reads/{basename}.auto_params.log" not in observed_paths
        if data_type in {"short_se", "short_pe"}:
            assert f"{old_root}/qc/processed/{basename}.fastp.log" not in observed_paths
            if auto_preprocess_enabled:
                assert f"{old_root}/reads/{basename}.auto_params.json" in rendered
        else:
            assert (
                f"{old_root}/qc/processed/{basename}.chopper.log" not in observed_paths
            )
            assert (
                f"{old_root}/qc/processed/{basename}.porechop_abi.log"
                not in observed_paths
            )


def test_mixed_groups_repeated_modes_and_generated_names_route_every_row(tmp_path):
    header = ["group", "sample", "data_type", "read1", "read2", "contam"]
    rows = [
        ["alpha", "", "short_se", "inputs/a1.fastq", "", ""],
        ["alpha", "", "short_se", "inputs/a2.fastq", "", ""],
        [
            "alpha",
            "short_se_1",
            "pacbio-raw",
            "inputs/a-long.fastq",
            "",
            "inputs/contam.fa",
        ],
        [
            "alpha",
            "",
            "short_pe",
            "inputs/a_R1.fastq",
            "inputs/a_R2.fastq",
            "",
        ],
        ["beta", "", "short_se", "inputs/b.fastq", "", ""],
        ["beta", "", "pacbio-raw", "inputs/b-long.fastq", "", ""],
    ]

    result = _snakemake(tmp_path, header, rows, "--printshellcmds")

    assert result.returncode == 0, _output(result)
    rendered = _output(result)
    expected_reads = {
        "results/alpha/reads/alpha_short_se_2.fastq",
        "results/alpha/reads/alpha_short_se_3.fastq",
        "results/alpha/reads/alpha_short_se_1.fastq",
        "results/alpha/reads/alpha_short_pe_1_R1.fastq",
        "results/alpha/reads/alpha_short_pe_1_R2.fastq",
        "results/beta/reads/beta_short_se_1.fastq",
        "results/beta/reads/beta_pacbio-raw_1.fastq",
    }
    for path in expected_reads:
        assert path in rendered
    expected_basenames = {
        ("alpha", "alpha_short_se_2"),
        ("alpha", "alpha_short_se_3"),
        ("alpha", "alpha_short_se_1"),
        ("alpha", "alpha_short_pe_1"),
        ("beta", "beta_short_se_1"),
        ("beta", "beta_pacbio-raw_1"),
    }
    for group, basename in expected_basenames:
        for stage in ("raw", "processed"):
            for extension in ("html", "json"):
                assert (
                    f"results/{group}/qc/{stage}/{basename}.sequali.{extension}"
                    in rendered
                )
    assert "process_short_se_fastp        3" in rendered
    assert "process_short_pe_fastp        1" in rendered
    assert "process_long_chopper          2" in rendered
    assert "AmbiguousRuleException" not in rendered
    assert "InputFunctionException" not in rendered


def test_same_sample_text_in_different_groups_routes_by_explicit_pair(tmp_path):
    result = _snakemake(
        tmp_path,
        ["group", "sample", "data_type", "read1"],
        [
            ["short_group", "reads", "short_se", "inputs/short.fastq"],
            ["long_group", "reads", "pacbio-hifi", "inputs/long.fastq"],
        ],
    )

    assert result.returncode == 0, _output(result)
    rendered = _output(result)
    assert "results/short_group/reads/short_group_reads.fastq" in rendered
    assert "results/long_group/reads/long_group_reads.fastq" in rendered
    assert "AmbiguousRuleException" not in rendered
    assert "InputFunctionException" not in rendered


def test_overlapping_group_sample_underscore_splits_are_unambiguous(tmp_path):
    result = _snakemake(
        tmp_path,
        ["group", "sample", "data_type", "read1"],
        [
            ["a", "b_c", "short_se", "inputs/first.fastq"],
            ["a_b", "c", "short_se", "inputs/second.fastq"],
        ],
    )

    assert result.returncode == 0, _output(result)
    rendered = _output(result)
    assert "results/a/reads/a_b_c.fastq" in rendered
    assert "results/a_b/reads/a_b_c.fastq" in rendered
    assert "AmbiguousRuleException" not in rendered
    assert "InputFunctionException" not in rendered


def test_220_character_normalized_basename_is_accepted_by_dag(tmp_path):
    generator = random.Random(1)
    alphabet = string.ascii_letters + string.digits
    group = "g" + "".join(generator.choices(alphabet, k=109))
    sample = "s" + "".join(generator.choices(alphabet, k=108))
    result = _snakemake(
        tmp_path,
        ["group", "sample", "data_type", "read1"],
        [[group, sample, "short_se", "inputs/reads.fastq"]],
    )

    assert result.returncode == 0, _output(result)
    assert f"results/{group}/reads/{group}_{sample}.fastq" in _output(result)


def test_parser_error_is_wrapped_for_workflow_users(tmp_path):
    result = _snakemake(
        tmp_path,
        ["group", "sample", "data_type", "read1"],
        [["group1", "sample1", "short_se", "inputs/reads.fastq", "extra"]],
    )

    assert result.returncode != 0
    assert "Run sheet row 2 has more fields than the header" in _output(result)


def test_schema_error_preserves_snakemake_validation_diagnostic(tmp_path):
    result = _snakemake(
        tmp_path,
        ["group", "sample", "data_type", "read1"],
        [["", "sample1", "short_se", "inputs/reads.fastq"]],
    )

    assert result.returncode != 0
    rendered = _output(result)
    assert "Error validating config file" in rendered
    assert "group" in rendered
    assert "non-empty" in rendered or "too short" in rendered


def test_module_skip_validation_controls_samples_schema_validation(tmp_path):
    header = ["group", "sample", "data_type", "read1", "module_owned_extra"]
    rows = [["group1", "sample1", "short_se", "inputs/reads.fastq", "ignored"]]
    workdir, config_path = _write_case(tmp_path, header, rows)
    harness = tmp_path / "module.smk"
    harness.write_text(
        f"""configfile: {str(config_path)!r}

module preprocessing:
    snakefile:
        {str(SNAKEFILE)!r}
    config:
        config
    skip_validation:
        True

use rule * from preprocessing exclude all as preprocessing_*
use rule all from preprocessing as preprocessing_all
""",
        encoding="utf-8",
    )

    result = subprocess.run(
        [
            "snakemake",
            "--snakefile",
            str(harness),
            "--directory",
            str(workdir),
            "--nolock",
            "--cores",
            "1",
            "--dry-run",
            "preprocessing_all",
        ],
        cwd=ROOT,
        check=False,
        capture_output=True,
        text=True,
    )

    assert result.returncode == 0, _output(result)
    assert "results/group1/reads/group1_sample1.fastq" in _output(result)
