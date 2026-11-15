import json
from pathlib import Path

import pytest


def _populate_inventory(result_root, inventory):
    for relative in sorted(
        inventory.directories, key=lambda value: (value.count("/"), value)
    ):
        (result_root / relative).mkdir(parents=True, exist_ok=True)
    for relative in inventory.files:
        (result_root / relative).write_text("fixture\n", encoding="utf-8")


def _case_tree(tmp_path, assertions, case_id):
    case = assertions.CASES[case_id]
    inventory = assertions.CASE_INVENTORIES[case_id]
    workdir = tmp_path / "workdir"
    result_root = workdir / "results" / case.sample
    result_root.mkdir(parents=True)
    _populate_inventory(result_root, inventory)
    return case, inventory, workdir, result_root


@pytest.mark.parametrize(
    "case_id",
    [
        "short-se",
        "short-pe-disabled",
        "short-pe-infer",
        "long",
        "nano-raw",
        "pacbio-raw",
        "hybrid",
    ],
)
def test_case_inventory_routes_logs_to_exact_tool_directories(
    workflow_assertions_module, case_id
):
    assertions = workflow_assertions_module
    case = assertions.CASES[case_id]
    inventory = assertions.CASE_INVENTORIES[case_id]
    expected_log_directories = {
        "logs",
        "logs/sequali",
        "logs/sequali/raw",
        "logs/sequali/processed",
    }
    if case.short_readsets:
        expected_log_directories.add("logs/fastp")
    if case.long_readsets:
        expected_log_directories.add("logs/chopper")
    if case.porechop_routed_readsets:
        expected_log_directories.add("logs/porechop_abi")
    if case.inference is not None:
        expected_log_directories.add("logs/auto_preprocess")

    assert {path for path in inventory.directories if path.startswith("logs")} == (
        expected_log_directories
    )
    for readset in case.readsets:
        for stage in ("raw", "processed"):
            assert f"logs/sequali/{stage}/{readset}.sequali.log" in inventory.files
            assert f"qc/{stage}/{readset}.sequali.log" in inventory.forbidden
    for readset in case.short_readsets:
        assert f"logs/fastp/{readset}.fastp.log" in inventory.files
        assert f"qc/processed/{readset}.fastp.log" in inventory.forbidden
    for readset in case.long_readsets:
        assert f"logs/chopper/{readset}.chopper.log" in inventory.files
        assert f"qc/processed/{readset}.chopper.log" in inventory.forbidden
    for readset in case.porechop_routed_readsets:
        assert f"logs/porechop_abi/{readset}.porechop_abi.log" in inventory.files
        assert f"qc/processed/{readset}.porechop_abi.log" in inventory.forbidden
    for readset in case.porechop_bypassed_readsets:
        assert f"logs/porechop_abi/{readset}.porechop_abi.log" in inventory.forbidden
        assert f"qc/processed/{readset}.porechop_abi.log" in inventory.forbidden

    inferred_readset = case.inference.readset if case.inference is not None else None
    for readset in case.readsets:
        new_log = f"logs/auto_preprocess/{readset}.auto_params.log"
        old_log = f"reads/{readset}.auto_params.log"
        assert old_log in inventory.forbidden
        if readset == inferred_readset:
            assert new_log in inventory.files
        else:
            assert new_log in inventory.forbidden


@pytest.mark.parametrize(
    "mutation",
    ["missing", "empty", "extra", "wrong-type", "nested", "wrongly-routed"],
)
def test_artifact_inventory_rejects_tree_regressions(
    tmp_path, workflow_assertions_module, mutation
):
    assertions = workflow_assertions_module
    case, inventory, _, result_root = _case_tree(tmp_path, assertions, "short-se")
    processed_log = result_root / "logs/fastp/test_short_se_short.fastp.log"

    if mutation == "missing":
        processed_log.unlink()
    elif mutation == "empty":
        processed_log.write_bytes(b"")
    elif mutation == "extra":
        (result_root / "unexpected.txt").write_text("unexpected\n", encoding="utf-8")
    elif mutation == "wrong-type":
        processed_log.unlink()
        processed_log.mkdir()
    elif mutation == "nested":
        nested = result_root / "reads/nested"
        nested.mkdir()
        (nested / "artifact.fastq").write_text("nested\n", encoding="utf-8")
    elif mutation == "wrongly-routed":
        processed_log.replace(
            result_root / "qc/processed/test_short_se_short.fastp.log"
        )

    with pytest.raises(AssertionError):
        assertions.assert_artifact_inventory(case, inventory, result_root)


def test_artifact_inventory_rejects_forbidden_porechop_temporary(
    tmp_path, workflow_assertions_module
):
    assertions = workflow_assertions_module
    case, inventory, _, result_root = _case_tree(tmp_path, assertions, "long")
    (result_root / "reads/test_long_long.porechop.fastq").write_text(
        "temporary\n", encoding="utf-8"
    )

    with pytest.raises(AssertionError, match="Forbidden workflow artifact"):
        assertions.assert_artifact_inventory(case, inventory, result_root)


def test_golden_fastq_comparison_detects_altered_output(
    tmp_path, workflow_assertions_module
):
    assertions = workflow_assertions_module
    case, inventory, workdir, result_root = _case_tree(tmp_path, assertions, "short-se")
    expected_root = workdir / "expected"
    expected_root.mkdir()
    fastq = b"@read\nAC\n+\n!!\n"
    expected = expected_root / case.read_comparisons[0][1]
    output = result_root / case.read_comparisons[0][0]
    expected.write_bytes(fastq)
    output.write_bytes(fastq)
    assertions.assert_artifact_inventory(case, inventory, result_root)
    assertions.assert_reads(case, workdir, result_root)

    output.write_bytes(b"@read\nAG\n+\n!!\n")

    with pytest.raises(AssertionError, match="Processed reads differ"):
        assertions.assert_reads(case, workdir, result_root)


@pytest.mark.parametrize("mutation", ["payload", "log"])
def test_inference_assertions_detect_exact_content_regressions(
    tmp_path, workflow_assertions_module, mutation
):
    assertions = workflow_assertions_module
    case, inventory, _, result_root = _case_tree(tmp_path, assertions, "short-pe-infer")
    expectation = case.inference
    assert expectation is not None
    json_path = result_root / "reads" / f"{expectation.readset}.auto_params.json"
    log_path = (
        result_root / "logs/auto_preprocess" / f"{expectation.readset}.auto_params.log"
    )
    json_path.write_text(json.dumps(expectation.payload), encoding="utf-8")
    log_path.write_text(expectation.log, encoding="utf-8")
    assertions.assert_artifact_inventory(case, inventory, result_root)
    assertions.assert_inference(case, result_root)

    if mutation == "payload":
        altered = json.loads(json.dumps(expectation.payload))
        altered["tool_args"]["fastp"]["--average_qual"] = 31
        json_path.write_text(json.dumps(altered), encoding="utf-8")
        message = "Unexpected inference payload"
    else:
        log_path.write_text(expectation.log + "altered\n", encoding="utf-8")
        message = "Unexpected inference log content"

    with pytest.raises(AssertionError, match=message):
        assertions.assert_inference(case, result_root)
