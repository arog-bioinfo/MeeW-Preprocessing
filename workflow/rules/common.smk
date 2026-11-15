import importlib.util
import json
import re
from pathlib import Path

from jsonschema.exceptions import ValidationError
from snakemake.exceptions import WorkflowError
from snakemake.utils import validate


def _load_dag_helper(module_name, source):
    """Load a module-owned DAG helper from Snakemake's resolved source path."""
    spec = importlib.util.spec_from_file_location(module_name, source)
    if spec is None or spec.loader is None:
        raise ImportError(f"Cannot load Preprocessing helper from {source}")
    helper = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(helper)
    return helper


_run_sheet_helpers = _load_dag_helper(
    "_meew_preprocessing_run_sheet",
    Path(workflow.source_path("../scripts/common/run_sheet.py")),
)
expand_configured_path = _run_sheet_helpers.expand_configured_path

# --- Data-type and run-sheet constants ---------------------------------------

FLYE_LONG_READ_TYPES = {
    "nano-raw",
    "nano-corr",
    "nano-hq",
    "pacbio-raw",
    "pacbio-corr",
    "pacbio-hifi",
}
SHORT_READ_TYPES = {"short_se", "short_pe"}
LONG_PREPROCESS_TYPES = {"long"} | FLYE_LONG_READ_TYPES
PORECHOP_ABI_LONG_READ_TYPES = {"long", "nano-raw", "nano-corr", "nano-hq"}


# --- Raw config validation and path normalization ---------------------------


validate(config, schema="../schemas/config.schema.yaml")

config["run_sheet"] = expand_configured_path(config["run_sheet"])
for section, key in [
    ("porechop_abi", "ab_initio_config"),
    ("porechop_abi", "custom_adapters"),
]:
    if key in config.get(section, {}):
        config[section][key] = expand_configured_path(config[section][key])


# --- Run-sheet parsing and row validation -----------------------------------


if workflow.modifier.skip_validation:
    _SAMPLES_ROW_VALIDATOR = None
else:
    _SAMPLES_ROW_VALIDATOR = _run_sheet_helpers.build_samples_row_validator(
        Path(workflow.source_path("../schemas/samples.schema.yaml"))
    )


def _validate_sample_row(row):
    if _SAMPLES_ROW_VALIDATOR is None:
        return
    try:
        _SAMPLES_ROW_VALIDATOR.validate(row)
    except ValidationError as error:
        raise WorkflowError("Error validating config file.", error)


def _load_run_sheet(run_sheet):
    try:
        return _run_sheet_helpers.load_run_sheet(run_sheet, _validate_sample_row)
    except _run_sheet_helpers.RunSheetError as error:
        raise WorkflowError(str(error)) from error


# --- Normalized run state ----------------------------------------------------


RUN_SHEET_ROWS = _load_run_sheet(config["run_sheet"])
ROW_BY_KEY = {(row["group"], row["sample"]): row for row in RUN_SHEET_ROWS}
ROW_KEYS = tuple(ROW_BY_KEY)
GROUPS = tuple(dict.fromkeys(row["group"] for row in RUN_SHEET_ROWS))
GROUP_PATTERN = "|".join(re.escape(group) for group in GROUPS)
RESULTS_ROOT = f"results/{{group,{GROUP_PATTERN}}}"
LOGS_ROOT = "<results>/logs"
HAS_SHORT_SE = any(row["data_type"] == "short_se" for row in RUN_SHEET_ROWS)
HAS_SHORT_PE = any(row["data_type"] == "short_pe" for row in RUN_SHEET_ROWS)
HAS_LONG = any(row["data_type"] in LONG_PREPROCESS_TYPES for row in RUN_SHEET_ROWS)


def _row_sample_pattern(data_types=None):
    """Constrain a sample wildcard to exact configured (group, sample) pairs."""
    alternatives = [
        f"(?<=/{re.escape(row['group'])}_){re.escape(row['sample'])}"
        for row in RUN_SHEET_ROWS
        if data_types is None or row["data_type"] in data_types
    ]
    return "(?:" + "|".join(alternatives) + ")"


ALL_ROW_SAMPLE_PATTERN = _row_sample_pattern()
SHORT_SE_ROW_SAMPLE_PATTERN = _row_sample_pattern({"short_se"})
SHORT_PE_ROW_SAMPLE_PATTERN = _row_sample_pattern({"short_pe"})
SHORT_ROW_SAMPLE_PATTERN = _row_sample_pattern(SHORT_READ_TYPES)
LONG_ROW_SAMPLE_PATTERN = _row_sample_pattern(LONG_PREPROCESS_TYPES)


# --- Readset and path helpers ------------------------------------------------


def row_record(group, sample):
    try:
        return ROW_BY_KEY[(group, sample)]
    except KeyError as error:
        raise WorkflowError(
            f"Unknown Preprocessing row: group={group!r}, sample={sample!r}"
        ) from error


def row_basename(group, sample):
    row_record(group, sample)
    return f"{group}_{sample}"


def _row_record_for_wildcards(wildcards):
    return row_record(wildcards.group, wildcards.sample)


def _require_row_type(wildcards, allowed_types, route):
    row = _row_record_for_wildcards(wildcards)
    if row["data_type"] not in allowed_types:
        raise WorkflowError(
            f"Preprocessing row ({wildcards.group!r}, {wildcards.sample!r}) has "
            f"data_type={row['data_type']!r}, which is not valid for {route}."
        )
    return row


def processed_read_paths(group, sample, results_root=None):
    data_type = row_record(group, sample)["data_type"]
    basename = row_basename(group, sample)
    root = results_root or "<results>"
    if data_type == "short_se":
        return f"{root}/reads/{basename}.fastq"
    if data_type == "short_pe":
        return [
            f"{root}/reads/{basename}_R1.fastq",
            f"{root}/reads/{basename}_R2.fastq",
        ]
    if data_type in LONG_PREPROCESS_TYPES:
        return f"{root}/reads/{basename}.fastq"
    raise WorkflowError(f"Unsupported data_type for processed reads: {data_type}")


def chopper_inputs(wildcards):
    row = _require_row_type(
        wildcards,
        LONG_PREPROCESS_TYPES,
        "the Chopper route",
    )
    data_type = row["data_type"]
    basename = row_basename(wildcards.group, wildcards.sample)

    if data_type in PORECHOP_ABI_LONG_READ_TYPES:
        fq = f"<results>/reads/{basename}.porechop.fastq"
    elif data_type.startswith("pacbio-"):
        fq = row["read1"]
    else:
        raise WorkflowError(f"Unsupported data_type for Chopper input: {data_type}")

    inputs = {"fq": fq}
    if row.get("contam"):
        inputs["contam"] = row["contam"]
    return inputs


def sequali_qc_reads(wildcards):
    row = _row_record_for_wildcards(wildcards)

    if wildcards.qc_stage == "raw":
        if row["data_type"] == "short_pe":
            return [row["read1"], row["read2"]]
        return row["read1"]

    return processed_read_paths(wildcards.group, wildcards.sample)


# --- Tool configuration ------------------------------------------------------

FASTP_ADAPTERS = config.get("fastp", {}).get("adapters", "")
FASTP_EXTRA = config.get("fastp", {}).get("extra", "")
SEQUALI_EXTRA = config.get("sequali", {}).get("extra", "")
CHOPPER_EXTRA = config.get("chopper", {}).get("extra", "")
PORECHOP_ABI_EXTRA = config.get("porechop_abi", {}).get("extra", "")
PORECHOP_ABI_CONFIG = config.get("porechop_abi", {}).get("ab_initio_config", "")
PORECHOP_ABI_CUSTOM_ADAPTERS = config.get("porechop_abi", {}).get("custom_adapters", "")
AUTO_PREPROCESS_CONFIG = config["auto_preprocess"]
AUTO_PREPROCESS_ENABLED = AUTO_PREPROCESS_CONFIG["enabled"] and any(
    row["data_type"] in SHORT_READ_TYPES for row in RUN_SHEET_ROWS
)


def auto_preprocess_data_type(wildcards):
    return _require_row_type(
        wildcards,
        SHORT_READ_TYPES,
        "automatic short-read inference",
    )["data_type"]


def _short_fastp_inputs(wildcards, expected_type):
    row = _require_row_type(
        wildcards,
        {expected_type},
        f"the {expected_type} fastp route",
    )
    sample = [row["read1"]]
    if row["data_type"] == "short_pe":
        sample.append(row["read2"])

    inputs = {"sample": sample}
    if AUTO_PREPROCESS_ENABLED:
        basename = row_basename(wildcards.group, wildcards.sample)
        inputs["auto_params"] = f"<results>/reads/{basename}.auto_params.json"
    return inputs


def short_se_fastp_inputs(wildcards):
    return _short_fastp_inputs(wildcards, "short_se")


def short_pe_fastp_inputs(wildcards):
    return _short_fastp_inputs(wildcards, "short_pe")


def long_porechop_inputs(wildcards):
    row = _require_row_type(
        wildcards,
        PORECHOP_ABI_LONG_READ_TYPES,
        "the Porechop_ABI route",
    )
    return {
        "fq": row["read1"],
        **({"ab_initio_config": PORECHOP_ABI_CONFIG} if PORECHOP_ABI_CONFIG else {}),
        **(
            {"custom_adapters": PORECHOP_ABI_CUSTOM_ADAPTERS}
            if PORECHOP_ABI_CUSTOM_ADAPTERS
            else {}
        ),
    }


def _render_cli_args(arg_map):
    parts = []
    for key, value in arg_map.items():
        if isinstance(value, bool):
            if value:
                parts.append(str(key))
            continue
        if value is None:
            continue
        parts.extend([str(key), str(value)])
    return " ".join(parts)


def compose_fastp_extra(wildcards, input):
    manual = (FASTP_EXTRA or "").strip()
    if not AUTO_PREPROCESS_ENABLED:
        return manual

    auto_params_path = input.auto_params
    if not auto_params_path:
        raise WorkflowError(
            "Auto preprocess is enabled but no auto-params path was provided for fastp."
        )

    params_path = Path(auto_params_path)
    if not params_path.exists():
        # During DAG construction/dry-run, this file may not yet exist.
        return manual

    with params_path.open(encoding="utf-8") as handle:
        payload = json.load(handle)

    tool_args = payload.get("tool_args", {})
    if not isinstance(tool_args, dict):
        return manual
    inferred_map = tool_args.get("fastp", {})
    if not isinstance(inferred_map, dict):
        return manual

    inferred = _render_cli_args(inferred_map).strip()
    combined = " ".join(part for part in [inferred, manual] if part)
    return combined


# --- Final-target construction ----------------------------------------------


def qc_targets(qc_stage):
    targets = []
    for row in RUN_SHEET_ROWS:
        wildcards = {
            "group": row["group"],
            "sample": row["sample"],
            "qc_stage": qc_stage,
        }
        targets.extend(
            str(output.apply_wildcards(wildcards))
            for output in [
                rules.sequali_qc.output.html,
                rules.sequali_qc.output.json,
            ]
        )
    return targets


def _processed_read_targets():
    targets = []
    for row in RUN_SHEET_ROWS:
        wildcards = {"group": row["group"], "sample": row["sample"]}
        if row["data_type"] == "short_se":
            outputs = [rules.process_short_se_fastp.output.trimmed]
        elif row["data_type"] == "short_pe":
            outputs = rules.process_short_pe_fastp.output.trimmed
        elif row["data_type"] in LONG_PREPROCESS_TYPES:
            outputs = [rules.process_long_chopper.output.out_reads]
        else:
            raise WorkflowError(
                f"Unsupported data_type for final targets: {row['data_type']}"
            )
        targets.extend(str(output.apply_wildcards(wildcards)) for output in outputs)
    return targets


def FINAL_TARGETS(wildcards):
    """Resolve aggregate targets after producer rules and pathvars are available."""
    return qc_targets("raw") + _processed_read_targets() + qc_targets("processed")
