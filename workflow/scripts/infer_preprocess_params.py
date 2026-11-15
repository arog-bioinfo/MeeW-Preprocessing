"""Infer short-read preprocessing parameters from Sequali reports.

The script writes a machine-readable JSON contract and prints a concise summary;
it does not modify reads or config files directly.
"""

import argparse
import json
import math
import sys
import traceback
from contextlib import redirect_stderr, redirect_stdout
from pathlib import Path
from typing import Any

import yaml
from jsonschema import Draft202012Validator
from jsonschema.exceptions import SchemaError

# ---------------------------------------------------------------------------
# Input/config loading
# ---------------------------------------------------------------------------


def load_schema(path: Path) -> dict[str, Any]:
    if not path.exists():
        raise ValueError(f"Schema file does not exist: {path}")
    try:
        with path.open(encoding="utf-8") as handle:
            schema = yaml.safe_load(handle)
    except yaml.YAMLError as exc:
        raise ValueError(f"Schema YAML parse failed for {path}: {exc}") from exc
    if not isinstance(schema, dict):
        raise ValueError(f"Schema root must be a mapping: {path}")
    try:
        Draft202012Validator.check_schema(schema)
    except SchemaError as exc:
        raise ValueError(f"Invalid JSON Schema {path}: {exc.message}") from exc
    return schema


def schema_error_location(error) -> str:
    location = "$"
    for part in error.absolute_path:
        location += f"[{part}]" if isinstance(part, int) else f".{part}"
    return location


def validate_against_schema(
    instance: Any, schema: dict[str, Any], description: str
) -> None:
    errors = sorted(
        Draft202012Validator(schema).iter_errors(instance),
        key=lambda error: tuple(str(part) for part in error.absolute_path),
    )
    if errors:
        error = errors[0]
        raise ValueError(
            f"{description} failed schema validation at "
            f"{schema_error_location(error)}: {error.message}"
        )


def reject_nonfinite_json_number(value: str) -> None:
    raise ValueError(f"JSON contains unsupported non-finite number: {value}")


def load_json(path: Path, schema: dict[str, Any]) -> dict[str, Any]:
    if not path.exists():
        raise ValueError(f"Sequali JSON does not exist: {path}")
    with path.open(encoding="utf-8") as handle:
        data = json.load(handle, parse_constant=reject_nonfinite_json_number)
    validate_against_schema(data, schema, "Sequali JSON")
    return data


def load_config(path: Path, schema: dict[str, Any]) -> dict[str, Any]:
    if not path.exists():
        raise ValueError(f"Config file does not exist: {path}")
    try:
        with path.open(encoding="utf-8") as handle:
            cfg = yaml.safe_load(handle)
    except yaml.YAMLError as exc:
        raise ValueError(f"Config YAML parse failed: {exc}") from exc
    validate_against_schema(cfg, schema, "Workflow config")
    return cfg["auto_preprocess"]


def load_inline_config(
    serialized: str, config_schema: dict[str, Any]
) -> dict[str, Any]:
    try:
        cfg = json.loads(
            serialized,
            parse_constant=reject_nonfinite_json_number,
        )
    except json.JSONDecodeError as exc:
        raise ValueError(
            f"Inline auto-preprocess config JSON parse failed: {exc}"
        ) from exc
    auto_preprocess_schema = config_schema["properties"]["auto_preprocess"]
    validate_against_schema(
        cfg,
        auto_preprocess_schema,
        "Inline auto-preprocess config",
    )
    return cfg


def load_native_config(
    cfg: dict[str, Any], config_schema: dict[str, Any]
) -> dict[str, Any]:
    """Validate a native config mapping with the inline JSON contract."""
    try:
        serialized = json.dumps(cfg)
    except (TypeError, ValueError) as exc:
        raise ValueError(
            f"Native auto-preprocess config JSON serialization failed: {exc}"
        ) from exc
    return load_inline_config(serialized, config_schema)


# ---------------------------------------------------------------------------
# Shared utility functions
# ---------------------------------------------------------------------------


def parse_inclusive_label(label: str) -> tuple[int, int]:
    if "-" in label:
        start, end = (int(value) for value in label.split("-", 1))
    else:
        start = end = int(label)
    if start < 1 or end < start:
        raise ValueError(
            "Invalid Sequali position label; base position 1 is the first base and "
            f"range labels are inclusive: {label}"
        )
    return start, end


def observed_end_mean_quality(
    data: dict[str, Any],
    quality_section_key: str,
    nominal_length: int,
) -> list[float]:
    mean_series = dict(data[quality_section_key]["end_percentiles"])["mean"]
    return mean_series[-nominal_length:]


def infer_headcrop_from_base_imbalance(
    base_content: dict[str, Any], max_headcrop: int, imbalance_threshold: float
) -> int:
    labels = base_content["x_labels"]
    a_vals = base_content["A"]
    c_vals = base_content["C"]
    g_vals = base_content["G"]
    t_vals = base_content["T"]

    n_positions = min(len(labels), len(a_vals), len(c_vals), len(g_vals), len(t_vals))
    start_idx = n_positions // 2

    i = start_idx
    while i > 0 and (
        abs(a_vals[i] - t_vals[i]) < imbalance_threshold
        and abs(g_vals[i] - c_vals[i]) < imbalance_threshold
    ):
        i -= 1

    at_diff = abs(a_vals[i] - t_vals[i])
    gc_diff = abs(g_vals[i] - c_vals[i])

    # Parsing every examined label preserves validation of the scanned bins.
    for j in range(start_idx, i - 1, -1):
        parse_inclusive_label(labels[j])

    no_trigger = (
        i == 0 and at_diff < imbalance_threshold and gc_diff < imbalance_threshold
    )
    if no_trigger:
        return 0

    raw_headcrop = parse_inclusive_label(labels[i])[1]
    max_allowed = parse_inclusive_label(labels[start_idx])[1]
    return min(raw_headcrop, max_allowed, max_headcrop)


def weighted_mott_contribution(
    mean_quality: float,
    quality_cutoff: float,
    positive_weight: float,
    negative_weight: float,
) -> float:
    difference = quality_cutoff - mean_quality
    if difference > 0:
        return difference * positive_weight
    if difference < 0:
        return difference * negative_weight
    return 0.0


# ---------------------------------------------------------------------------
# Short-read inference
# ---------------------------------------------------------------------------


def infer_short(
    data: dict[str, Any], cfg: dict[str, Any], data_type: str
) -> dict[str, Any]:
    tailcrop_mott_quality_cutoff = cfg["tailcrop_mott_quality_cutoff"]
    mott_positive_weight = cfg.get("mott_positive_weight", 1.0)
    mott_negative_weight = cfg.get("mott_negative_weight", 1.0)
    composition_headcrop_enabled = cfg["composition_headcrop_enabled"]
    tailcrop_min_retained_fraction = cfg["tailcrop_min_retained_fraction"]
    average_quality = cfg["average_quality"]
    minimum_length = cfg["minimum_length"]

    def infer_short_mate(
        quality_section_key: str,
        base_section_key: str,
        nominal_length: int,
    ) -> tuple[dict[str, int | None], bool, bool]:
        headcrop = (
            infer_headcrop_from_base_imbalance(
                data[base_section_key],
                cfg["max_headcrop"],
                cfg["headcrop_base_imbalance_threshold"],
            )
            if composition_headcrop_enabled
            else 0
        )
        tailcrop_retained_length_floor = max(
            math.ceil(nominal_length * tailcrop_min_retained_fraction), 1
        )

        observed_mean_quality = observed_end_mean_quality(
            data,
            quality_section_key,
            nominal_length,
        )

        cumulative_score = 0.0
        best_score = 0.0
        best_tailcrop_span = 0
        # Tail-constrained aggregate modified-Mott: every reversed prefix is a
        # 3' suffix candidate, and cumulative score intentionally never resets.
        for tailcrop_span, mean_quality in enumerate(
            reversed(observed_mean_quality), start=1
        ):
            cumulative_score += weighted_mott_contribution(
                mean_quality,
                tailcrop_mott_quality_cutoff,
                mott_positive_weight,
                mott_negative_weight,
            )
            if cumulative_score > best_score:
                best_score = cumulative_score
                best_tailcrop_span = tailcrop_span

        raw_tailcrop_cutoff = (
            nominal_length - best_tailcrop_span if best_tailcrop_span else None
        )
        safeguarded_tailcrop_cutoff = (
            max(raw_tailcrop_cutoff, tailcrop_retained_length_floor)
            if raw_tailcrop_cutoff is not None
            else None
        )
        safeguard_adjusted = (
            raw_tailcrop_cutoff is not None
            and safeguarded_tailcrop_cutoff != raw_tailcrop_cutoff
        )
        tailcrop_suppressed = (
            safeguarded_tailcrop_cutoff is not None
            and safeguarded_tailcrop_cutoff >= nominal_length
        )
        tailcrop_cutoff = None if tailcrop_suppressed else safeguarded_tailcrop_cutoff
        maximum_length_after_headcrop_bases = (
            max(tailcrop_cutoff - headcrop, 1) if tailcrop_cutoff is not None else None
        )
        recommendation = {
            "headcrop_bases": headcrop,
            "tailcrop_base_position": tailcrop_cutoff,
            "maximum_length_after_headcrop_bases": (
                maximum_length_after_headcrop_bases
            ),
        }
        return recommendation, safeguard_adjusted, tailcrop_suppressed

    read1_nominal_length = data["summary"]["maximum_length"]
    read1_recommendation, r1_safeguard_adjusted, r1_tailcrop_suppressed = (
        infer_short_mate(
            "per_position_mean_quality_and_spread",
            "per_position_base_content",
            read1_nominal_length,
        )
    )
    if data_type == "short_pe":
        read2_nominal_length = data["summary_read2"]["maximum_length"]
        read2_recommendation, r2_safeguard_adjusted, r2_tailcrop_suppressed = (
            infer_short_mate(
                "per_position_mean_quality_and_spread_read2",
                "per_position_base_content_read2",
                read2_nominal_length,
            )
        )
    else:
        read2_recommendation = None
        r2_safeguard_adjusted = False
        r2_tailcrop_suppressed = False

    recommendations = {
        "read1": read1_recommendation,
        "read2": read2_recommendation,
    }

    fastp_args: dict[str, Any] = {}
    maximum_length_args: dict[str, Any] = {}
    for read_number, recommendation in enumerate(
        [read1_recommendation, read2_recommendation], start=1
    ):
        if recommendation is not None and recommendation["headcrop_bases"] > 0:
            fastp_args[f"--trim_front{read_number}"] = recommendation["headcrop_bases"]
        if (
            recommendation is not None
            and recommendation["maximum_length_after_headcrop_bases"] is not None
        ):
            maximum_length_args[f"--max_len{read_number}"] = recommendation[
                "maximum_length_after_headcrop_bases"
            ]

    fastp_args.update(maximum_length_args)
    fastp_args["--average_qual"] = average_quality
    fastp_args["--length_required"] = minimum_length

    tool_args = {"fastp": fastp_args}

    warnings = []
    safeguard_adjusted_mates = [
        mate
        for mate, adjusted in [
            ("read1", r1_safeguard_adjusted),
            ("read2", r2_safeguard_adjusted),
        ]
        if adjusted
    ]
    tailcrop_suppressed_mates = [
        mate
        for mate, suppressed in [
            ("read1", r1_tailcrop_suppressed),
            ("read2", r2_tailcrop_suppressed),
        ]
        if suppressed
    ]

    def format_mates(mates: list[str]) -> str:
        return " and ".join(mates)

    if safeguard_adjusted_mates:
        warnings.append(
            {
                "code": "tailcrop_safeguard_adjusted",
                "message": (
                    "Retained-length safeguard adjusted a raw tailcrop cutoff for "
                    f"{format_mates(safeguard_adjusted_mates)}."
                ),
            }
        )
    if tailcrop_suppressed_mates:
        warnings.append(
            {
                "code": "tailcrop_suppressed_as_no_op",
                "message": (
                    "Retained-length safeguard suppressed the tailcrop as a no-op "
                    f"for {format_mates(tailcrop_suppressed_mates)}."
                ),
            }
        )
    if composition_headcrop_enabled:
        composition_headcrop_mates = ["read1"]
        if data_type == "short_pe":
            composition_headcrop_mates.append("read2")
        warnings.append(
            {
                "code": "composition_headcrop_enabled",
                "message": (
                    "Composition-derived 5' headcrop is explicitly enabled for "
                    f"{format_mates(composition_headcrop_mates)}."
                ),
            }
        )

    return {
        "recommendations": recommendations,
        "tool_args": tool_args,
        "warnings": warnings,
    }


# ---------------------------------------------------------------------------
# Output rendering and entrypoint/error handling
# ---------------------------------------------------------------------------


def render_tool_arguments(tool_args: dict[str, dict[str, Any]]) -> str:
    rendered_tools = []
    for tool, arg_map in tool_args.items():
        parts = [tool]
        for key, value in arg_map.items():
            if isinstance(value, bool):
                if value:
                    parts.append(str(key))
                continue
            if value is None:
                continue
            parts.extend([str(key), str(value)])
        rendered_tools.append(" ".join(parts))
    return "; ".join(rendered_tools)


def render_recommendations(
    recommendations: dict[str, dict[str, int | None] | None],
) -> str:
    rendered_mates = []
    for mate in ["read1", "read2"]:
        recommendation = recommendations[mate]
        if recommendation is None:
            rendered_mates.append(f"{mate}=null")
            continue
        rendered_values = ", ".join(
            f"{key}={value}" for key, value in recommendation.items()
        )
        rendered_mates.append(f"{mate}({rendered_values})")
    return "; ".join(rendered_mates)


def print_summary(payload: dict[str, Any]) -> None:
    warnings = payload["warnings"]

    print("Auto-preprocess inference complete")
    print(f"Recommendations: {render_recommendations(payload['recommendations'])}")
    print(f"Tool arguments: {render_tool_arguments(payload['tool_args'])}")
    if warnings:
        print(
            "Warnings: "
            + " | ".join(
                f"{warning['code']}: {warning['message']}" for warning in warnings
            )
        )
    else:
        print("Warnings: none")


def run_inference(
    *,
    sequali_json: Path,
    data_type: str,
    config_schema_path: Path,
    sequali_schema_path: Path,
    output_path: Path,
    config_path: Path | None = None,
    config_json: str | None = None,
    auto_preprocess_config: dict[str, Any] | None = None,
) -> dict[str, Any]:
    print("Auto-preprocess inference started", flush=True)

    config_schema = load_schema(config_schema_path)
    sequali_schema = load_schema(sequali_schema_path)
    data = load_json(sequali_json, sequali_schema)
    if data_type == "short_pe" and "summary_read2" not in data:
        raise ValueError(
            "Sequali JSON does not contain the complete read2 section set required "
            "for data_type=short_pe"
        )

    config_sources = sum(
        source is not None
        for source in (config_path, config_json, auto_preprocess_config)
    )
    if config_sources != 1:
        raise ValueError("Exactly one auto-preprocess config source is required")
    if config_json is not None:
        cfg = load_inline_config(config_json, config_schema)
    elif config_path is not None:
        cfg = load_config(config_path, config_schema)
    else:
        cfg = load_native_config(auto_preprocess_config, config_schema)

    if cfg.get("enabled", True) is False:
        raise ValueError(
            "Auto-preprocess inference cannot run when auto_preprocess.enabled is false"
        )

    payload = infer_short(data, cfg, data_type)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w") as handle:
        json.dump(payload, handle, indent=2)
        handle.write("\n")

    print_summary(payload)
    return payload


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Infer short-read preprocessing parameters from Sequali JSON"
    )
    parser.add_argument("--sequali-json", required=True)
    parser.add_argument("--data-type", required=True, choices=["short_se", "short_pe"])
    config_group = parser.add_mutually_exclusive_group(required=True)
    config_group.add_argument("--config")
    config_group.add_argument("--config-json")
    parser.add_argument("--config-schema", required=True)
    parser.add_argument("--sequali-schema", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    run_inference(
        sequali_json=Path(args.sequali_json),
        data_type=args.data_type,
        config_path=Path(args.config) if args.config is not None else None,
        config_json=args.config_json,
        config_schema_path=Path(args.config_schema),
        sequali_schema_path=Path(args.sequali_schema),
        output_path=Path(args.output),
    )


def snakemake_main(snakemake_object: Any) -> None:
    log_path = Path(str(snakemake_object.log[0]))
    log_path.parent.mkdir(parents=True, exist_ok=True)

    with log_path.open("w", encoding="utf-8") as log_handle:
        with redirect_stdout(log_handle), redirect_stderr(log_handle):
            try:
                run_inference(
                    sequali_json=Path(str(snakemake_object.input.sequali_json)),
                    data_type=str(snakemake_object.params.data_type),
                    auto_preprocess_config=dict(
                        snakemake_object.params.auto_preprocess_config
                    ),
                    config_schema_path=Path(str(snakemake_object.input.config_schema)),
                    sequali_schema_path=Path(
                        str(snakemake_object.input.sequali_schema)
                    ),
                    output_path=Path(str(snakemake_object.output.auto_params)),
                )
            except ValueError as exc:
                print(f"ERROR: {exc}", file=sys.stderr, flush=True)
                raise
            except Exception:
                traceback.print_exc()
                raise


if __name__ == "__main__":
    if "snakemake" in globals():
        snakemake_main(snakemake)
    else:
        try:
            main()
        except ValueError as exc:
            raise SystemExit(f"ERROR: {exc}")
