from __future__ import annotations

import builtins
import json
import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace

import pytest
import yaml

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "workflow" / "scripts" / "infer_preprocess_params.py"
CONFIG_SCHEMA = ROOT / "workflow" / "schemas" / "config.schema.yaml"
SEQUALI_SCHEMA = ROOT / "workflow" / "schemas" / "sequali_short_read.schema.yaml"


def run_native_script(tmp_path, data_type, sequali, auto_preprocess_config):
    workdir = tmp_path / "native"
    workdir.mkdir()
    sequali_path = workdir / "sequali.json"
    output_path = workdir / "output.json"
    log_path = workdir / "inference.log"
    sequali_path.write_text(json.dumps(sequali), encoding="utf-8")

    serialized_config = json.dumps(auto_preprocess_config)
    snakefile = tmp_path / "inference.smk"
    snakefile.write_text(
        f"""import json

AUTO_PREPROCESS_CONFIG = json.loads({serialized_config!r})

rule native_inference:
    input:
        sequali_json="sequali.json",
        config_schema={str(CONFIG_SCHEMA)!r},
        sequali_schema={str(SEQUALI_SCHEMA)!r},
    output:
        auto_params="output.json",
    log:
        "inference.log",
    params:
        data_type={data_type!r},
        auto_preprocess_config=AUTO_PREPROCESS_CONFIG,
    script:
        {str(SCRIPT)!r}
""",
        encoding="utf-8",
    )
    result = subprocess.run(
        [
            "snakemake",
            "--snakefile",
            str(snakefile),
            "--directory",
            str(workdir),
            "--cores",
            "1",
            "--nolock",
            "output.json",
        ],
        cwd=ROOT,
        check=False,
        capture_output=True,
        text=True,
    )
    return result, output_path, log_path


def run_cli(
    tmp_path,
    data_type,
    sequali,
    auto_preprocess_config,
    *,
    config_mode="config-json",
):
    workdir = tmp_path / "cli"
    workdir.mkdir()
    sequali_path = workdir / "sequali.json"
    output_path = workdir / "output.json"
    sequali_path.write_text(json.dumps(sequali), encoding="utf-8")

    command = [
        sys.executable,
        str(SCRIPT),
        "--sequali-json",
        str(sequali_path),
        "--data-type",
        data_type,
    ]
    if config_mode == "config-json":
        command.extend(["--config-json", json.dumps(auto_preprocess_config)])
    else:
        config_path = workdir / "config.yaml"
        config_path.write_text(
            yaml.safe_dump(
                {
                    "run_sheet": "config/run.tsv",
                    "auto_preprocess": auto_preprocess_config,
                }
            ),
            encoding="utf-8",
        )
        command.extend(["--config", str(config_path)])
    command.extend(
        [
            "--config-schema",
            str(CONFIG_SCHEMA),
            "--sequali-schema",
            str(SEQUALI_SCHEMA),
            "--output",
            str(output_path),
        ]
    )
    result = subprocess.run(
        command,
        cwd=ROOT,
        check=False,
        capture_output=True,
        text=True,
    )
    return result, output_path, command


def end_mean_values(data, section="per_position_mean_quality_and_spread"):
    return next(
        values for name, values in data[section]["end_percentiles"] if name == "mean"
    )


def recommendation(tailcrop, headcrop=0):
    return {
        "headcrop_bases": headcrop,
        "tailcrop_base_position": tailcrop,
        "maximum_length_after_headcrop_bases": (
            max(tailcrop - headcrop, 1) if tailcrop is not None else None
        ),
    }


def test_standard_100_slot_anchored_graph_uses_exact_l_minus_n_cutoff(
    inference_module, make_config, make_sequali
):
    data = make_sequali(read1_length=150, read1_bad=(-4, -3, -2, -1))

    assert len(end_mean_values(data)) == 100
    observed = inference_module.observed_end_mean_quality(
        data,
        "per_position_mean_quality_and_spread",
        data["summary"]["maximum_length"],
    )
    assert observed == end_mean_values(data)

    payload = inference_module.infer_short(
        data,
        make_config(tailcrop_min_retained_fraction=0.0),
        "short_se",
    )

    assert payload["recommendations"] == {
        "read1": recommendation(146),
        "read2": None,
    }
    assert payload["tool_args"]["fastp"]["--max_len1"] == 146


def test_short_reads_observe_only_final_maximum_length_mean_anchors(
    inference_module, make_config, make_sequali
):
    data = make_sequali(read1_length=10)

    observed = inference_module.observed_end_mean_quality(
        data,
        "per_position_mean_quality_and_spread",
        data["summary"]["maximum_length"],
    )
    assert observed == [30.0] * 10
    assert data["per_position_mean_quality_and_spread"]["percentiles"] == [
        ["mean", [1.0]]
    ]
    assert "per_position_quality_distribution" not in data

    payload = inference_module.infer_short(data, make_config(), "short_se")

    assert payload == {
        "recommendations": {
            "read1": recommendation(None),
            "read2": None,
        },
        "tool_args": {"fastp": {"--average_qual": 30, "--length_required": 2}},
        "warnings": [],
    }


@pytest.mark.parametrize(
    ("mean_quality", "expected_score_sign", "expected_cutoff"),
    [
        (19.0, "positive", 9),
        (20.0, "zero", None),
        (21.0, "negative", None),
    ],
)
def test_mott_cutoff_score_sign_controls_terminal_candidate(
    inference_module,
    make_config,
    make_sequali,
    mean_quality,
    expected_score_sign,
    expected_cutoff,
):
    payload = inference_module.infer_short(
        make_sequali(read1_mean_overrides={-1: mean_quality}),
        make_config(tailcrop_min_retained_fraction=0.0),
        "short_se",
    )

    score = 20.0 - mean_quality
    assert {
        "positive": score > 0,
        "zero": score == 0,
        "negative": score < 0,
    }[expected_score_sign]
    assert payload["recommendations"]["read1"]["tailcrop_base_position"] == (
        expected_cutoff
    )


@pytest.mark.parametrize(
    (
        "mean_quality",
        "positive_weight",
        "negative_weight",
        "expected_contribution",
    ),
    [
        (18.0, 3.0, 7.0, 6.0),
        (22.0, 3.0, 7.0, -14.0),
        (20.0, 3.0, 7.0, 0.0),
    ],
)
def test_weighted_mott_contribution_uses_unweighted_difference_sign(
    inference_module,
    mean_quality,
    positive_weight,
    negative_weight,
    expected_contribution,
):
    assert (
        inference_module.weighted_mott_contribution(
            mean_quality,
            20.0,
            positive_weight,
            negative_weight,
        )
        == expected_contribution
    )


@pytest.mark.parametrize("data_type", ["short_se", "short_pe"])
def test_explicit_unit_mott_weights_match_omitted_fallback(
    inference_module, make_config, make_sequali, data_type
):
    data = make_sequali(
        read1_mean_overrides={-3: 19.0, -2: 21.0, -1: 19.0},
        read2_mean_overrides={-4: 0.0, -3: 30.0, -2: 18.0, -1: 19.0},
        read2_length=12,
        read2_bad=() if data_type == "short_pe" else None,
    )
    explicit_config = make_config(tailcrop_min_retained_fraction=0.0)
    omitted_config = dict(explicit_config)
    del omitted_config["mott_positive_weight"]
    del omitted_config["mott_negative_weight"]

    explicit = inference_module.infer_short(data, explicit_config, data_type)
    omitted = inference_module.infer_short(data, omitted_config, data_type)

    assert explicit == omitted
    assert explicit["recommendations"]["read1"] == recommendation(9)
    if data_type == "short_pe":
        assert explicit["recommendations"]["read2"] == recommendation(8)
    else:
        assert explicit["recommendations"]["read2"] is None


def test_changing_positive_mott_weight_alone_changes_poor_quality_evidence(
    inference_module, make_config, make_sequali
):
    data = make_sequali(
        read1_mean_overrides={-3: 19.0, -2: 21.0, -1: 19.0},
    )
    unit_weight = inference_module.infer_short(
        data,
        make_config(tailcrop_min_retained_fraction=0.0),
        "short_se",
    )
    increased_positive_weight = inference_module.infer_short(
        data,
        make_config(
            mott_positive_weight=2.0,
            tailcrop_min_retained_fraction=0.0,
        ),
        "short_se",
    )

    assert unit_weight["recommendations"]["read1"] == recommendation(9)
    assert increased_positive_weight["recommendations"]["read1"] == recommendation(7)


def test_changing_negative_mott_weight_alone_changes_good_quality_offset(
    inference_module, make_config, make_sequali
):
    data = make_sequali(
        read1_mean_overrides={-3: 19.0, -2: 21.0, -1: 19.0},
    )
    unit_weight = inference_module.infer_short(
        data,
        make_config(tailcrop_min_retained_fraction=0.0),
        "short_se",
    )
    decreased_negative_weight = inference_module.infer_short(
        data,
        make_config(
            mott_negative_weight=0.5,
            tailcrop_min_retained_fraction=0.0,
        ),
        "short_se",
    )

    assert unit_weight["recommendations"]["read1"] == recommendation(9)
    assert decreased_negative_weight["recommendations"]["read1"] == recommendation(7)


def test_zero_mott_weights_remove_only_their_matching_evidence(
    inference_module, make_config, make_sequali
):
    positive_removed = inference_module.infer_short(
        make_sequali(read1_mean_overrides={-1: 19.0}),
        make_config(
            mott_positive_weight=0.0,
            tailcrop_min_retained_fraction=0.0,
        ),
        "short_se",
    )
    negative_removed = inference_module.infer_short(
        make_sequali(
            read1_mean_overrides={-3: 19.0, -2: 30.0, -1: 19.0},
        ),
        make_config(
            mott_negative_weight=0.0,
            tailcrop_min_retained_fraction=0.0,
        ),
        "short_se",
    )
    both_removed = inference_module.infer_short(
        make_sequali(read1_mean_overrides={-1: 19.0}),
        make_config(
            mott_positive_weight=0.0,
            mott_negative_weight=0.0,
            tailcrop_min_retained_fraction=0.0,
        ),
        "short_se",
    )

    assert positive_removed["recommendations"]["read1"] == recommendation(None)
    assert negative_removed["recommendations"]["read1"] == recommendation(7)
    assert both_removed["recommendations"]["read1"] == recommendation(None)


def test_weighted_paired_mates_remain_independent(
    inference_module, make_config, make_sequali
):
    payload = inference_module.infer_short(
        make_sequali(
            read1_mean_overrides={-3: 19.0, -2: 21.0, -1: 19.0},
            read2_bad=(),
            read2_mean_overrides={-1: 19.0},
        ),
        make_config(
            mott_positive_weight=2.0,
            tailcrop_min_retained_fraction=0.0,
        ),
        "short_pe",
    )

    assert payload["recommendations"] == {
        "read1": recommendation(7),
        "read2": recommendation(9),
    }
    assert payload["tool_args"]["fastp"]["--max_len1"] == 7
    assert payload["tool_args"]["fastp"]["--max_len2"] == 9


@pytest.mark.parametrize(
    "mean_quality",
    [20.0, 30.0],
)
def test_all_good_or_cutoff_equal_observed_tails_do_not_crop(
    inference_module,
    make_config,
    make_sequali,
    mean_quality,
):
    payload = inference_module.infer_short(
        make_sequali(
            read1_mean_overrides={offset: mean_quality for offset in range(-10, 0)},
        ),
        make_config(tailcrop_min_retained_fraction=0.0),
        "short_se",
    )

    assert payload["recommendations"]["read1"]["tailcrop_base_position"] is None
    assert payload["warnings"] == []


def test_cumulative_good_positions_offset_poor_evidence(
    inference_module, make_config, make_sequali
):
    payload = inference_module.infer_short(
        make_sequali(
            read1_mean_overrides={-3: 19.0, -2: 30.0, -1: 19.0},
        ),
        make_config(tailcrop_min_retained_fraction=0.0),
        "short_se",
    )

    # Reverse scores are +1, -10, +1; the second poor position cannot recover.
    assert payload["recommendations"]["read1"]["tailcrop_base_position"] == 9


def test_maximum_positive_cumulative_suffix_is_selected(
    inference_module, make_config, make_sequali
):
    payload = inference_module.infer_short(
        make_sequali(
            read1_mean_overrides={-4: 0.0, -3: 30.0, -2: 18.0, -1: 19.0},
        ),
        make_config(tailcrop_min_retained_fraction=0.0),
        "short_se",
    )

    # Reverse cumulative scores are 1, 3, -7, and 13; span four is maximal.
    assert payload["recommendations"]["read1"]["tailcrop_base_position"] == 6


def test_strict_best_score_update_keeps_shortest_equal_maximum(
    inference_module, make_config, make_sequali
):
    payload = inference_module.infer_short(
        make_sequali(
            read1_mean_overrides={-3: 19.0, -2: 21.0, -1: 19.0},
        ),
        make_config(tailcrop_min_retained_fraction=0.0),
        "short_se",
    )

    # Reverse cumulative scores are 1, 0, and 1; the first one-base maximum wins.
    assert payload["recommendations"]["read1"]["tailcrop_base_position"] == 9


def test_all_poor_observed_tail_selects_complete_tail_before_safeguards(
    inference_module, make_config, make_sequali
):
    payload = inference_module.infer_short(
        make_sequali(
            read1_length=150,
            read1_bad=tuple(range(-100, 0)),
        ),
        make_config(tailcrop_min_retained_fraction=0.0),
        "short_se",
    )

    assert payload["recommendations"] == {
        "read1": recommendation(50),
        "read2": None,
    }
    assert payload["warnings"] == []


def test_short_all_poor_read_uses_only_observed_anchors(
    inference_module, make_config, make_sequali
):
    payload = inference_module.infer_short(
        make_sequali(read1_length=3, read1_bad=(-3, -2, -1)),
        make_config(tailcrop_min_retained_fraction=0.0),
        "short_se",
    )

    assert payload["recommendations"]["read1"]["tailcrop_base_position"] == 1
    assert [warning["code"] for warning in payload["warnings"]] == [
        "tailcrop_safeguard_adjusted",
    ]


def test_repeated_q30_q30_q30_q19_pattern_crops_only_terminal_q19(
    inference_module, make_config, make_sequali
):
    payload = inference_module.infer_short(
        make_sequali(
            read1_length=8,
            read1_mean_overrides={-5: 19.0, -1: 19.0},
        ),
        make_config(tailcrop_min_retained_fraction=0.0),
        "short_se",
    )

    # 5' to 3': Q30,Q30,Q30,Q19,Q30,Q30,Q30,Q19. Retain the first seven.
    assert payload["recommendations"] == {
        "read1": recommendation(7),
        "read2": None,
    }


def test_more_frequent_poor_positions_produce_longer_crop(
    inference_module, make_config, make_sequali
):
    one_poor = inference_module.infer_short(
        make_sequali(read1_mean_overrides={-1: 19.0}),
        make_config(tailcrop_min_retained_fraction=0.0),
        "short_se",
    )
    two_poor = inference_module.infer_short(
        make_sequali(read1_mean_overrides={-2: 19.0, -1: 19.0}),
        make_config(tailcrop_min_retained_fraction=0.0),
        "short_se",
    )

    assert one_poor["recommendations"]["read1"]["tailcrop_base_position"] == 9
    assert two_poor["recommendations"]["read1"]["tailcrop_base_position"] == 8


def test_tailcrop_safeguard_adjustment(inference_module, make_config, make_sequali):
    payload = inference_module.infer_short(
        make_sequali(read1_bad=(-4, -3, -2, -1)),
        make_config(tailcrop_min_retained_fraction=0.8),
        "short_se",
    )

    assert payload["recommendations"] == {
        "read1": recommendation(8),
        "read2": None,
    }
    assert payload["tool_args"]["fastp"]["--max_len1"] == 8
    assert payload["warnings"] == [
        {
            "code": "tailcrop_safeguard_adjusted",
            "message": (
                "Retained-length safeguard adjusted a raw tailcrop cutoff for read1."
            ),
        }
    ]


def test_one_base_technical_floor(inference_module, make_config, make_sequali):
    payload = inference_module.infer_short(
        make_sequali(read1_length=2, read1_bad=(-2, -1)),
        make_config(tailcrop_min_retained_fraction=0.0),
        "short_se",
    )

    assert payload["recommendations"] == {
        "read1": recommendation(1),
        "read2": None,
    }
    assert payload["tool_args"]["fastp"]["--max_len1"] == 1
    assert [warning["code"] for warning in payload["warnings"]] == [
        "tailcrop_safeguard_adjusted",
    ]


def test_tailcrop_safeguard_no_op_suppression(
    inference_module, make_config, make_sequali
):
    payload = inference_module.infer_short(
        make_sequali(read1_bad=(-4, -3, -2, -1)),
        make_config(tailcrop_min_retained_fraction=1.0),
        "short_se",
    )

    assert payload["recommendations"]["read1"]["tailcrop_base_position"] is None
    assert (
        payload["recommendations"]["read1"]["maximum_length_after_headcrop_bases"]
        is None
    )
    assert "--max_len1" not in payload["tool_args"]["fastp"]
    assert [warning["code"] for warning in payload["warnings"]] == [
        "tailcrop_safeguard_adjusted",
        "tailcrop_suppressed_as_no_op",
    ]


def test_paired_no_op_is_mate_specific_and_does_not_suppress_other_headcrop(
    inference_module, make_config, make_sequali
):
    payload = inference_module.infer_short(
        make_sequali(
            read1_bad=(-4, -3, -2, -1),
            read2_bad=(),
            read2_imbalance=2,
        ),
        make_config(
            composition_headcrop_enabled=True,
            tailcrop_min_retained_fraction=1.0,
        ),
        "short_pe",
    )

    assert payload["recommendations"] == {
        "read1": recommendation(None),
        "read2": recommendation(None, headcrop=3),
    }
    assert payload["tool_args"]["fastp"] == {
        "--trim_front2": 3,
        "--average_qual": 30,
        "--length_required": 2,
    }
    assert payload["warnings"][1] == {
        "code": "tailcrop_suppressed_as_no_op",
        "message": (
            "Retained-length safeguard suppressed the tailcrop as a no-op for " "read1."
        ),
    }
    assert payload["warnings"][2] == {
        "code": "composition_headcrop_enabled",
        "message": (
            "Composition-derived 5' headcrop is explicitly enabled for read1 and "
            "read2."
        ),
    }


def test_paired_mates_emit_independent_tailcrop_arguments(
    inference_module, make_config, make_sequali
):
    payload = inference_module.infer_short(
        make_sequali(read1_bad=(-3, -2, -1), read2_bad=(-2, -1)),
        make_config(tailcrop_min_retained_fraction=0.0),
        "short_pe",
    )

    assert payload["recommendations"] == {
        "read1": recommendation(7),
        "read2": recommendation(8),
    }
    assert payload["tool_args"]["fastp"] == {
        "--max_len1": 7,
        "--max_len2": 8,
        "--average_qual": 30,
        "--length_required": 2,
    }
    assert payload["warnings"] == []
    assert "--trim_tail1" not in payload["tool_args"]["fastp"]
    assert "--trim_tail2" not in payload["tool_args"]["fastp"]


def test_paired_null_tailcrop_does_not_affect_other_mate(
    inference_module, make_config, make_sequali
):
    payload = inference_module.infer_short(
        make_sequali(
            read1_length=10,
            read1_bad=(-3, -2, -1),
            read2_length=12,
            read2_bad=(),
        ),
        make_config(tailcrop_min_retained_fraction=0.75),
        "short_pe",
    )

    assert payload["recommendations"] == {
        "read1": recommendation(8),
        "read2": recommendation(None),
    }
    assert payload["tool_args"]["fastp"]["--max_len1"] == 8
    assert "--max_len2" not in payload["tool_args"]["fastp"]
    assert payload["warnings"][0]["code"] == "tailcrop_safeguard_adjusted"
    assert "read1" in payload["warnings"][0]["message"]
    assert "read2" not in payload["warnings"][0]["message"]


def test_paired_mates_use_independent_floors_and_one_warning_per_category(
    inference_module, make_config, make_sequali
):
    payload = inference_module.infer_short(
        make_sequali(
            read1_length=10,
            read1_bad=tuple(range(-10, 0)),
            read2_length=12,
            read2_bad=tuple(range(-12, 0)),
        ),
        make_config(tailcrop_min_retained_fraction=0.5),
        "short_pe",
    )

    assert payload["recommendations"] == {
        "read1": recommendation(5),
        "read2": recommendation(6),
    }
    assert payload["tool_args"]["fastp"] == {
        "--max_len1": 5,
        "--max_len2": 6,
        "--average_qual": 30,
        "--length_required": 2,
    }
    assert payload["warnings"] == [
        {
            "code": "tailcrop_safeguard_adjusted",
            "message": (
                "Retained-length safeguard adjusted a raw tailcrop cutoff for "
                "read1 and read2."
            ),
        },
    ]


def test_paired_headcrop_and_post_headcrop_maximum_are_mate_specific(
    inference_module, make_config, make_sequali
):
    payload = inference_module.infer_short(
        make_sequali(
            read1_bad=(-3, -2, -1),
            read2_bad=(),
            read1_imbalance=2,
        ),
        make_config(
            composition_headcrop_enabled=True,
            tailcrop_min_retained_fraction=0.0,
        ),
        "short_pe",
    )

    assert payload["recommendations"] == {
        "read1": recommendation(7, headcrop=3),
        "read2": recommendation(None),
    }
    assert payload["tool_args"]["fastp"] == {
        "--trim_front1": 3,
        "--max_len1": 4,
        "--average_qual": 30,
        "--length_required": 2,
    }
    assert payload["warnings"] == [
        {
            "code": "composition_headcrop_enabled",
            "message": (
                "Composition-derived 5' headcrop is explicitly enabled for read1 "
                "and read2."
            ),
        }
    ]


def test_paired_fastp_arguments_preserve_grouped_rendering_order(
    inference_module, make_config, make_sequali
):
    payload = inference_module.infer_short(
        make_sequali(
            read1_bad=(-3, -2, -1),
            read2_bad=(-2, -1),
            read1_imbalance=2,
            read2_imbalance=2,
        ),
        make_config(
            composition_headcrop_enabled=True,
            tailcrop_min_retained_fraction=0.0,
        ),
        "short_pe",
    )

    assert payload["recommendations"] == {
        "read1": recommendation(7, headcrop=3),
        "read2": recommendation(8, headcrop=3),
    }
    fastp_args = payload["tool_args"]["fastp"]
    assert list(fastp_args) == [
        "--trim_front1",
        "--trim_front2",
        "--max_len1",
        "--max_len2",
        "--average_qual",
        "--length_required",
    ]
    assert inference_module.render_tool_arguments(payload["tool_args"]) == (
        "fastp --trim_front1 3 --trim_front2 3 --max_len1 4 --max_len2 5 "
        "--average_qual 30 --length_required 2"
    )


def test_composition_headcrop_requires_opt_in(
    inference_module, make_config, make_sequali
):
    data = make_sequali(read1_imbalance=2)
    disabled = inference_module.infer_short(data, make_config(), "short_se")
    enabled = inference_module.infer_short(
        data, make_config(composition_headcrop_enabled=True), "short_se"
    )

    assert disabled["recommendations"]["read1"]["headcrop_bases"] == 0
    assert enabled["recommendations"]["read1"]["headcrop_bases"] == 3
    assert enabled["tool_args"]["fastp"]["--trim_front1"] == 3
    assert enabled["warnings"] == [
        {
            "code": "composition_headcrop_enabled",
            "message": "Composition-derived 5' headcrop is explicitly enabled for read1.",
        }
    ]


def test_malformed_sequali_json(tmp_path, inference_module, inference_schemas):
    path = tmp_path / "sequali.json"
    path.write_text("{", encoding="utf-8")

    with pytest.raises(json.JSONDecodeError):
        inference_module.load_json(path, inference_schemas["sequali"])


def test_missing_sequali_file_is_rejected(
    tmp_path, inference_module, inference_schemas
):
    with pytest.raises(ValueError, match="Sequali JSON does not exist"):
        inference_module.load_json(
            tmp_path / "missing.json", inference_schemas["sequali"]
        )


def test_complete_workflow_config_extracts_auto_preprocess(
    tmp_path, inference_module, inference_schemas, make_config
):
    path = tmp_path / "config.yaml"
    auto_preprocess = make_config()
    config = {
        "run_sheet": "config/run_sheet.tsv",
        "auto_preprocess": auto_preprocess,
        "fastp": {"extra": "--qualified_quality_phred 20"},
    }
    path.write_text(yaml.safe_dump(config), encoding="utf-8")
    assert (
        inference_module.load_config(path, inference_schemas["config"])
        == auto_preprocess
    )


def test_cli_rejects_long_data_type():
    result = subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            "--sequali-json",
            "unused.json",
            "--data-type",
            "long",
            "--config",
            "unused.yaml",
            "--config-schema",
            str(CONFIG_SCHEMA),
            "--sequali-schema",
            str(SEQUALI_SCHEMA),
            "--output",
            "unused-output.json",
        ],
        check=False,
        capture_output=True,
        text=True,
    )

    assert result.returncode != 0
    assert "invalid choice: 'long'" in result.stderr


def test_structured_warnings_have_stable_order_and_content(
    inference_module, make_config, make_sequali
):
    payload = inference_module.infer_short(
        make_sequali(
            read1_bad=tuple(range(-10, 0)),
            read1_imbalance=2,
        ),
        make_config(
            composition_headcrop_enabled=True,
            tailcrop_min_retained_fraction=1.0,
        ),
        "short_se",
    )

    assert payload["warnings"] == [
        {
            "code": "tailcrop_safeguard_adjusted",
            "message": (
                "Retained-length safeguard adjusted a raw tailcrop cutoff for read1."
            ),
        },
        {
            "code": "tailcrop_suppressed_as_no_op",
            "message": (
                "Retained-length safeguard suppressed the tailcrop as a no-op for "
                "read1."
            ),
        },
        {
            "code": "composition_headcrop_enabled",
            "message": "Composition-derived 5' headcrop is explicitly enabled for read1.",
        },
    ]


def test_cli_start_is_flushed_and_completion_is_printed(
    tmp_path, monkeypatch, inference_module, make_config, make_sequali
):
    sequali_path = tmp_path / "sequali.json"
    config_path = tmp_path / "config.yaml"
    output_path = tmp_path / "output.json"
    sequali_path.write_text(json.dumps(make_sequali()), encoding="utf-8")
    auto_preprocess = make_config()
    del auto_preprocess["mott_positive_weight"]
    del auto_preprocess["mott_negative_weight"]
    config_path.write_text(
        yaml.safe_dump(
            {
                "run_sheet": "config/run_sheet.tsv",
                "auto_preprocess": auto_preprocess,
            }
        ),
        encoding="utf-8",
    )
    monkeypatch.setattr(
        sys,
        "argv",
        [
            str(SCRIPT),
            "--sequali-json",
            str(sequali_path),
            "--data-type",
            "short_se",
            "--config",
            str(config_path),
            "--config-schema",
            str(CONFIG_SCHEMA),
            "--sequali-schema",
            str(SEQUALI_SCHEMA),
            "--output",
            str(output_path),
        ],
    )

    calls = []

    def record_print(*args, **kwargs):
        calls.append((args, kwargs))

    monkeypatch.setattr(builtins, "print", record_print)
    inference_module.main()

    assert calls[0] == (("Auto-preprocess inference started",), {"flush": True})
    assert calls[1] == (("Auto-preprocess inference complete",), {})
    assert calls[-1] == (("Warnings: none",), {})
    assert json.loads(output_path.read_text(encoding="utf-8"))["warnings"] == []


def test_cli_omitted_mott_weights_execute_with_unit_fallback(
    tmp_path, monkeypatch, inference_module, make_config, make_sequali
):
    sequali_path = tmp_path / "sequali.json"
    output_path = tmp_path / "output.json"
    sequali_path.write_text(json.dumps(make_sequali()), encoding="utf-8")
    auto_preprocess = make_config(average_quality=17, minimum_length=4)
    del auto_preprocess["mott_positive_weight"]
    del auto_preprocess["mott_negative_weight"]
    monkeypatch.setattr(
        sys,
        "argv",
        [
            str(SCRIPT),
            "--sequali-json",
            str(sequali_path),
            "--data-type",
            "short_se",
            "--config-json",
            json.dumps(auto_preprocess),
            "--config-schema",
            str(CONFIG_SCHEMA),
            "--sequali-schema",
            str(SEQUALI_SCHEMA),
            "--output",
            str(output_path),
        ],
    )

    inference_module.main()

    assert json.loads(output_path.read_text(encoding="utf-8"))["tool_args"][
        "fastp"
    ] == {
        "--average_qual": 17,
        "--length_required": 4,
    }


@pytest.mark.parametrize(
    ("unsupported_key", "unsupported_value"),
    [
        ("config", "legacy.yaml"),
        ("short_read", {}),
        ("crop_start_position", 30),
        ("min_crop_length", 30),
        ("unknown", True),
    ],
)
def test_cli_config_json_rejects_unsupported_auto_preprocess_keys(
    tmp_path,
    monkeypatch,
    inference_module,
    make_config,
    make_sequali,
    unsupported_key,
    unsupported_value,
):
    sequali_path = tmp_path / "sequali.json"
    output_path = tmp_path / "output.json"
    sequali_path.write_text(json.dumps(make_sequali()), encoding="utf-8")
    auto_preprocess = (
        {"short_read": make_config()}
        if unsupported_key == "short_read"
        else {**make_config(), unsupported_key: unsupported_value}
    )
    monkeypatch.setattr(
        sys,
        "argv",
        [
            str(SCRIPT),
            "--sequali-json",
            str(sequali_path),
            "--data-type",
            "short_se",
            "--config-json",
            json.dumps(auto_preprocess),
            "--config-schema",
            str(CONFIG_SCHEMA),
            "--sequali-schema",
            str(SEQUALI_SCHEMA),
            "--output",
            str(output_path),
        ],
    )

    with pytest.raises(
        ValueError,
        match="Inline auto-preprocess config failed schema validation",
    ):
        inference_module.main()


def test_cli_config_file_is_validated_with_workflow_schema(
    tmp_path, monkeypatch, inference_module, make_config, make_sequali
):
    sequali_path = tmp_path / "sequali.json"
    config_path = tmp_path / "config.yaml"
    output_path = tmp_path / "output.json"
    sequali_path.write_text(json.dumps(make_sequali()), encoding="utf-8")
    config_path.write_text(
        yaml.safe_dump(
            {
                "run_sheet": "config/run_sheet.tsv",
                "auto_preprocess": make_config(minimum_length=0),
            }
        ),
        encoding="utf-8",
    )
    monkeypatch.setattr(
        sys,
        "argv",
        [
            str(SCRIPT),
            "--sequali-json",
            str(sequali_path),
            "--data-type",
            "short_se",
            "--config",
            str(config_path),
            "--config-schema",
            str(CONFIG_SCHEMA),
            "--sequali-schema",
            str(SEQUALI_SCHEMA),
            "--output",
            str(output_path),
        ],
    )

    with pytest.raises(
        ValueError,
        match=r"Workflow config failed schema validation at "
        r"\$\.auto_preprocess\.minimum_length",
    ):
        inference_module.main()


@pytest.mark.parametrize("data_type", ["short_se", "short_pe"])
@pytest.mark.parametrize("config_mode", ["config-json", "config-file"])
def test_native_script_and_cli_modes_have_equivalent_outputs_and_logs(
    tmp_path, make_config, make_sequali, data_type, config_mode
):
    sequali = make_sequali(
        read1_bad=(-3, -2, -1),
        read2_bad=(-2, -1) if data_type == "short_pe" else None,
    )
    auto_preprocess_config = make_config(tailcrop_min_retained_fraction=0.0)

    native_result, native_output, native_log = run_native_script(
        tmp_path,
        data_type,
        sequali,
        auto_preprocess_config,
    )
    cli_result, cli_output, _ = run_cli(
        tmp_path,
        data_type,
        sequali,
        auto_preprocess_config,
        config_mode=config_mode,
    )

    assert native_result.returncode == 0, native_result.stdout + native_result.stderr
    assert cli_result.returncode == 0, cli_result.stdout + cli_result.stderr
    assert cli_result.stderr == ""
    assert native_output.read_bytes() == cli_output.read_bytes()
    assert native_log.read_text(encoding="utf-8") == cli_result.stdout
    assert cli_result.stdout.startswith("Auto-preprocess inference started\n")
    assert "Auto-preprocess inference complete\n" in cli_result.stdout


@pytest.mark.parametrize(
    ("case_name", "data_type", "sequali_update", "config_update", "diagnostic"),
    [
        (
            "invalid-sequali",
            "short_se",
            lambda data: data["meta"].update(sequali_version="unsupported"),
            lambda cfg: None,
            "Sequali JSON failed schema validation",
        ),
        (
            "missing-pe-read2",
            "short_pe",
            lambda data: None,
            lambda cfg: None,
            "complete read2 section set required for data_type=short_pe",
        ),
        (
            "nonfinite-sequali",
            "short_se",
            lambda data: data["per_position_mean_quality_and_spread"][
                "end_percentiles"
            ][5][1].__setitem__(-1, float("nan")),
            lambda cfg: None,
            "JSON contains unsupported non-finite number: NaN",
        ),
        (
            "invalid-config",
            "short_se",
            lambda data: None,
            lambda cfg: cfg.update(minimum_length=0),
            "Inline auto-preprocess config failed schema validation",
        ),
        (
            "nonfinite-config",
            "short_se",
            lambda data: None,
            lambda cfg: cfg.update(tailcrop_mott_quality_cutoff=float("inf")),
            "JSON contains unsupported non-finite number: Infinity",
        ),
        (
            "disabled-inference",
            "short_se",
            lambda data: None,
            lambda cfg: cfg.update(enabled=False),
            "auto_preprocess.enabled is false",
        ),
    ],
)
def test_native_and_cli_failures_have_equivalent_error_diagnostics(
    tmp_path,
    make_config,
    make_sequali,
    case_name,
    data_type,
    sequali_update,
    config_update,
    diagnostic,
):
    sequali = make_sequali()
    auto_preprocess_config = make_config()
    sequali_update(sequali)
    config_update(auto_preprocess_config)

    native_result, native_output, native_log = run_native_script(
        tmp_path,
        data_type,
        sequali,
        auto_preprocess_config,
    )
    cli_result, cli_output, _ = run_cli(
        tmp_path,
        data_type,
        sequali,
        auto_preprocess_config,
    )

    assert case_name
    assert native_result.returncode != 0
    assert cli_result.returncode != 0
    assert not native_output.exists()
    assert not cli_output.exists()
    native_diagnostic = native_log.read_text(encoding="utf-8")
    assert native_diagnostic.startswith("Auto-preprocess inference started\n")
    assert "ERROR:" in native_diagnostic
    assert diagnostic in native_diagnostic
    assert "ERROR:" in cli_result.stderr
    assert diagnostic in cli_result.stderr


@pytest.mark.parametrize(
    ("missing_argument", "diagnostic"),
    [
        ("sequali-json", "Sequali JSON does not exist"),
        ("config", "Config file does not exist"),
        ("config-schema", "Schema file does not exist"),
        ("sequali-schema", "Schema file does not exist"),
    ],
)
def test_cli_missing_input_assets_keep_error_contract(
    tmp_path,
    make_config,
    make_sequali,
    missing_argument,
    diagnostic,
):
    sequali = make_sequali()
    _, _, command = run_cli(
        tmp_path,
        "short_se",
        sequali,
        make_config(),
        config_mode="config-file",
    )
    argument = f"--{missing_argument}"
    value_index = command.index(argument) + 1
    command[value_index] = str(tmp_path / f"missing-{missing_argument}")

    result = subprocess.run(
        command,
        cwd=ROOT,
        check=False,
        capture_output=True,
        text=True,
    )

    assert result.returncode != 0
    assert result.stdout == "Auto-preprocess inference started\n"
    assert result.stderr.startswith("ERROR: ")
    assert diagnostic in result.stderr


@pytest.mark.parametrize(
    ("schema_text", "diagnostic"),
    [
        ("{", "Schema YAML parse failed"),
        ("type: not-a-json-schema-type\n", "Invalid JSON Schema"),
    ],
)
def test_cli_invalid_schema_keeps_error_contract(
    tmp_path, make_config, make_sequali, schema_text, diagnostic
):
    _, _, command = run_cli(
        tmp_path,
        "short_se",
        make_sequali(),
        make_config(),
    )
    invalid_schema = tmp_path / "invalid.schema.yaml"
    invalid_schema.write_text(schema_text, encoding="utf-8")
    value_index = command.index("--config-schema") + 1
    command[value_index] = str(invalid_schema)

    result = subprocess.run(
        command,
        cwd=ROOT,
        check=False,
        capture_output=True,
        text=True,
    )

    assert result.returncode != 0
    assert result.stdout == "Auto-preprocess inference started\n"
    assert result.stderr.startswith("ERROR: ")
    assert diagnostic in result.stderr


def test_cli_help_retains_standalone_interface():
    result = subprocess.run(
        [sys.executable, str(SCRIPT), "--help"],
        cwd=ROOT,
        check=False,
        capture_output=True,
        text=True,
    )

    assert result.returncode == 0
    assert result.stderr == ""
    for option in [
        "--sequali-json",
        "--data-type {short_se,short_pe}",
        "--config CONFIG",
        "--config-json CONFIG_JSON",
        "--config-schema",
        "--sequali-schema",
        "--output",
    ]:
        assert option in result.stdout


def test_native_config_validation_rejects_nested_nonfinite_numbers(
    inference_module, inference_schemas, make_config
):
    config = make_config()
    config["nested"] = {"value": float("-inf")}

    with pytest.raises(
        ValueError,
        match="JSON contains unsupported non-finite number: -Infinity",
    ):
        inference_module.load_native_config(config, inference_schemas["config"])


def test_unexpected_native_exception_writes_traceback_to_log(
    tmp_path, monkeypatch, inference_module, make_config
):
    log_path = tmp_path / "native.log"
    snakemake_object = SimpleNamespace(
        input=SimpleNamespace(
            sequali_json="sequali.json",
            config_schema=CONFIG_SCHEMA,
            sequali_schema=SEQUALI_SCHEMA,
        ),
        output=SimpleNamespace(auto_params="output.json"),
        params=SimpleNamespace(
            data_type="short_se",
            auto_preprocess_config=make_config(),
        ),
        log=[log_path],
    )

    def raise_unexpected(**kwargs):
        raise RuntimeError("unexpected native failure")

    monkeypatch.setattr(inference_module, "run_inference", raise_unexpected)

    with pytest.raises(RuntimeError, match="unexpected native failure"):
        inference_module.snakemake_main(snakemake_object)

    diagnostic = log_path.read_text(encoding="utf-8")
    assert "Traceback (most recent call last):" in diagnostic
    assert "RuntimeError: unexpected native failure" in diagnostic
