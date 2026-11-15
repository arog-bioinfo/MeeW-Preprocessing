from copy import deepcopy
from pathlib import Path

import pytest
import yaml
from jsonschema import Draft202012Validator, ValidationError

ROOT = Path(__file__).resolve().parents[1]
CONFIG_SCHEMA_PATH = ROOT / "workflow" / "schemas" / "config.schema.yaml"
SAMPLES_SCHEMA_PATH = ROOT / "workflow" / "schemas" / "samples.schema.yaml"
SEQUALI_SCHEMA_PATH = ROOT / "workflow" / "schemas" / "sequali_short_read.schema.yaml"
OUTPUT_SCHEMA_PATH = (
    ROOT / "workflow" / "schemas" / "auto_preprocess_output.schema.yaml"
)
REQUIRED_INFERENCE_KEYS = [
    "composition_headcrop_enabled",
    "max_headcrop",
    "headcrop_base_imbalance_threshold",
    "tailcrop_mott_quality_cutoff",
    "tailcrop_min_retained_fraction",
    "average_quality",
    "minimum_length",
]
OPTIONAL_INFERENCE_KEYS = [
    "mott_positive_weight",
    "mott_negative_weight",
]
INFERENCE_KEYS = REQUIRED_INFERENCE_KEYS + OPTIONAL_INFERENCE_KEYS
LEGACY_INFERENCE_KEYS = [
    "tailcrop_min_consecutive_good_bases",
    "tailcrop_mean_quality_threshold",
    "tailcrop_min_bad_suffix_span_bases",
    "head_base_imbalance_threshold",
    "quality_crop_mean_threshold",
    "quality_crop_min_bad_suffix_span_bases",
    "quality_crop_min_retained_fraction",
    "quality_crop_low_quality_bucket_boundary",
    "quality_crop_low_quality_fraction_threshold",
]
LEGACY_RECOMMENDATION_KEYS = [
    "front_trim_bases",
    "quality_crop_base_position",
    "maximum_length_after_front_trim_bases",
    "headcrop_bases",
    "tailcrop_base_position",
    "maximum_length_after_headcrop_bases",
]
LEGACY_WARNING_CODES = [
    "tailcrop_good_run_not_found",
    "quality_crop_safeguard_adjusted",
    "quality_crop_suppressed_as_no_op",
]


def load_yaml(path):
    return yaml.safe_load(path.read_text(encoding="utf-8"))


CONFIG_SCHEMA = load_yaml(CONFIG_SCHEMA_PATH)
SAMPLES_SCHEMA = load_yaml(SAMPLES_SCHEMA_PATH)
SEQUALI_SCHEMA = load_yaml(SEQUALI_SCHEMA_PATH)
OUTPUT_SCHEMA = load_yaml(OUTPUT_SCHEMA_PATH)
for schema in (CONFIG_SCHEMA, SAMPLES_SCHEMA, SEQUALI_SCHEMA, OUTPUT_SCHEMA):
    Draft202012Validator.check_schema(schema)
CONFIG_VALIDATOR = Draft202012Validator(CONFIG_SCHEMA)
SAMPLES_VALIDATOR = Draft202012Validator(SAMPLES_SCHEMA)
SEQUALI_VALIDATOR = Draft202012Validator(SEQUALI_SCHEMA)
OUTPUT_VALIDATOR = Draft202012Validator(OUTPUT_SCHEMA)
BASE_WORKFLOW_CONFIG = load_yaml(ROOT / "config" / "config.yaml")


def fresh_workflow_config():
    return deepcopy(BASE_WORKFLOW_CONFIG)


@pytest.mark.parametrize(
    "config_path",
    [ROOT / "config" / "config.yaml"]
    + sorted((ROOT / ".test" / "config").glob("config.*.yaml")),
    ids=lambda path: path.name,
)
def test_workflow_configs_match_schema(config_path):
    CONFIG_VALIDATOR.validate(load_yaml(config_path))


@pytest.mark.parametrize(
    "legacy_key",
    ["config", "short_read", "crop_start_position", "min_crop_length"],
)
def test_workflow_schema_rejects_legacy_auto_preprocess_keys(legacy_key):
    config = fresh_workflow_config()
    if legacy_key == "short_read":
        config["auto_preprocess"] = {
            "enabled": True,
            "short_read": {
                key: config["auto_preprocess"][key] for key in INFERENCE_KEYS
            },
        }
    else:
        config["auto_preprocess"][legacy_key] = "legacy.yaml"

    with pytest.raises(ValidationError):
        CONFIG_VALIDATOR.validate(config)


def test_workflow_schema_requires_auto_preprocess_root():
    config = fresh_workflow_config()
    del config["auto_preprocess"]

    with pytest.raises(ValidationError):
        CONFIG_VALIDATOR.validate(config)


def test_workflow_schema_rejects_removed_stage_sheets_setting():
    config = fresh_workflow_config()
    config["stage_sheets"] = {"enabled": False}

    with pytest.raises(ValidationError):
        CONFIG_VALIDATOR.validate(config)


@pytest.mark.parametrize("field", ["group", "sample"])
def test_run_sheet_identifier_components_allow_128_safe_characters(field):
    row = {
        "group": "group",
        "sample": "sample",
        "data_type": "short_se",
        "read1": "reads.fastq",
    }
    row[field] = "A" + "b" * 127

    SAMPLES_VALIDATOR.validate(row)


@pytest.mark.parametrize("field", ["group", "sample"])
@pytest.mark.parametrize("value", ["A" + "b" * 128, "bad/name", ".hidden"])
def test_run_sheet_identifier_components_reject_unsafe_or_oversized_values(
    field, value
):
    row = {
        "group": "group",
        "sample": "sample",
        "data_type": "short_se",
        "read1": "reads.fastq",
    }
    row[field] = value

    with pytest.raises(ValidationError):
        SAMPLES_VALIDATOR.validate(row)


def test_enabled_auto_preprocess_contract_has_exactly_nine_scientific_settings():
    config = fresh_workflow_config()

    assert len(INFERENCE_KEYS) == 9
    assert set(config["auto_preprocess"]) == {"enabled", *INFERENCE_KEYS}
    assert config["auto_preprocess"]["tailcrop_mott_quality_cutoff"] == 20
    assert config["auto_preprocess"]["mott_positive_weight"] == 1.0
    assert config["auto_preprocess"]["mott_negative_weight"] == 1.0


@pytest.mark.parametrize(
    "missing_key",
    REQUIRED_INFERENCE_KEYS,
)
@pytest.mark.parametrize("enabled", [True, None], ids=["enabled", "default-enabled"])
def test_workflow_schema_requires_inference_fields_unless_disabled(
    missing_key, enabled
):
    config = fresh_workflow_config()
    del config["auto_preprocess"][missing_key]
    if enabled is None:
        del config["auto_preprocess"]["enabled"]
    else:
        config["auto_preprocess"]["enabled"] = enabled

    with pytest.raises(ValidationError):
        CONFIG_VALIDATOR.validate(config)


def test_workflow_schema_allows_disabled_inference_settings_to_be_omitted():
    config = fresh_workflow_config()
    config["auto_preprocess"] = {"enabled": False}

    CONFIG_VALIDATOR.validate(config)


@pytest.mark.parametrize("legacy_key", LEGACY_INFERENCE_KEYS)
@pytest.mark.parametrize("disabled", [False, True])
def test_workflow_schema_rejects_legacy_inference_fields(legacy_key, disabled):
    config = fresh_workflow_config()
    if disabled:
        config["auto_preprocess"] = {"enabled": False}
    config["auto_preprocess"][legacy_key] = 20

    with pytest.raises(ValidationError) as error:
        CONFIG_VALIDATOR.validate(config)
    assert error.value.validator == "additionalProperties"


@pytest.mark.parametrize(
    ("key", "invalid_value"),
    [
        ("composition_headcrop_enabled", "false"),
        ("max_headcrop", 1.5),
        ("max_headcrop", -1),
        ("headcrop_base_imbalance_threshold", "0.1"),
        ("headcrop_base_imbalance_threshold", -0.1),
        ("headcrop_base_imbalance_threshold", 1.1),
        ("tailcrop_mott_quality_cutoff", "20"),
        ("tailcrop_mott_quality_cutoff", -1),
        ("tailcrop_mott_quality_cutoff", 93.1),
        ("mott_positive_weight", "1.0"),
        ("mott_positive_weight", -0.1),
        ("mott_negative_weight", "1.0"),
        ("mott_negative_weight", -0.1),
        ("tailcrop_min_retained_fraction", "0.75"),
        ("tailcrop_min_retained_fraction", -0.1),
        ("tailcrop_min_retained_fraction", 1.1),
        ("average_quality", "30"),
        ("average_quality", -1),
        ("average_quality", 94),
        ("minimum_length", "30"),
        ("minimum_length", 0),
    ],
)
def test_workflow_schema_rejects_invalid_inference_fields(key, invalid_value):
    config = fresh_workflow_config()
    config["auto_preprocess"]["enabled"] = False
    config["auto_preprocess"][key] = invalid_value

    with pytest.raises(ValidationError):
        CONFIG_VALIDATOR.validate(config)


@pytest.mark.parametrize("quality_cutoff", [0, 93])
def test_workflow_schema_accepts_mott_quality_cutoff_boundaries(quality_cutoff):
    config = fresh_workflow_config()
    config["auto_preprocess"]["tailcrop_mott_quality_cutoff"] = quality_cutoff

    CONFIG_VALIDATOR.validate(config)


@pytest.mark.parametrize("enabled", [True, None], ids=["enabled", "default-enabled"])
def test_workflow_schema_allows_mott_weights_to_be_omitted(enabled):
    config = fresh_workflow_config()
    for key in OPTIONAL_INFERENCE_KEYS:
        del config["auto_preprocess"][key]
    if enabled is None:
        del config["auto_preprocess"]["enabled"]
    else:
        config["auto_preprocess"]["enabled"] = enabled

    CONFIG_VALIDATOR.validate(config)


@pytest.mark.parametrize("key", OPTIONAL_INFERENCE_KEYS)
@pytest.mark.parametrize("value", [0, 0.5, 1.0, 100])
def test_workflow_schema_accepts_nonnegative_mott_weights(key, value):
    config = fresh_workflow_config()
    config["auto_preprocess"][key] = value

    CONFIG_VALIDATOR.validate(config)


@pytest.mark.parametrize("paired", [False, True], ids=["single-end", "paired-end"])
def test_sequali_short_read_reports_match_schema(make_sequali, paired):
    data = make_sequali(read2_bad=() if paired else None)
    data["unrelated_report_section"] = {"producer_owned": True}
    assert "per_position_quality_distribution" not in data
    assert "per_position_quality_distribution_read2" not in data

    SEQUALI_VALIDATOR.validate(data)


@pytest.mark.parametrize(
    "mutation",
    [
        "wrong_version",
        "missing_r1_section",
        "nonpositive_maximum_length",
        "malformed_percentile_pair",
        "missing_mean",
        "duplicate_mean",
        "short_mean",
        "mean_out_of_range",
        "missing_composition_series",
        "invalid_composition_label",
        "composition_fraction_out_of_range",
    ],
)
def test_sequali_schema_rejects_malformed_r1_structures(make_sequali, mutation):
    data = make_sequali()
    quality = data["per_position_mean_quality_and_spread"]
    composition = data["per_position_base_content"]

    if mutation == "wrong_version":
        data["meta"]["sequali_version"] = "1.1.0"
    elif mutation == "missing_r1_section":
        del data["per_position_mean_quality_and_spread"]
    elif mutation == "nonpositive_maximum_length":
        data["summary"]["maximum_length"] = 0
    elif mutation == "malformed_percentile_pair":
        quality["end_percentiles"][0] = ["bottom 1%"]
    elif mutation == "missing_mean":
        quality["end_percentiles"][5][0] = "average"
    elif mutation == "duplicate_mean":
        quality["end_percentiles"].append(
            ["mean", list(quality["end_percentiles"][5][1])]
        )
    elif mutation == "short_mean":
        quality["end_percentiles"][5][1].pop()
    elif mutation == "mean_out_of_range":
        quality["end_percentiles"][5][1][-1] = 94
    elif mutation == "missing_composition_series":
        del composition["A"]
    elif mutation == "invalid_composition_label":
        composition["x_labels"][0] = "0"
    elif mutation == "composition_fraction_out_of_range":
        composition["A"][0] = 1.1

    with pytest.raises(ValidationError):
        SEQUALI_VALIDATOR.validate(data)


@pytest.mark.parametrize(
    "missing_section",
    [
        "summary_read2",
        "per_position_mean_quality_and_spread_read2",
        "per_position_base_content_read2",
    ],
)
def test_sequali_schema_requires_a_complete_r2_section_set(
    make_sequali, missing_section
):
    data = make_sequali(read2_bad=())
    del data[missing_section]

    with pytest.raises(ValidationError):
        SEQUALI_VALIDATOR.validate(data)


def test_sequali_schema_does_not_constrain_unused_percentile_values(make_sequali):
    data = make_sequali()
    mean = data["per_position_mean_quality_and_spread"]["end_percentiles"][5]
    data["per_position_mean_quality_and_spread"]["end_percentiles"] = [
        ["unused producer series", {"shape": "not consumed"}],
        mean,
    ]

    SEQUALI_VALIDATOR.validate(data)


@pytest.mark.parametrize("data_type", ["short_se", "short_pe"])
def test_generated_short_payloads_match_output_schema(
    inference_module, make_config, make_sequali, data_type
):
    data = (
        make_sequali(read1_bad=(-3, -2, -1), read2_bad=(-2, -1))
        if data_type == "short_pe"
        else make_sequali(read1_bad=(-3, -2, -1))
    )
    payload = inference_module.infer_short(data, make_config(), data_type)

    OUTPUT_VALIDATOR.validate(payload)
    assert (payload["recommendations"]["read2"] is not None) == (
        data_type == "short_pe"
    )


@pytest.mark.parametrize(
    "required_key",
    ["read1", "read2"],
)
def test_output_schema_requires_both_mate_recommendation_keys(
    inference_module, make_config, make_sequali, required_key
):
    payload = inference_module.infer_short(make_sequali(), make_config(), "short_se")
    del payload["recommendations"][required_key]

    with pytest.raises(ValidationError) as error:
        OUTPUT_VALIDATOR.validate(payload)
    assert error.value.validator == "required"


@pytest.mark.parametrize(
    "required_key",
    [
        "headcrop_bases",
        "tailcrop_base_position",
        "maximum_length_after_headcrop_bases",
    ],
)
def test_output_schema_requires_each_mate_recommendation_field(
    inference_module, make_config, make_sequali, required_key
):
    payload = inference_module.infer_short(make_sequali(), make_config(), "short_se")
    del payload["recommendations"]["read1"][required_key]

    with pytest.raises(ValidationError) as error:
        OUTPUT_VALIDATOR.validate(payload)
    assert error.value.validator == "required"


@pytest.mark.parametrize("legacy_key", LEGACY_RECOMMENDATION_KEYS)
def test_output_schema_rejects_legacy_recommendation_fields(
    inference_module, make_config, make_sequali, legacy_key
):
    payload = inference_module.infer_short(make_sequali(), make_config(), "short_se")
    payload["recommendations"][legacy_key] = 1

    with pytest.raises(ValidationError) as error:
        OUTPUT_VALIDATOR.validate(payload)
    assert error.value.validator == "additionalProperties"


@pytest.mark.parametrize(
    "legacy_key",
    [
        "front_trim_bases",
        "quality_crop_base_position",
        "maximum_length_after_front_trim_bases",
    ],
)
def test_output_schema_rejects_legacy_fields_inside_mate_recommendation(
    inference_module, make_config, make_sequali, legacy_key
):
    payload = inference_module.infer_short(make_sequali(), make_config(), "short_se")
    payload["recommendations"]["read1"][legacy_key] = 1

    with pytest.raises(ValidationError) as error:
        OUTPUT_VALIDATOR.validate(payload)
    assert error.value.validator == "additionalProperties"


@pytest.mark.parametrize("legacy_code", LEGACY_WARNING_CODES)
def test_output_schema_rejects_legacy_warning_codes(
    inference_module, make_config, make_sequali, legacy_code
):
    payload = inference_module.infer_short(
        make_sequali(read1_bad=(-4, -3, -2, -1)),
        make_config(tailcrop_min_retained_fraction=0.8),
        "short_se",
    )
    payload["warnings"][0]["code"] = legacy_code

    with pytest.raises(ValidationError) as error:
        OUTPUT_VALIDATOR.validate(payload)
    assert error.value.validator == "enum"


@pytest.mark.parametrize(
    ("retained_fraction", "expected_codes"),
    [
        (0.8, ["tailcrop_safeguard_adjusted"]),
        (
            1.0,
            ["tailcrop_safeguard_adjusted", "tailcrop_suppressed_as_no_op"],
        ),
    ],
)
def test_output_schema_accepts_current_warning_codes(
    inference_module,
    make_config,
    make_sequali,
    retained_fraction,
    expected_codes,
):
    payload = inference_module.infer_short(
        make_sequali(read1_bad=(-4, -3, -2, -1)),
        make_config(tailcrop_min_retained_fraction=retained_fraction),
        "short_se",
    )

    OUTPUT_VALIDATOR.validate(payload)
    assert [warning["code"] for warning in payload["warnings"]] == expected_codes


@pytest.mark.parametrize(
    "current_code",
    ["composition_headcrop_enabled"],
)
def test_output_schema_accepts_each_current_warning_code(
    inference_module, make_config, make_sequali, current_code
):
    payload = inference_module.infer_short(make_sequali(), make_config(), "short_se")
    payload["warnings"] = [{"code": current_code, "message": "Current warning."}]

    OUTPUT_VALIDATOR.validate(payload)
