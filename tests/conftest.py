from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
INFERENCE_SCRIPT = ROOT / "workflow" / "scripts" / "infer_preprocess_params.py"


def _load_module(module_name, path):
    spec = importlib.util.spec_from_file_location(module_name, path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = module
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="session")
def inference_module():
    return _load_module("infer_preprocess_params", INFERENCE_SCRIPT)


@pytest.fixture(scope="session")
def inference_schemas(inference_module):
    schema_root = ROOT / "workflow" / "schemas"
    return {
        "config": inference_module.load_schema(schema_root / "config.schema.yaml"),
        "sequali": inference_module.load_schema(
            schema_root / "sequali_short_read.schema.yaml"
        ),
    }


@pytest.fixture(scope="session")
def workflow_assertions_module():
    return _load_module(
        "assert_workflow_case_outputs",
        ROOT / "tests" / "assert_workflow_case_outputs.py",
    )


@pytest.fixture
def make_config():
    def factory(**overrides):
        auto_preprocess = {
            "composition_headcrop_enabled": False,
            "max_headcrop": 3,
            "headcrop_base_imbalance_threshold": 0.1,
            "tailcrop_mott_quality_cutoff": 20,
            "mott_positive_weight": 1.0,
            "mott_negative_weight": 1.0,
            "tailcrop_min_retained_fraction": 0.5,
            "average_quality": 30,
            "minimum_length": 2,
        }
        auto_preprocess.update(overrides)
        return auto_preprocess

    return factory


@pytest.fixture
def make_sequali():
    percentile_names = (
        "bottom 1%",
        "bottom 5%",
        "bottom 10%",
        "bottom 25%",
        "bottom 50%",
        "mean",
        "top 50%",
        "top 25%",
        "top 10%",
        "top 5%",
        "top 1%",
    )

    def mate_sections(
        nominal_length,
        bad_offsets,
        imbalance_index=None,
        mean_overrides=None,
    ):
        slot_count = 100
        observed_count = min(nominal_length, slot_count)
        first_observed = slot_count - observed_count
        means = [0.0] * first_observed + [30.0] * observed_count

        def slot_index(offset):
            assert -slot_count <= offset <= -1
            index = slot_count + offset
            assert index >= first_observed
            return index

        for offset in bad_offsets:
            index = slot_index(offset)
            means[index] = 10.0
        for offset, value in (mean_overrides or {}).items():
            means[slot_index(offset)] = value

        labels = [str(index) for index in range(1, 11)]
        bases = {
            "x_labels": labels,
            "A": [0.25] * len(labels),
            "C": [0.25] * len(labels),
            "G": [0.25] * len(labels),
            "T": [0.25] * len(labels),
        }
        if imbalance_index is not None:
            bases["A"][imbalance_index] = 0.45

        return {
            "quality": {
                # Deliberately bad whole-read bins must not drive tail inference.
                "x_labels": ["1-10"],
                "percentiles": [["mean", [1.0]]],
                "end_percentiles": [[name, list(means)] for name in percentile_names],
            },
            "bases": bases,
        }

    def factory(
        *,
        read1_bad=(),
        read2_bad=None,
        read1_length=10,
        read2_length=10,
        read1_imbalance=None,
        read2_imbalance=None,
        read1_mean_overrides=None,
        read2_mean_overrides=None,
    ):
        read1 = mate_sections(
            read1_length,
            read1_bad,
            read1_imbalance,
            read1_mean_overrides,
        )
        data = {
            "meta": {"sequali_version": "1.0.2"},
            "summary": {"maximum_length": read1_length},
            "per_position_mean_quality_and_spread": read1["quality"],
            "per_position_base_content": read1["bases"],
        }

        if read2_bad is not None:
            read2 = mate_sections(
                read2_length,
                read2_bad,
                read2_imbalance,
                read2_mean_overrides,
            )
            data.update(
                {
                    "summary_read2": {"maximum_length": read2_length},
                    "per_position_mean_quality_and_spread_read2": read2["quality"],
                    "per_position_base_content_read2": read2["bases"],
                }
            )
        return data

    return factory
