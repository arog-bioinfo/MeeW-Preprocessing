"""Pure run-sheet parsing, normalization, and validation setup helpers."""

from __future__ import annotations

import csv
import os
from collections.abc import Callable, Mapping
from pathlib import Path
from typing import Any

import yaml
from jsonschema import Draft202012Validator

REQUIRED_HEADERS = frozenset({"group", "sample", "data_type", "read1"})
PATH_FIELDS = ("read1", "read2", "contam")
MAX_NORMALIZED_BASENAME_LENGTH = 220


class RunSheetError(ValueError):
    """Report a run-sheet structure or normalized naming error."""


def expand_configured_path(value: object) -> str:
    """Expand environment variables and ``~`` without rebasing relative paths."""
    if value is None or str(value).strip() == "":
        return ""
    return os.path.expanduser(os.path.expandvars(str(value)))


def build_samples_row_validator(schema_path: str | Path) -> Draft202012Validator:
    """Construct the validator for the current, self-contained samples schema."""
    schema = yaml.safe_load(Path(schema_path).read_text(encoding="utf-8"))
    Draft202012Validator.check_schema(schema)
    return Draft202012Validator(schema)


def load_run_sheet(
    run_sheet: str | Path,
    validate_row: Callable[[Mapping[str, str]], Any] | None,
) -> list[dict[str, str]]:
    """Parse and normalize a TSV run sheet using one caller-owned validator."""
    run_sheet_path = Path(run_sheet)

    if not run_sheet_path.exists():
        raise RunSheetError(f"Run sheet does not exist: {run_sheet}")

    try:
        with run_sheet_path.open(newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle, delimiter="\t", strict=True)
            fieldnames = [field.strip() for field in (reader.fieldnames or [])]
            if not fieldnames or any(not field for field in fieldnames):
                raise RunSheetError(
                    f"Run sheet must have a non-empty TSV header: {run_sheet}"
                )

            duplicate_headers = sorted(
                {field for field in fieldnames if fieldnames.count(field) > 1}
            )
            if duplicate_headers:
                raise RunSheetError(
                    "Run sheet contains duplicate header(s): "
                    + ", ".join(duplicate_headers)
                )
            reader.fieldnames = fieldnames

            missing_headers = sorted(REQUIRED_HEADERS - set(fieldnames))
            if missing_headers:
                raise RunSheetError(
                    "Run sheet is missing required header(s): "
                    + ", ".join(missing_headers)
                )

            rows = []
            for row in reader:
                if None in row:
                    raise RunSheetError(
                        f"Run sheet row {reader.line_num} has more fields than the header"
                    )
                if any(value is None for value in row.values()):
                    raise RunSheetError(
                        f"Run sheet row {reader.line_num} has fewer fields than the header"
                    )

                raw_row = {key: value.strip() for key, value in row.items()}
                if validate_row is not None:
                    validate_row(raw_row)

                normalized = dict(raw_row)
                for path_field in PATH_FIELDS:
                    if normalized.get(path_field):
                        normalized[path_field] = expand_configured_path(
                            normalized[path_field]
                        )
                rows.append(normalized)
    except csv.Error as error:
        raise RunSheetError(f"Malformed TSV run sheet {run_sheet}: {error}") from error

    if not rows:
        raise RunSheetError(
            f"Run sheet must contain at least one data row: {run_sheet}"
        )

    explicit_by_group: dict[str, set[str]] = {}
    for row in rows:
        if row["sample"]:
            explicit_by_group.setdefault(row["group"], set()).add(row["sample"])

    counters: dict[tuple[str, str], int] = {}
    used_by_group: dict[str, set[str]] = {}
    normalized_rows = []
    for row in rows:
        normalized = dict(row)
        group = normalized["group"]
        sample = normalized["sample"]
        used = used_by_group.setdefault(group, set())
        if not sample:
            counter_key = (group, normalized["data_type"])
            number = counters.get(counter_key, 0) + 1
            reserved = explicit_by_group.get(group, set())
            while (
                f"{normalized['data_type']}_{number}" in reserved
                or f"{normalized['data_type']}_{number}" in used
            ):
                number += 1
            counters[counter_key] = number
            sample = f"{normalized['data_type']}_{number}"
            normalized["sample"] = sample

        if sample in used:
            raise RunSheetError(
                "Run sheet contains duplicate normalized (group, sample) pair: "
                f"({group!r}, {sample!r})"
            )
        used.add(sample)

        basename = f"{group}_{sample}"
        if len(basename) > MAX_NORMALIZED_BASENAME_LENGTH:
            raise RunSheetError(
                "Normalized row basename exceeds 220 characters for "
                f"group={group!r}, sample={sample!r}: {len(basename)} characters. "
                "Shorten the group or sample identifier."
            )
        normalized_rows.append(normalized)

    return normalized_rows
