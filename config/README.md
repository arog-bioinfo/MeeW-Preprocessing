## Workflow configuration

The workflow processes any positive number of independent rows across any number of groups.

`config/config.yaml` is organized into sectioned blocks for workflow inputs,
automatic parameter inference, and tool parameters. Set these fields:

- `run_sheet`: path to a TSV file containing one or more preprocessing rows
- `auto_preprocess.enabled`: toggle workflow-integrated auto-inference of short-read preprocess args from raw Sequali JSON (default `true`); every short row is eligible independently, while long rows ignore this setting
- `auto_preprocess.<setting>`: inline short-read inference and configured fastp filtering settings; the seven baseline scientific settings are required unless `auto_preprocess.enabled` is explicitly `false`, while the two Mott weights are optional and default to `1.0`
- `sequali.extra`: optional extra options string passed to the sequali wrapper
- `fastp.adapters`: optional adapter-related options string passed to the fastp wrapper
- `fastp.extra`: optional extra options string passed to the fastp wrapper
- `chopper.extra`: optional manual extra options string passed to the chopper wrapper; Chopper does not consume workflow-inferred arguments
- `porechop_abi.extra`: optional extra options string passed to porechop_abi for generic `long` and Nanopore reads
- `porechop_abi.ab_initio_config`: optional path to an ab initio config file for porechop_abi
- `porechop_abi.custom_adapters`: optional path to a custom adapters file for porechop_abi (`--custom_adapters`, 3-line format: name/start/end)

All configured paths and run-sheet path fields expand `~` and environment
variables. Relative paths are intentionally left relative to the Snakemake
working directory. For example, paths in the consolidated fixture sheets
`.test/config/run.inference.tsv` and `.test/config/run.no_inference.tsv` remain
relative to `.test` when using `--directory .test`.

The tracked workflow fixtures use exactly two suites. Both sheets have the
superset header `group`, `sample`, `data_type`, `read1`, `read2`, `contam` and
contain all rows for their suite in one file. `config.inference.yaml` covers
`test_infer_pe` and the short-plus-long `test_hybrid` group with inference
enabled for short rows only. `config.no_inference.yaml` covers `test_short_se`,
`test_short_pe`, `test_long`, `test_nano_raw`, and `test_pacbio_raw` with
inference disabled. Both fixture configs share their Sequali, fastp, Chopper,
and Porechop_ABI settings.

### Run sheet format (TSV)

Required columns:

- `group`: nonempty result-root identifier
- `sample`: required header whose cell may be blank; blank cells receive deterministic generated identifiers
- `data_type`: one of `short_se`, `short_pe`, `long`, `nano-raw`, `nano-corr`, `nano-hq`, `pacbio-raw`, `pacbio-corr`, `pacbio-hifi`
- `read1`: path to input read file

Optional columns:

- `read2`: second read file path (required only for `short_pe`)
- `contam`: contaminant reference supported by long-read preprocessing

`group` and every nonblank `sample` use `^[A-Za-z0-9][A-Za-z0-9._-]{0,127}$`. A blank sample becomes `<data_type>_<n>` in original row order, with a separate counter for each group and data type. All explicit sample values in a group are reserved before generation, so generated values skip them. Normalized `(group, sample)` pairs must be unique, and the resulting `<group>_<sample>` basename may not exceed 220 characters. Group-less sheets and sheets missing the `sample` header are rejected without inference.

### Mode-specific validation rules

- `short_se`: requires `read1`, forbids `read2`
- `short_pe`: requires `read1` and `read2`
- `short_se` and `short_pe`: forbid `contam`; contaminant filtering is supported only by the long-read preprocessing path
- `long` and Flye-compatible long-read values require `read1` and forbid `read2`
- groups may contain any valid mix or repetition of these row types; there are no standalone group-level hybrid/cardinality restrictions

For long-read preprocessing, generic `long` and the Nanopore types (`nano-raw`, `nano-corr`, `nano-hq`) use raw FASTQ -> porechop_abi -> Chopper. Explicit PacBio types (`pacbio-raw`, `pacbio-corr`, `pacbio-hifi`) bypass porechop_abi and use raw FASTQ -> Chopper. The optional `contam` input enables Chopper contaminant filtering on either route and is available for every long-read type.

Each `results/{group}` root contains `reads/`, `qc/`, and `logs/`. Every row-specific filename starts with the exact `<group>_<sample>` basename. Final FASTQs and optional auto-parameter JSON files are placed in `reads/`; raw and processed HTML/JSON reports are under `qc/`; text logs are separated by tool under `logs/`. Aggregate targets follow the producing rules' resolved paths, including module pathvar overrides. Existing logs are not moved, and cached jobs do not automatically rerun solely to recreate them at the new paths.

### Example files

`config/config.yaml`:

```yaml
# ====================
# Workflow Inputs
# ====================

# --------------------
# Run Sheet
# --------------------
run_sheet: config/run_sheet.tsv

# ====================
# Automatic Parameter Inference
# ====================

# --------------------
# Auto Preprocess
# --------------------
auto_preprocess:
  enabled: true
  # Composition-derived 5' headcrop is disabled by default. When explicitly
  # enabled, the legacy A/T and G/C imbalance scan and maximum remain in use.
  composition_headcrop_enabled: false
  max_headcrop: 12
  headcrop_base_imbalance_threshold: 0.10

  # Infer each mate's 3' tailcrop from its Sequali end-anchor mean qualities.
  # Starting at the 3' end, choose the shortest suffix with the greatest positive
  # cumulative weighted sum of (tailcrop_mott_quality_cutoff - mean_quality).
  tailcrop_mott_quality_cutoff: 20
  mott_positive_weight: 1.0
  mott_negative_weight: 1.0
  # Safeguard a raw tailcrop trigger by retaining at least this fraction of
  # the corresponding nominal maximum raw-read length (rounded up).
  tailcrop_min_retained_fraction: 0.75

  # Short-read read-level filtering
  # Sequali/NanoPlot average quality uses expected-error rates,
  # while FastQC/Falco/Fastp average phred scores; Q>=X counts
  # can therefore differ substantially and are not strictly comparable.
  average_quality: 30
  minimum_length: 30

# ====================
# Tool Parameters
# ====================

# --------------------
# Quality Control
# --------------------
sequali:
  extra: ""

# --------------------
# Short-Read Preprocessing
# --------------------
fastp:
  adapters: ""
  extra: ""

# --------------------
# Long-Read Preprocessing
# --------------------
chopper:
  extra: ""

# --------------------
# Long-Read Adapter Trimming
# --------------------
porechop_abi:
  extra: ""
  ab_initio_config: "config/porechop_abi.config"
  custom_adapters: ""
```

`config/run_sheet.tsv` (short PE example):

```tsv
group	sample	data_type	read1	read2
my_short_pe	short	short_pe	reads/my_short_pe_R1.fastq.gz	reads/my_short_pe_R2.fastq.gz
```

`config/run_sheet.tsv` (Flye-compatible long-read example):

```tsv
group	sample	data_type	read1
my_nano_raw	long	nano-raw	reads/my_nano_raw.fastq.gz
```

`config/run_sheet.tsv` (hybrid example):

```tsv
group	sample	data_type	read1	read2
my_hybrid	short	short_pe	reads/my_hybrid_R1.fastq.gz	reads/my_hybrid_R2.fastq.gz
my_hybrid	long	nano-raw	reads/my_hybrid_long.fastq.gz	""
```

## Inline short-read auto-preprocess inference config

The `auto_preprocess` section of the complete workflow `config/config.yaml` configures the short-read-only Sequali-based inference script:

- script: `workflow/scripts/infer_preprocess_params.py`

For direct standalone script use, `--config` now takes the complete workflow
config and extracts its `auto_preprocess` section. The removed inference-only
configuration layout is intentionally unsupported. For example:

```bash
python workflow/scripts/infer_preprocess_params.py \
  --sequali-json path/to/report.json \
  --data-type short_pe \
  --config config/config.yaml \
  --config-schema workflow/schemas/config.schema.yaml \
  --sequali-schema workflow/schemas/sequali_short_read.schema.yaml \
  --output path/to/auto_params.json
```

The script validates its config and Sequali inputs against those schemas and
emits a canonical recommendation JSON containing exactly:

- `recommendations`
- `tool_args`
- `warnings`

`workflow/schemas/auto_preprocess_output.schema.yaml` is the tested generated-output contract; the script does not load the output schema at runtime.

Detailed inline validation for every supplied inference setting under
`auto_preprocess` is part of `workflow/schemas/config.schema.yaml` and runs
during workflow setup in `workflow/rules/common.smk` (Snakemake `validate`),
before any jobs execute. The complete workflow config must contain
`auto_preprocess`; its seven baseline scientific settings may be omitted only
when `enabled` is explicitly `false`. The optional `mott_positive_weight` and
`mott_negative_weight` settings may be omitted in enabled configurations and
then execute as `1.0`/`1.0`. The canonical enabled config explicitly shows all
nine controls. Pure-long configurations that omit the baseline settings should
therefore set `enabled: false`; long-read routing remains manual either way.

`workflow/schemas/sequali_short_read.schema.yaml` is the runtime input contract
for the Sequali fields consumed by inference. It pins `meta.sequali_version` to
`1.0.2`, validates the R1 structures and any complete R2 set, and permits
unrelated report sections. The rule declares both schema files as inputs via
`workflow.source_path()`, so local, cached-source, and imported-module execution
use the module-owned contracts. Upgrading Sequali requires reviewing and updating
this schema and the inference tests before changing the QC wrapper version.

Terminology and units:

| Concept                        | Contract                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      |
| ------------------------------ | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Composition positions and bins | Base position 1 is the 5' start. A Sequali composition label such as `31` is one bin, and a range such as `31-35` is one inclusive bin with start 31 and end 35. Internal series indices are explicitly 0-based. When composition headcrop is opted into, its unchanged scan starts at the middle available series bin by index, not at the physical half-read base position.                                                                                                                                                                                                                                                                                                                                                                                 |
| Tailcrop anchors               | Tailcrop inference uses only the mean series in `per_position_mean_quality_and_spread.end_percentiles`. Sequali 1.0.2 emits exactly 100 slots ordered from anchor `-100` through the 3'-most base at `-1`. For each mate, exactly the final `min(summary.maximum_length, 100)` slots are observed; any earlier slots are ignored.                                                                                                                                                                                                                                                                                                                                                                                                                             |
| Lengths and crops              | Bases. `headcrop_bases` is the number removed from the 5' start; `tailcrop_base_position` is the absolute original-read retained-length cutoff before the removed 3' suffix; `maximum_length_after_headcrop_bases` is relative to the post-headcrop read.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     |
| Tailcrop suffix score          | For each observed anchor, let `difference = tailcrop_mott_quality_cutoff - mean_quality` (default cutoff 20). Its contribution is `difference * mott_positive_weight` when `difference > 0`, `difference * mott_negative_weight` when `difference < 0`, and exactly `0` when `difference == 0`. Both weights are nonnegative and default/fall back to `1.0`. Starting at anchor `-1`, cumulative scores are evaluated 3'-inward without reset. The greatest positive score determines the crop; ties select the shortest suffix, and a nonpositive maximum produces no candidate. Multiplying both weights by the same positive factor only rescales every cumulative score, so their relative ratio—not common positive scaling—is operationally meaningful. |
| Tailcrop safeguard             | For mate `i`, nominal length `Li` is the corresponding raw Sequali summary `maximum_length`, and `mate_floor_i = ceil(Li * tailcrop_min_retained_fraction)`. `ceil` always rounds upward to a whole base. The default fraction is 0.75, and a non-null fastp maximum retains at least one base.                                                                                                                                                                                                                                                                                                                                                                                                                                                               |
| Base imbalance                 | `abs(A-T)`, `abs(G-C)`, their maximum, and thresholds are dimensionless fractions in [0,1]. The opt-in legacy scan continues only while both differences are strictly less than the threshold; an imbalance equal to the threshold triggers.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  |
| Quality                        | Tailcrop calculations use Sequali-reported mean quality directly as Phred-score differences. `average_quality` and its fastp argument use a tool-specific Phred-scaled threshold; read-average calculations differ among tools and are not strictly comparable.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               |

Short-read inference behavior:

- short-read recommendations:
  - `recommendations.read1` always contains the R1 recommendation;
    `recommendations.read2` contains the independent R2 recommendation for
    `short_pe` and is null for `short_se`. There are no merged scalar
    recommendation aliases.
  - each mate recommendation has `headcrop_bases`, `tailcrop_base_position`,
    and `maximum_length_after_headcrop_bases`
  - `headcrop_bases` is the number of bases removed from that mate's 5' start;
    it is zero by default because `composition_headcrop_enabled: false` skips
    the composition scan. Setting the boolean to `true` explicitly opts into
    the previous composition algorithm independently for each mate: it scans
    backward from the series-bin midpoint, uses the inclusive trigger-bin end
    and `max_headcrop`, and stops at the first triggering bin. No quality-based
    headcrop is performed.
  - for each mate, scan observed exact anchors from `-1` toward `-100` with
    cumulative score initialized to zero. At each anchor, first compute
    `difference = tailcrop_mott_quality_cutoff - mean_quality`; add
    `difference * mott_positive_weight` when positive,
    `difference * mott_negative_weight` when negative, and `0` when equal to
    zero. Accumulation does not reset. Record a span only when its cumulative
    score is strictly greater than the previous best. If the best score remains
    nonpositive, do not create a tailcrop candidate; otherwise use
    `nominal_length - best_span` as the raw cutoff
  - both Mott weights are optional nonnegative numbers and omission behaves as
    `1.0`/`1.0`, preserving the unweighted result. A zero positive weight
    removes poor-quality evidence, a zero negative weight removes good-quality
    offset, and both zero produce no positive candidate. Only the relative
    positive-to-negative ratio changes candidate selection; multiplying both by
    the same positive factor does not
  - each mate's fractional floor is
    `ceil(mate_nominal_length * tailcrop_min_retained_fraction)`, using
    `summary.maximum_length` for R1 and `summary_read2.maximum_length` for R2
  - `tailcrop_min_retained_fraction` is the only configurable tail-retention
    floor. Each non-null mate candidate is independently the maximum of its raw
    cutoff, its own fractional floor, and the one-base technical floor
  - each mate's `tailcrop_base_position` is its safeguarded candidate, or null
    when inactive
  - each mate's `maximum_length_after_headcrop_bases` is its tailcrop cutoff
    minus its own 5' headcrop, with a one-base technical floor, or null when
    tailcrop is inactive
  - each mate uses the mean series in
    `per_position_mean_quality_and_spread.end_percentiles`; `short_pe` uses the
    corresponding `_read2` section without merging mate recommendations
  - `short_pe` emits `--trim_front1/2` and `--max_len1/2` only from the
    corresponding mate's nonzero/non-null values; one mate's null tailcrop does
    not suppress the other mate's argument
  - if a mate's safeguarded cutoff is **greater than or equal to** its own
    nominal length, that mate's tailcrop is suppressed as a no-op: its
    recommendation cutoff/maximum-length fields are null and its `--max_len`
    argument is omitted
  - configured, not inferred, read-level filters emit:
    `--average_qual` (`auto_preprocess.average_quality`) and
    `--length_required` (`auto_preprocess.minimum_length`)
- no adapter inference in this phase

Long-read preprocessing is not auto-inferred. Configure Chopper manually with `chopper.extra`; generic `long` and Nanopore routes also retain their manually configured porechop_abi behavior.

Warnings use concise `code` and `message` fields, name every affected mate, and
occur at most once per category. Safeguard adjustment precedes no-op
suppression, and an explicit opt-in composition-headcrop warning remains last.
An empty warning list is valid. Successful inference prints final
recommendations, rendered tool arguments, and warnings (or `none`) to stdout.
The workflow captures this compact summary, along with any failure stderr, in
`results/{group}/logs/auto_preprocess/{group}_{sample}.auto_params.log`.

Key knobs in the flat `auto_preprocess` mapping of the complete workflow config:

- `auto_preprocess.composition_headcrop_enabled`
- `auto_preprocess.max_headcrop` (used only when composition headcrop is enabled)
- `auto_preprocess.headcrop_base_imbalance_threshold`
- `auto_preprocess.tailcrop_mott_quality_cutoff`
- `auto_preprocess.mott_positive_weight` (optional; default/fallback `1.0`)
- `auto_preprocess.mott_negative_weight` (optional; default/fallback `1.0`)
- `auto_preprocess.tailcrop_min_retained_fraction`
- `auto_preprocess.average_quality`
- `auto_preprocess.minimum_length`

Workflow integration behavior:

- Auto-inference is short-read-only. Pure long modes
  (`long`, `nano-raw`, `nano-corr`, `nano-hq`, `pacbio-raw`, `pacbio-corr`,
  `pacbio-hifi`) do not schedule inference or create auto-parameter JSON even
  when `auto_preprocess.enabled: true`.
- For every eligible short row, `auto_preprocess.enabled: true` runs inference;
  long rows have no inference jobs or auto-parameter targets.
- For each eligible short row, the workflow runs
  `workflow/scripts/infer_preprocess_params.py` after raw Sequali QC and writes
  `results/{group}/reads/{group}_{sample}.auto_params.json`; its text log is
  `results/{group}/logs/auto_preprocess/{group}_{sample}.auto_params.log`.
- Inferred fastp args are prepended and manual `fastp.extra` is appended
  (**inferred first, manual last**).
- Chopper uses only manual `chopper.extra`; no inferred long-read arguments are
  passed by the workflow.
- No inferred args are applied to `porechop_abi`.
- Adapter auto-inference remains disabled in this phase.
