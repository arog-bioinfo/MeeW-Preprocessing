# Snakemake pipeline: MeeW-Preprocessing

[![Snakemake](https://img.shields.io/badge/Snakemake-%E2%89%A59.12.0-039475?logo=snakemake&logoColor=white)](https://snakemake.readthedocs.io/)
[![run with conda](https://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/projects/conda/en/latest/)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

<!-- TEMPORARY CATALOG OPT-OUT: later replace every standalone "pipeline" with the standard term formed by "work" + "flow", then remove this comment. -->

A Snakemake pipeline for quality control and preprocessing of short and long sequencing reads.

> [!WARNING]
> **Research preview:** This pipeline is under active development and is not a deployment-ready release. Validate its results and resource requirements for your use case.

## Usage

The pipeline accepts any positive number of run-sheet rows across any number of groups. Supported modes are single-end or paired-end short reads and generic, Nanopore, or PacBio long reads. Every route runs Sequali before and after preprocessing. Short reads use optional Sequali-based parameter inference followed by fastp; generic and Nanopore reads use porechop_abi and Chopper, while PacBio reads go directly to Chopper.

Install Snakemake `>=9.12.0` and Conda, then edit [`config/config.yaml`](config/config.yaml) and its referenced run sheet to use your reads. Relative input paths are resolved from the Snakemake working directory. Final reads and QC reports are grouped below `results/{group}/`.

See the [configuration guide](config/README.md) for all run-sheet modes, inference settings, algorithms, and output contracts. Run the commands below from the repository root.

## Deployment options

Conda is the supported software deployment method. After configuring `config/config.yaml`, inspect the planned jobs and then execute them with:

```bash
snakemake --cores 2 --dry-run
snakemake \
  --cores 2 \
  --software-deployment-method conda
```

Adjust cores and resource overrides for the data and execution environment. This pipeline does not supply a container deployment.

## Execution profiles

No pipeline-specific profile is supplied. A user-maintained [Snakemake profile](https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles) may provide executor and site-specific resource settings.

## License

The pipeline is distributed under the [MIT License](LICENSE).

## References

No module-specific DOI is declared. Cite the repository version used, the tools used by the selected route—[Sequali](https://github.com/rhpvorderman/sequali), [fastp](https://github.com/OpenGene/fastp), [porechop_abi](https://github.com/bonsai-team/Porechop_ABI), and [Chopper](https://github.com/wdecoster/chopper)—and Snakemake:

> Köster, J., Mölder, F., Jablonski, K. P., Letcher, B., Hall, M. B., Tomkins-Tinch, C. H., Sochat, V., Forster, J., Lee, S., Twardziok, S. O., Kanitz, A., Wilm, A., Holtgrewe, M., Rahmann, S., & Nahnsen, S. _Sustainable data analysis with Snakemake_. F1000Research, 10:33, **2021**. https://doi.org/10.12688/f1000research.29032.2.
