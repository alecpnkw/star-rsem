
# star-rsem

A [Snakemake](https://snakemake.readthedocs.io/) workflow for RNA-seq analysis using [STAR](https://github.com/alexdobin/STAR) for alignment and [RSEM](https://github.com/deweylab/RSEM) for quantification.

## Overview

This workflow automates the following steps:
- Quality control with [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) and [MultiQC](https://multiqc.info/)
- Genome index generation with [STAR](https://github.com/alexdobin/STAR)
- Paired-end alignment with [STAR](https://github.com/alexdobin/STAR)
- Gene and transcript quantification with [RSEM](https://github.com/deweylab/RSEM)

## Usage

1. Clone this repository.
2. Prepare your configuration files:
   - `config/config.yaml`: Main workflow configuration.
   - `config/samples.tsv`: Tab-separated sample sheet with columns for sample names and FASTQ file paths.
3. Edit the config file to specify genome reference files and sample sheet location.
4. Preview the workflow:
   ```bash
   snakemake --dry-run
   ```
5. Execute the workflow 

## Configuration

A default config file path is located at `config/config.yaml`. The sample sheet must have a `sample` column and columns for R1 and R2 FASTQ files.

## Output

- STAR genome index: `results/star-index/`
- STAR alignments: `results/star-pe/`
- RSEM quantification: `results/rsem-combined/`
- QC reports: `results/fastqc/`, `results/multiqc/`

## Notes

Environment management is supported via conda or environment modules (see comments in rule files). This workflow is still a work in progress.