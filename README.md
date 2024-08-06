# Snakemake workflow: `hg_genotype`

[![Snakemake](https://img.shields.io/badge/snakemake-≥8-brightgreen.svg)](https://snakemake.github.io)
[![Tests](https://github.com/3d-omics/hg_genotype/actions/workflows/main.yml/badge.svg)](https://github.com/3d-omics/hg_genotype/actions/workflows/main.yml)

A Snakemake workflow for  Short Variant Discovery in Host Genomes


## Usage

- Test that it works:
  - Make sure you have installed `snakemake>=8`
  - Run the pipeline: `snakemake --use-conda --profile profile/default --jobs 100`. It will download all the necesary software through conda. It should take less than 5 minutes.

- Run it with your own data:
  - Edit `config/samples.tsv` and add your samples and where are they located.
  - Edit `config/features.tsv` with information regarding the reference you are using.
  - Run the pipeline: `snakemake --use-conda --profile profile/default --jobs 8 all`.
  - If you are in a cluster with slurm, add `--executor slurm`.

## Features

- FASTQ processing with `fastp`.
- Mapping with `bwa-mem2`
- SAM/BAM/CRAM processing with `samtools` and `GATK`.
- SNP calling with `GATK4`.
- SNP annotation with `SNPEff` and `VEP`
- Sample swap detection with `somalier`.
- Reporting with `MultiQC`.

## DAG

![host_genomics_pipeline](./schema.svg?raw=true)

## References

- [`fastp`](https://github.com/OpenGene/fastp)
- [`bwa-mem2`](https://github.com/bwa-mem2/bwa-mem2)
- [`samtools`](https://github.com/samtools/samtools)
- [`GATK`](https://github.com/broadinstitute/gatk)
- [`SNPEff`](https://github.com/pcingola/SnpEff)
- [`VEP`](https://github.com/Ensembl/ensembl-vep)
- [`somalier`](https://github.com/brentp/somalier)
- [`MultiQC`](https://github.com/MultiQC/MultiQC)
