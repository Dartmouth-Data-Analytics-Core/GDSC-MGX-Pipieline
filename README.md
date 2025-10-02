# GDSC-MGX-Pipieline
Snakemake pipeline for preprocessing of metagenomics using Biobakery tools

## Introduction 
The pipeline is designed to provide efficient pre-processing of metagenomics samples (MGX) data on high performance computing clusters (HPCs) using the Torque/PBS scheduler or using a single high CPU, high RAM machine, for the *Genomic Data Science Core (GDSC)* of the *Center for Quantitative Biology (CQB)*, located at Dartmouth College. Both single- and paired-end datasets are supported. The pipeline has been built and tested using human data sets. Required software can be installed using Conda with the enrionment file (environment.yml), or specified as paths in the config.yaml file.

<img src="img/cqb_logo.jpg" width="250" height="140" >

## Pipeline summary:
This pipeline makes use of [The Biobakery](https://github.com/biobakery) The major steps implmented in the pipeline include: 

- FASTQ quality control assesment using [*FASTQC*](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- Host read filtering using [*kneaddata*](https://huttenhower.sph.harvard.edu/kneaddata/)
- Taxonomic classification using [*metaphlan*](https://huttenhower.sph.harvard.edu/metaphlan/)
- Function classification using [*humann*](https://huttenhower.sph.harvard.edu/humann/) 
