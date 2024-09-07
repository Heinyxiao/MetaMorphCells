# MetaMorphCells
MetaMorphCells is a bioinformatics research project focused on studying a small population of cancer cells that undergo dedifferentiation into cancer stem cells. This repository contains various scripts and notebooks for preprocessing, analyzing, and visualizing single-cell sequencing data to explore the molecular mechanisms underlying cancer cell dedifferentiation.

# Introduction
Cancer cells can undergo dedifferentiation, acquiring stem-like properties that enhance their ability to proliferate, resist therapies, and metastasize. This project uses bioinformatics tools to investigate these processes in ovarian cancer cell populations at a single-cell resolution. The repository includes scripts for data preprocessing, perturbation analysis, gene regulatory network inference, and more.

# Usage
Each script and notebook in this repository serves a specific purpose in the analysis pipeline. Below is a brief description of each file to help you understand the workflow.

# Methods

**ALRA**: Implements the ALRA (Adaptively-thresholded Low Rank Approximation) imputation method (https://github.com/KlugerLab/ALRA) for single-cell RNA-seq data.

**CytoTRACE2**: An R markdown file for running CytoTRACE2 (https://github.com/digitalcytometry/cytotrace2), a tool to predict the differentiation state of single cells based on transcriptional data. 

**scGPT**: Applying the scGPT model (https://github.com/bowang-lab/scGPT) to single-cell data for gene expression prediction and analysis, including gene regulatory network (GRN) inference and perturbation analysis tailored for ovarian cancer clinical single-cell RNA-seq datasets.

**scPopcorn**: Identification of unique cell clusters using the scPopcorn package (https://github.com/ncbi/scPopCorn).

**scTour**: Model training and lineage trajectory inference using the scTour framework (https://github.com/LiQian-XC/sctour).

**scanpy**: Scanpy (https://github.com/scverse/scanpy) is a toolkit for comprehensive single-cell RNA-seq analysis and spatial transcriptomic analysis. 

**scvelo**: Using scVelo (https://github.com/theislab/scvelo_notebooks) for RNA velocity analysis.

**velocyto**: Velocyto (https://github.com/velocyto-team/velocyto.py) is a tool for RNA velocity analysis in single-cell RNA-seq data.


