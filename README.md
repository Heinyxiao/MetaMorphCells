# MetaMorphCells: Unveiling Cancer Cell Dedifferentiation
**MetaMorphCells** is a bioinformatics research project focused on studying a small population of cancer cells that undergo dedifferentiation into cancer stem cells. This repository contains various scripts and notebooks for preprocessing, analyzing, and visualizing single-cell sequencing data to explore the molecular mechanisms underlying cancer cell dedifferentiation.

## Project Overview
Cancer cells can sometimes "rewind" their development, changing into a more primitive, stem-like state. This transformation gives them special abilities, like growing faster, resisting treatments, and spreading to other parts of the body. 

In this project, we use advanced data analysis techniques to study these changes in ovarian cancer cells. By examining individual cells, we aim to uncover how they transform and find clues that might help us develop better treatments in the future. We use single-cell genomic sequencing data to model and predict cancer cell behavior. The graph below illustrates the predicted pathways of how different cancer cell clusters differentiate, based on data from a clinical dataset (GSE165897).

<img src="https://github.com/Heinyxiao/MetaMorphCells/blob/main/pseudotime_vector_field.png" alt="Pseudotime Vector Field" width="300"/>


🔬 Key Topics:
- Cancer stem cells
- Dedifferentiation
- Single-cell RNA sequencing (scRNA-seq)
- Bioinformatics tools and workflows

## Repository Structure
This repository is organized into the following sections:

  **1. Data Preprocessing**
  
  - Scripts for data cleaning and normalization of scRNA-seq datasets.
  
  **2. Analysis Pipeline**
  
  - Includes the full pipeline from clustering, differential expression, imputation, to trajectory inference.
  
  **3. Gene Regulatory Network Inference**
  
  - Tools for building gene regulatory networks and identifying key regulators.
  
  **4. Perturbation Analysis**
  
  - Scripts to model and analyze perturbations in the dedifferentiation process.
  
  **5. Visualization**
  
  - Custom visualizations for cell clusters, gene expression, and trajectories.

## Key Tools & Methods

| **Tool/Method** | **Description** |
|-----------------|------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **ALRA**        | Implements the ALRA (Adaptively-thresholded Low Rank Approximation) imputation method [ALRA GitHub](https://github.com/KlugerLab/ALRA)                                         |
| **CytoTRACE2**  | Predicts differentiation state of single cells based on transcriptional data. [CytoTRACE2 GitHub](https://github.com/digitalcytometry/cytotrace2)                             |
| **scGPT**       | Applying the scGPT model for gene expression prediction, GRN inference, and perturbation analysis. [scGPT GitHub](https://github.com/bowang-lab/scGPT)                         |
| **scPopcorn**   | Identification of unique cell clusters using the scPopcorn package. [scPopcorn GitHub](https://github.com/ncbi/scPopCorn)                                                     |
| **scTour**      | Model training and lineage trajectory inference using scTour. [scTour GitHub](https://github.com/LiQian-XC/sctour)                                                            |
| **Scanpy**      | Comprehensive scRNA-seq analysis toolkit. [Scanpy GitHub](https://github.com/scverse/scanpy)                                                                                  |
| **scVelo**      | RNA velocity analysis. [scVelo GitHub](https://github.com/theislab/scvelo_notebooks)                                                                                          |
| **Velocyto**    | RNA velocity analysis in scRNA-seq data. [Velocyto GitHub](https://github.com/velocyto-team/velocyto.py)                                                                      |

## Key Results (Paper under preperation)
- Identification of dedifferentiation markers in ovarian cancer cells.
- RNA velocity maps showing trajectory of dedifferentiation.
- Gene regulatory networks highlighting potential therapeutic targets.

## Future Directions
This project is still evolving. Future updates may include:

- Integration with spatial transcriptomics data.
- Deeper perturbation studies using additional clinical datasets.
