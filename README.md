# MetaMorphCells: Unveiling Cancer Cell Dedifferentiation
**MetaMorphCells** is a bioinformatics research project focused on studying a small population of cancer cells that undergo dedifferentiation into cancer stem cells (CSCs). This repository contains scripts and notebooks for preprocessing, analyzing, and visualizing single-cell sequencing data to explore the molecular mechanisms underlying dedifferentiation, therapy resistance, and stemness in ovarian cancer.

## Project Overview
Cancer cells can sometimes "rewind" their development, reverting to a more primitive, stem-like state. This transformation gives them survival advantages such as therapy resistance, higher proliferative capacity, and metastatic potential.

In this project, we combine single-cell RNA-seq, spatial transcriptomics, gene regulatory network inference, and perturbation modeling to map and predict dedifferentiation trajectories in ovarian cancer. By integrating classical methods and foundation models, we aim to identify molecular drivers and therapeutic vulnerabilities.

<img src="https://github.com/Heinyxiao/MetaMorphCells/blob/main/pseudotime_vector_field.png" alt="Pseudotime Vector Field" width="300"/>


## Key Topics:
- Cancer stem cells & dedifferentiation
- Single-cell transcriptomics
- Gene regulatory network (GRN) inference
- Perturbation modeling & foundation models

## Repository Structure  

**1. Data Preprocessing**  
- `GSE222557_data_preprocess.ipynb`  
- `h5ad_file_process.py`  

**2. Analysis Pipeline**  
- `scanpy_OC.ipynb`  
- `scPopcorn.ipynb`  
- `scTour_infer.py`, `scTour_model_training.py`  

**3. Gene Regulatory Network Inference**  
- `GRN_in_house_human.py`  
- `SCimilarity_Gene_Attribution.ipynb`  

**4. Perturbation & Foundation Models**  
- `geneformer_perturbation.py`  
- `Tutorial_Perturbation.ipynb`  
- `ESM-2.ipynb`  
- `Cell_Type_Classification_Fine_Tuning.ipynb`  

**5. Visualization**  
- `Attention_Visualization.ipynb`  
- `pseudotime_vector_field.png`  

## Key Tools & Methods

| Tool/Method      | Description |  
|------------------|-------------|  
| **ALRA**         | Low-rank imputation of scRNA-seq data. |  
| **CytoTRACE2**   | Differentiation potential inference. |  
| **scGPT**        | Foundation model for scRNA-seq (expression prediction, GRN inference, perturbation). |  
| **Geneformer**   | Transformer-based biological foundation model for perturbation and fine-tuning. |  
| **SCimilarity**  | Foundation model for cross-dataset similarity & gene attribution. |  
| **ESM-2**        | Protein language model for structural biology and ligand–receptor inference. |  
| **scPopcorn**    | Rare/unique cluster identification. |  
| **scTour**       | Lineage trajectory inference with deep generative models. |  
| **Scanpy**       | Comprehensive scRNA-seq analysis toolkit. |  
| **scVelo / Velocyto** | RNA velocity analysis for dynamic inference. |  

## Key Results (Paper under preperation)
- Identification of dedifferentiation markers in ovarian cancer cells.  
- RNA velocity and pseudotime maps showing dedifferentiation trajectories.  
- Gene regulatory networks highlighting candidate therapeutic targets.  
- Perturbation analysis predicting vulnerabilities of WNT5A-CAF crosstalk.  
- Early integration with **spatial transcriptomics** confirming regional CSC enrichment.  


## Future Directions
This project is still evolving. Future updates may include:
- Large-scale perturbation analysis with foundation models (Geneformer, scGPT).  
- Integration of **multi-omics (ATAC, CUT&RUN, proteomics)** for regulatory inference.  
- Structural modeling of ligand–receptor pairs with **ESM-2 + AlphaFold2**.  
- Clinical dataset integration for biomarker discovery and patient stratification.  
