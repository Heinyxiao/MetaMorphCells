
import copy
import json
import os
from pathlib import Path
import sys
import warnings

#import torch
from anndata import AnnData
import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import networkx as nx
import pandas as pd
import tqdm
import gseapy as gp

#from torchtext.vocab import Vocab
#from torchtext._torchtext import (
#    Vocab as VocabPybind,
#)

sys.path.insert(0, "../")
import scgpt as scg
from scgpt.tasks import GeneEmbedding
from scgpt.tokenizer.gene_tokenizer import GeneVocab
from scgpt.model import TransformerModel
from scgpt.preprocess import Preprocessor
from scgpt.utils import set_seed 

os.environ["KMP_WARNINGS"] = "off"
warnings.filterwarnings('ignore')

def get_full_file_path(output_path, file_name):
    return f"{output_path}{file_name}"


set_seed(42)
pad_token = "<pad>"
special_tokens = [pad_token, "<cls>", "<eoc>"]
n_hvg = 50
n_bins = 20
mask_value = -1
pad_value = -2
n_input_bins = n_bins

## 1.2 Load dataset of interest
# Specify data path
adata = sc.read('/N/slate/xuexiao/combine_all_3/OC_cell_in_house.h5ad')
adata.shape
print(adata.obs.info())

sc.pp.filter_cells(adata, min_genes=100)
sc.pp.filter_genes(adata, min_cells=3)

# Normalize the data
sc.pp.normalize_total(adata, target_sum=1e4)

# Logarithmize the data
sc.pp.log1p(adata)

sc.pp.highly_variable_genes(adata, n_top_genes=2000, batch_key="orig.ident")

adata

# PCA
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)

# UMAP
sc.tl.umap(adata)

output_path = "/N/slate/xuexiao/combine_all_3/scTour_result/"

sc.pl.umap(adata, color=['cell_subtype', 'orig.ident'], legend_loc='on data')

batch_effect_file_name = "in_house_batch_effect.png"
batch_effect_output_path = get_full_file_path(output_path, batch_effect_file_name)
print(f"Saving batch effect plot to: {batch_effect_output_path}")
plt.savefig(batch_effect_output_path, dpi=300, bbox_inches=None)
plt.close()

adata.write('/N/slate/xuexiao/combine_all_3/preprocessed_OC_cell_in_house.h5ad')



