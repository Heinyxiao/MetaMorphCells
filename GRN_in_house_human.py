
import copy
import json
import os
from pathlib import Path
import sys
import warnings

import torch
from anndata import AnnData
import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import networkx as nx
import pandas as pd
import tqdm
import gseapy as gp

from torchtext.vocab import Vocab
from torchtext._torchtext import (
    Vocab as VocabPybind,
)

sys.path.insert(0, "../")
import scgpt as scg
from scgpt.tasks import GeneEmbedding
from scgpt.tokenizer.gene_tokenizer import GeneVocab
from scgpt.model import TransformerModel
from scgpt.preprocess import Preprocessor
from scgpt.utils import set_seed 

os.environ["KMP_WARNINGS"] = "off"
warnings.filterwarnings('ignore')

set_seed(42)
pad_token = "<pad>"
special_tokens = [pad_token, "<cls>", "<eoc>"]
n_hvg = 2000
n_bins = 15
mask_value = -1
pad_value = -2
n_input_bins = n_bins

# Specify data path
data_dir = Path("/N/slate/xuexiao/GSE165897")
adata = sc.read(
    str(data_dir / "Patient_EOC.h5ad"))  

ori_batch_col = "patient_id"

#adata.obs["celltype"] = adata.obs["cell_subtype"].astype(str)
data_is_raw = True

# Specify output path
output_path = "/N/slate/xuexiao/GSE165897/scGPT_result/"
gene_program_heatmap_name = "gene_program_activation_scores_human_GSE165897.png"
cos_similarity_file_name = "cosine_similarity_network_human_GSE165897.png"

def get_full_file_path(output_path, file_name):
    return f"{output_path}{file_name}"

# Step 1: Load pre-trained model and dataset
## 1.1 Load pre-trained model
# Specify model path; here we load the pre-trained scGPT blood model
model_dir = Path("/N/u/xuexiao/Quartz/scripts/AI_tools/scGPT_human")
model_config_file = model_dir / "args.json"
model_file = model_dir / "best_model.pt"
vocab_file = model_dir / "vocab.json"

vocab = GeneVocab.from_file(vocab_file)
for s in special_tokens:
    if s not in vocab:
        vocab.append_token(s)

# Retrieve model parameters from config files
with open(model_config_file, "r") as f:
    model_configs = json.load(f)
print(
    f"Resume model from {model_file}, the model args will override the "
    f"config {model_config_file}."
)
embsize = model_configs["embsize"]
nhead = model_configs["nheads"]
d_hid = model_configs["d_hid"]
nlayers = model_configs["nlayers"]
n_layers_cls = model_configs["n_layers_cls"]

gene2idx = vocab.get_stoi()

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

ntokens = len(vocab)  # size of vocabulary
model = TransformerModel(
    ntokens,
    embsize,
    nhead,
    d_hid,
    nlayers,
    vocab=vocab,
    pad_value=pad_value,
    n_input_bins=n_input_bins,
)

try:
    model.load_state_dict(torch.load(model_file))
    print(f"Loading all model params from {model_file}")
except:
    # only load params that are in the model and match the size
    model_dict = model.state_dict()
    pretrained_dict = torch.load(model_file)
    pretrained_dict = {
        k: v
        for k, v in pretrained_dict.items()
        if k in model_dict and v.shape == model_dict[k].shape
    }
    for k, v in pretrained_dict.items():
        print(f"Loading params {k} with shape {v.shape}")
        model_dict.update(pretrained_dict)
        model.load_state_dict(model_dict)

model.to(device)

## 1.2 Load dataset of interest

# Preprocess the data following the scGPT data pre-processing pipeline
preprocessor = Preprocessor(
    use_key="X",  # the key in adata.layers to use as raw data
    filter_gene_by_counts=3,  # step 1
    filter_cell_by_counts=False,  # step 2
    normalize_total=1e4,  # 3. whether to normalize the raw data and to what sum
    result_normed_key="X_normed",  # the key in adata.layers to store the normalized data
    log1p=data_is_raw,  # 4. whether to log1p the normalized data
    result_log1p_key="X_log1p",
    subset_hvg=n_hvg,  # 5. whether to subset the raw data to highly variable genes
    hvg_flavor="seurat_v3" if data_is_raw else "cell_ranger",
    binning=n_bins,  # 6. whether to bin the raw data and to what number of bins
    result_binned_key="X_binned",  # the key in adata.layers to store the binned data
)
preprocessor(adata, batch_key = ori_batch_col)

# Step 2: Retrieve scGPT's gene embeddings

# Retrieve the data-independent gene embeddings from scGPT
gene_ids = np.array([id for id in gene2idx.values()])
gene_embeddings = model.encoder(torch.tensor(gene_ids, dtype=torch.long).to(device))
gene_embeddings = gene_embeddings.detach().cpu().numpy()

# Filter on the intersection between the Immune Human HVGs found in step 1.2 and scGPT's 30+K foundation model vocab
gene_embeddings = {gene: gene_embeddings[i] for i, gene in enumerate(gene2idx.keys()) if gene in adata.var.index.tolist()}
print('Retrieved gene embeddings for {} genes.'.format(len(gene_embeddings)))

# Construct gene embedding network
embed = GeneEmbedding(gene_embeddings)

# Step 3: Extract gene programs from gene embedding network
##3.1 Perform Louvain clustering on the gene embedding network
# Perform Louvain clustering with desired resolution; here we specify resolution=40
gdata = embed.get_adata(resolution=40)
# Retrieve the gene clusters
metagenes = embed.get_metagenes(gdata)

##3.2 Filter on clusters with 5 or more genes
# Obtain the set of gene programs from clusters with #genes >= 5
mgs = dict()
for mg, genes in metagenes.items():
    if len(genes) > 4:
        mgs[mg] = genes
# Here are the gene programs identified
print(mgs)

# Step 4: Visualize gene program activation on the Immune Human dataset
sns.set(font_scale=0.35)
embed.score_metagenes(adata, metagenes)
embed.plot_metagenes_scores(adata, mgs, "cell_subtype")

fig = plt.gcf()
ax = plt.gca()

# Customize yticklabels for more margin between rows
plt.setp(ax.get_yticklabels(), rotation=0, ha="right", rotation_mode="anchor", fontsize=10)

# Add more margin by setting labelpad
ax.tick_params(axis='y', pad=10)  # Adjust pad value as needed


save_file_path = get_full_file_path(output_path, gene_program_heatmap_name)
print(f"Saving plot to: {save_file_path}")
plt.savefig(save_file_path, dpi=300, bbox_inches=None)
plt.close()

# Step 5: Visualize network connectivity within desired gene program
# Retrieve gene program 3 which contains the CD3 gene set
CD_genes = mgs['6']
print(CD_genes)

# Initialize an empty list to collect DataFrames
df_list = []

# Compute cosine similarities among genes in this gene program
#df_CD = pd.DataFrame(columns=['Gene', 'Similarity', 'Gene1'])

for i in tqdm.tqdm(CD_genes):
    df = embed.compute_similarities(i, CD_genes)
    df['Gene1'] = i
    df_list.append(df)
df_CD = pd.concat(df_list, ignore_index=True)
df_CD_sub = df_CD[df_CD['Similarity']<0.99].sort_values(by='Gene') # Filter out edges from each gene to itself
print("df_CD_sub is", df_CD_sub)

# Creates a graph from the cosine similarity network
input_node_weights = [(row['Gene'], row['Gene1'], round(row['Similarity'], 2)) for i, row in df_CD_sub.iterrows()]
G = nx.Graph()
G.add_weighted_edges_from(input_node_weights)
# Plot the cosine similarity network; strong edges (> select threshold) are highlighted
thresh = 0.3
plt.figure(figsize=(20, 20))
widths = nx.get_edge_attributes(G, 'weight')

elarge = [(u, v) for (u, v, d) in G.edges(data=True) if d["weight"] > thresh]
esmall = [(u, v) for (u, v, d) in G.edges(data=True) if d["weight"] <= thresh]

pos = nx.spring_layout(G, k=0.4, iterations=15, seed=3)

width_large = {}
width_small = {}
for i, v in enumerate(list(widths.values())):
    if v > thresh:
        width_large[list(widths.keys())[i]] = v*10
    else:
        width_small[list(widths.keys())[i]] = max(v, 0)*10

nx.draw_networkx_edges(G, pos,
                       edgelist = width_small.keys(),
                       width=list(width_small.values()),
                       edge_color='lightblue',
                       alpha=0.8)
nx.draw_networkx_edges(G, pos, 
                       edgelist = width_large.keys(), 
                       width = list(width_large.values()), 
                       alpha = 0.5, 
                       edge_color = "blue", 
                      )
# node labels
nx.draw_networkx_labels(G, pos, font_size=25, font_family="sans-serif")
# edge weight labels
d = nx.get_edge_attributes(G, "weight")
edge_labels = {k: d[k] for k in elarge}
nx.draw_networkx_edge_labels(G, pos, edge_labels, font_size=15)

ax = plt.gca()
ax.margins(0.08)
plt.axis("off")

# Save the plot to a file
cos_simolarity_output_path = get_full_file_path(output_path, cos_similarity_file_name)
print(f"Saving second plot to: {cos_simolarity_output_path}")
plt.savefig(cos_simolarity_output_path, dpi=300, bbox_inches=None)
plt.close()


# Step 6: Reactome pathway analysis
# Meta info about the number of terms (tests) in the databases
df_database = pd.DataFrame(
data = [['GO_Biological_Process_2021', 6036],
['GO_Molecular_Function_2021', 1274],
['Reactome_2022', 1818]],
columns = ['dataset', 'term'])
# Select desired database for query; here use Reactome as an example
databases = ['GO_Biological_Process_2021']
m = df_database[df_database['dataset'].isin(databases)]['term'].sum()
# p-value correction for total number of tests done
p_thresh = 0.05/m
# Perform pathway enrichment analysis using the gseapy package in the Reactome database
df_list = []

# Perform pathway enrichment analysis using the gseapy package in the Reactome database
enr_Reactome = gp.enrichr(gene_list=CD_genes,
                          gene_sets=databases,
                          organism='Human', 
                          outdir='test/GO_Biological_Process_2021',
                          cutoff=0.5)

# Get the results and filter based on the p-value threshold
out = enr_Reactome.results
out_filtered = out[out['P-value'] < p_thresh]

# Append the filtered results to the list
df_list.append(out_filtered)

# Concatenate all DataFrames in the list
df = pd.concat(df_list, ignore_index=True)

# Print the resulting DataFrame to verify
print(df)
