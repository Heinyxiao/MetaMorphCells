## Running using scGPT environment
## "/N/slate/xuexiao/scgpt_env"

import sctour as sct
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

n_hvg = 2000
n_bins = 15

## Load data
adata = sc.read('/N/slate/xuexiao/GSE165897/Patient_EOC_with_ptime.h5ad')
#adata.obs['celltype'] = 'cell_subtype'
print(adata.shape)
adata.X = adata.X.astype(np.float32)
# Display the first few rows of the metadata
print(adata.obs.head())

# Display the columns in the metadata
print(adata.obs.columns)

# Summary of the metadata
print(adata.obs.info())
# Load model
tnode = sct.predict.load_model('/N/slate/xuexiao/GSE165897/scTour_model_GSE165897_EOC.pth')


## Infer cellular dynamics
# Pseudotime
adata.obs['ptime'] = tnode.get_time()

# Latent space
#zs represents the latent z from variational inference, and pred_zs represents the latent z from ODE solver
#mix_zs represents the weighted combination of the two, which is used for downstream analysis
mix_zs, zs, pred_zs = tnode.get_latentsp(alpha_z=0.5, alpha_predz=0.5)
adata.obsm['X_TNODE'] = mix_zs

# Vector field
adata.obsm['X_VF'] = tnode.get_vector_field(adata.obs['ptime'].values, adata.obsm['X_TNODE'])

# Save data
#adata.write("/N/slate/xuexiao/GSE165897/Patient_EOC_with_ptime.h5ad")

## Visualization
adata = adata[np.argsort(adata.obs['ptime'].values), :]
sc.pp.neighbors(adata, use_rep='X_TNODE', n_neighbors=15)
sc.tl.umap(adata, min_dist=0.1)

# Create subplots
fig, axs = plt.subplots(ncols=3, nrows=2, figsize=(20, 15))

# Plot UMAP colored by cell type
sc.pl.umap(adata, color='cell_subtype', ax=axs[0, 0], legend_loc='on data', show=False, frameon=False)
axs[0, 0].set_title("UMAP by Cell Type")
fig.savefig("scTour_model_GSE165897_EOC_umap_celltype.png")

# Plot UMAP colored by sample batch
sc.pl.umap(adata, color='patient_id', ax=axs[0, 1], show=False, frameon=False)
axs[0, 1].set_title("UMAP by Sample Batch")
fig.savefig("scTour_model_GSE165897_EOC_umap_sample_batch.png")

# Plot UMAP colored by chemo status
sc.pl.umap(adata, color='treatment_phase', ax=axs[0, 2], show=False, frameon=False)
axs[0, 2].set_title("UMAP by Treatment Phase")
fig.savefig("scTour_model_GSE165897_EOC_umap_treatment_phase.png")

# Plot UMAP colored by pseudotime
sc.pl.umap(adata, color='ptime', ax=axs[1, 0], show=False, frameon=False)
axs[1, 0].set_title("UMAP by Pseudotime")
fig.savefig("scTour_model_GSE165897_EOC_umap_pseudotime.png")

# Plot vector field
sct.vf.plot_vector_field(adata, zs_key='X_TNODE', vf_key='X_VF', use_rep_neigh='X_TNODE', color='cell_subtype', show=False, ax=axs[1, 1], legend_loc='none', frameon=False, size=100, alpha=0.2)
axs[1, 1].set_title("Vector Field")
fig.savefig("scTour_model_GSE165897_EOC_vector_field.png")

# Save the entire figure
plt.tight_layout()
plt.savefig("scTour_model_GSE165897_EOC_combined_plots.png")
plt.show()

## Use model to predict other dataset


# Load new dataset
adata_test1 = sc.read('/N/slate/xuexiao/combine_all_3/OC_cell_in_house.h5ad')
print(adata_test1.shape)
adata.X = adata_test1.X.astype(np.float32)
# Display the first few rows of the metadata
print(adata_test1.obs.head())

# Display the columns in the metadata
print(adata_test1.obs.columns)

# Summary of the metadata
print(adata_test1.obs.info())


#test dataset 1
#the first parameter is the trained model
adata_test1.obs['ptime'] = sct.predict.predict_time(tnode, adata_test1)
mix_zs, zs, pred_zs = sct.predict.predict_latentsp(tnode, adata_test1, alpha_z=0.4, alpha_predz=0.6, mode='coarse')
adata_test1.obsm['X_TNODE'] = mix_zs
adata_test1.obsm['X_VF'] = sct.predict.predict_vector_field(tnode, adata_test1.obs['ptime'].values, adata_test1.obsm['X_TNODE'])

# View result
fig, axs = plt.subplots(ncols=3, nrows=2, figsize=(15, 10))

# Plot UMAP colored by cell type
sc.pl.umap(adata_test1, color='cell_subtype', legend_loc='on data', show=False, ax=axs[0, 0], frameon=False)
axs[0, 0].set_title("UMAP by Cell Type")
fig.savefig("scTour_model_GSE165897_EOC_umap_celltype_test1.png")

# Plot UMAP colored by sample batch
sc.pl.umap(adata_test1, color='orig.ident', show=False, ax=axs[0, 1], frameon=False)
axs[0, 1].set_title("UMAP by Sample Batch")
fig.savefig("scTour_model_GSE165897_EOC_umap_sample_batch_test1_in_house.png")

# Plot UMAP colored by chemo status
sc.pl.umap(adata_test1, color='CytoTRACE2_Potency', ax=axs[0, 2], show=False, frameon=False)
axs[0, 2].set_title("UMAP by Treatment Phase")
fig.savefig("scTour_model_GSE165897_EOC_umap_CytoTRACE2_Potency_test1_in_house.png")

# Plot UMAP colored by pseudotime
sc.pl.umap(adata_test1, color='ptime', show=False, ax=axs[1, 0], frameon=False)
axs[1, 0].set_title("UMAP by Pseudotime")
fig.savefig("scTour_model_GSE165897_EOC_umap_pseudotime_test1_in_house.png")

# Plot vector field with cell type color
sct.vf.plot_vector_field(adata_test1, zs_key='X_TNODE', vf_key='X_VF', use_rep_neigh='X_TNODE', show=False, ax=axs[1, 1], color='cell_subtype', t_key='ptime', frameon=False, size=80, alpha=0.05)
axs[1, 1].set_title("Vector Field by Cell Type")
fig.savefig("scTour_model_GSE165897_EOC_vector_field_celltype_test1_in_house.png")

# Plot vector field with treatment phase color
sct.vf.plot_vector_field(adata_test1, zs_key='X_TNODE', vf_key='X_VF', use_rep_neigh='X_TNODE', show=False, ax=axs[1, 2], color='treatment_phase', t_key='ptime', frameon=False, size=80, alpha=0.05)
axs[1, 2].set_title("Vector Field by Treatment Phase")
fig.savefig("scTour_model_GSE165897_EOC_vector_field_treatment_phase_test1_in_house.png")

# Save the entire figure
plt.tight_layout()
plt.savefig("scTour_model_GSE165897_EOC_combined_plots_test1_in_house.png")
plt.show()
