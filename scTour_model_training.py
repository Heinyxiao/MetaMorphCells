## Running using scGPT environment
## "/N/slate/xuexiao/scgpt_env"

import sctour as sct
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

n_hvg = 3000
n_bins = 15

## Load data
adata = sc.read('/N/slate/xuexiao/GSE165897/preprocessed_Patient_EOC.h5ad')
#adata.obs['celltype'] = 'cell_subtype'
print(adata.shape)
adata.X = adata.X.astype(np.float32)

## Model training
# QC matrics
sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)
# Select highly variable genes 
sc.pp.highly_variable_genes(adata, n_top_genes=3000, batch_key="patient_id")
# Model training
tnode = sct.train.Trainer(adata, loss_mode='mse', alpha_recon_lec=0.5, alpha_recon_lode=0.5)
tnode.train()
print("Finish model training")

##The first parameter is the directory where you want to save the model, and the second parameter is the prefix for the model name.
tnode.save_model(save_dir='/N/slate/xuexiao/GSE165897/', save_prefix='scTour_model_GSE165897_EOC_3khvg')

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
adata.write("/N/slate/xuexiao/GSE165897/Patient_EOC_with_ptime.h5ad")

## Visualization
adata = adata[np.argsort(adata.obs['ptime'].values), :]
sc.pp.neighbors(adata, use_rep='X_TNODE', n_neighbors=20)
sc.tl.umap(adata, min_dist=0.1)

# Plot UMAP colored by cell type
fig, ax = plt.subplots()
sc.pl.umap(adata, color='cell_subtype', ax=ax, legend_loc='on data', show=False, frameon=False)
ax.set_title("UMAP by Cell Type")
fig.savefig("scTour_model_GSE165897_EOC_umap_celltype.png")
plt.close(fig)

# Plot UMAP colored by sample batch
fig, ax = plt.subplots()
sc.pl.umap(adata, color='patient_id', ax=ax, show=False, frameon=False)
ax.set_title("UMAP by Sample Batch")
fig.savefig("scTour_model_GSE165897_EOC_umap_sample_batch.png")
plt.close(fig)

# Plot UMAP colored by chemo status
fig, ax = plt.subplots()
sc.pl.umap(adata, color='treatment_phase', ax=ax, show=False, frameon=False)
ax.set_title("UMAP by treatment_phase")
fig.savefig("scTour_model_GSE165897_EOC_umap_treatment_phase.png")
plt.close(fig)

# Plot UMAP colored by pseudotime
fig, ax = plt.subplots()
sc.pl.umap(adata, color='ptime', ax=ax, show=False, frameon=False)
ax.set_title("UMAP by Pseudotime")
fig.savefig("scTour_model_GSE165897_EOC_umap_pseudotime.png")
plt.close(fig)

# Plot vector field
fig, ax = plt.subplots()
sct.vf.plot_vector_field(adata, zs_key='X_TNODE', vf_key='X_VF', use_rep_neigh='X_TNODE', color='cell_subtype', show=False, ax=ax, legend_loc='none', frameon=False, size=100, alpha=0.2)
ax.set_title("Vector Field")
fig.savefig("scTour_model_GSE165897_EOC_vector_field.png")
plt.close(fig)

# Create combined plot
fig, axs = plt.subplots(ncols=3, nrows=2, figsize=(15, 10))
sc.pl.umap(adata, color='cell_subtype', ax=axs[0, 0], legend_loc='on data', show=False, frameon=False)
axs[0, 0].set_title("UMAP by Cell Type")

sc.pl.umap(adata, color='patient_id', ax=axs[0, 1], show=False, frameon=False)
axs[0, 1].set_title("UMAP by Sample Batch")

sc.pl.umap(adata, color='treatment_phase', ax=axs[0, 2], show=False, frameon=False)
axs[0, 2].set_title("UMAP by treatment_phase")

sc.pl.umap(adata, color='ptime', ax=axs[1, 0], show=False, frameon=False)
axs[1, 0].set_title("UMAP by Pseudotime")

sct.vf.plot_vector_field(adata, zs_key='X_TNODE', vf_key='X_VF', use_rep_neigh='X_TNODE', color='cell_subtype', show=False, ax=axs[1, 1], legend_loc='none', frameon=False, size=100, alpha=0.2)
axs[1, 1].set_title("Vector Field")

fig.delaxes(axs[1, 2])  # Remove the empty subplot
plt.tight_layout()
fig.savefig("scTour_model_GSE165897_EOC_combined_plots.png")
plt.show()


## Use model to predict other dataset

# Load model
#tnode = sct.predict.load_model('/N/slate/xuexiao/GSE165897/scTour_model_GSE165897_EOC.pth')

# Load new dataset
adata_test1 = sc.read('/N/slate/xuexiao/combine_all_3/OC_cell_in_house.h5ad')
adata_test1.shape
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
sc.pl.umap(adata_test1, color='patient_id', show=False, ax=axs[0, 1], frameon=False)
axs[0, 1].set_title("UMAP by Sample Batch")
fig.savefig("scTour_model_GSE165897_EOC_umap_sample_batch_test1_in_house.png")

# Plot UMAP colored by chemo status
sc.pl.umap(adata_test1, color='treatment_phase', ax=axs[0, 2], show=False, frameon=False)
axs[0, 2].set_title("UMAP by Treatment Phase")
fig.savefig("scTour_model_GSE165897_EOC_umap_treatment_phase_test1_in_house.png")

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
