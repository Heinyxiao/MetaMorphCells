import anndata as ad
import pandas as pd
import os

directory_path = '/Users/xuexiao/Library/CloudStorage/GoogleDrive-heinyxiao@gmail.com/My Drive/Lab/Projects/Dedifferentiation/Data'

# Load the .h5ad file
adata = ad.read_h5ad('/Users/xuexiao/Library/CloudStorage/GoogleDrive-heinyxiao@gmail.com/My Drive/Lab/Projects/Dedifferentiation/Data/Patient_EOC_with_ptime.h5ad')
#print(adata)

# Save the main expression matrix and metadata
#adata.to_df().to_csv(os.path.join(directory_path, 'GSE165897_EOC_expression_matrix.csv'))
#adata.obs.to_csv(os.path.join(directory_path, 'GSE165897_EOC_metadata.csv'))



# Load the new metadata
new_metadata = pd.read_csv(os.path.join(directory_path, "GSE165897_EOC_new_metadata.csv"), index_col=0)

# Update the metadata
for col in new_metadata.columns:
    adata.obs[col] = new_metadata[col]
print(adata)
adata.write_h5ad(os.path.join(directory_path,"GSE165897_EOC_with_ptime_cytotrace2.h5ad"))
quit()