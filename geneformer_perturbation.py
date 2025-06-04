# Geneformer In Silico Perturbation Analysis
# Using 30M model with DCM example data

import os
import pickle
import numpy as np
import pandas as pd
from datasets import load_dataset
import torch
from transformers import BertForMaskedLM, BertForSequenceClassification
from geneformer import InSilicoPerturber, EmbExtractor
import matplotlib.pyplot as plt
import seaborn as sns

# Set device
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
print(f"Using device: {device}")

# 1. Load the DCM example dataset from HuggingFace
print("Loading DCM dataset from HuggingFace...")
dataset = load_dataset("ctheodoris/Genecorpus-30M", 
                      data_files="example_input_files/cell_classification/disease_classification/human_dcm_hcm_nf.dataset")

# Extract the dataset
train_dataset = dataset["train"]
print(f"Dataset loaded with {len(train_dataset)} samples")

# 2. Load the Geneformer 30M model
print("Loading Geneformer 30M model...")
model_name = "ctheodoris/Geneformer"  # This is the 30M model
model = BertForMaskedLM.from_pretrained(model_name, 
                                       output_attentions=False, 
                                       output_hidden_states=False)
model.to(device)
model.eval()

# 3. Extract all genes from the dataset for comprehensive perturbation
print("Extracting all genes from dataset for comprehensive perturbation screen...")

# Get all unique genes from the dataset
all_gene_ids = set()
sample_size = min(1000, len(train_dataset))  # Sample to avoid memory issues

for i in range(sample_size):
    input_ids = train_dataset[i]['input_ids']
    all_gene_ids.update(input_ids)

# Remove special tokens (padding, cls, sep, etc.)
# Geneformer uses token IDs, where gene IDs typically start from a certain number
# Special tokens are usually < 10, so filter those out
genes_to_perturb = [gene_id for gene_id in all_gene_ids if gene_id >= 10]

print(f"Found {len(genes_to_perturb)} unique genes in the dataset")
print(f"First 10 gene IDs: {sorted(genes_to_perturb)[:10]}")

# For very large gene sets, you might want to limit the analysis
# Uncomment the next lines to limit to most frequently expressed genes
# from collections import Counter
# gene_frequency = Counter()
# for i in range(sample_size):
#     gene_frequency.update(train_dataset[i]['input_ids'])
# 
# # Get top N most frequent genes (excluding special tokens)
# top_genes = [gene_id for gene_id, freq in gene_frequency.most_common() if gene_id >= 10]
# genes_to_perturb = top_genes[:500]  # Limit to top 500 genes
# print(f"Limited analysis to top {len(genes_to_perturb)} most expressed genes")

# 4. Initialize the InSilicoPerturber
print("Initializing InSilicoPerturber...")
isp = InSilicoPerturber(
    perturb_type="delete",  # Options: "delete", "overexpress", "inhibit"
    perturb_rank_shift=None,  # For overexpress/inhibit, how much to shift
    genes_to_perturb=genes_to_perturb,
    combos=0,  # 0 for single gene perturbations, >0 for combinations
    anchor_gene=None,  # Gene to anchor combinations to
    model_type="Pretrained",  # Since we're using pretrained model
    num_cls=0,  # Number of classes (0 for pretrained)
    emb_mode="cell",  # Embedding mode
    cell_emb_style="mean_pool",  # How to pool cell embeddings
    filter_data={"cell_type": ["cardiomyocyte"]},  # Filter for specific cell types
    cell_states_to_model={"state_key": "disease", "start_state": "dcm", "goal_state": "nf"},
    max_ncells=500,  # Reduced for comprehensive screen
    emb_layer=-1,  # Which layer to extract embeddings from
    forward_batch_size=50,  # Reduced batch size for comprehensive analysis
    nproc=4  # Number of processes
)

# 5. Run the perturbation analysis
print("Running in silico perturbation analysis...")
print("This may take several minutes depending on your hardware...")

# Prepare output directory
output_dir = "./geneformer_perturbation_results"
os.makedirs(output_dir, exist_ok=True)

try:
    # Run perturbation
    perturbation_stats = isp.perturb_data(
        model_directory=model_name,
        input_data_file=train_dataset,  # Use the loaded dataset
        output_directory=output_dir,
        output_prefix="dcm_perturbation"
    )
    
    print("Perturbation analysis completed successfully!")
    
    # 6. Load and analyze results - Comprehensive Analysis
    print("Loading and analyzing comprehensive perturbation results...")
    
    # Load the perturbation results
    results_file = os.path.join(output_dir, "dcm_perturbation_stats.pickle")
    if os.path.exists(results_file):
        with open(results_file, 'rb') as f:
            results = pickle.load(f)
        
        print("Results loaded successfully!")
        print(f"Number of perturbations analyzed: {len(results)}")
        
        # Create comprehensive analysis of all gene perturbations
        gene_effects = []
        
        for gene_id, stats in results.items():
            effect_data = {
                'gene_id': gene_id,
                'mean_shift': stats.get('mean_shift', 0),
                'mean_shift_magnitude': abs(stats.get('mean_shift', 0)),
                'n_cells': stats.get('n_cells', 0),
                'std_shift': stats.get('std_shift', 0),
                'significance_score': abs(stats.get('mean_shift', 0)) / (stats.get('std_shift', 1) + 1e-8)  # Effect size / variability
            }
            gene_effects.append(effect_data)
        
        # Convert to DataFrame for easier analysis
        results_df = pd.DataFrame(gene_effects)
        
        # Sort by significance score (combination of effect size and consistency)
        results_df = results_df.sort_values('significance_score', ascending=False)
        
        print("\nTop 20 Most Impactful Gene Perturbations:")
        print("="*60)
        top_20 = results_df.head(20)
        for idx, row in top_20.iterrows():
            print(f"Gene ID: {row['gene_id']}")
            print(f"  Mean Shift: {row['mean_shift']:.4f}")
            print(f"  Magnitude: {row['mean_shift_magnitude']:.4f}")
            print(f"  Significance Score: {row['significance_score']:.4f}")
            print(f"  Cells Analyzed: {row['n_cells']}")
            print("-" * 40)
        
        # Save comprehensive results
        results_df.to_csv(os.path.join(output_dir, 'comprehensive_perturbation_results.csv'), index=False)
        print(f"\nFull results saved to: {os.path.join(output_dir, 'comprehensive_perturbation_results.csv')}")
    
    # 7. Enhanced Visualizations for comprehensive analysis
    if 'results_df' in locals():
        print("Creating comprehensive visualizations...")
        
        # 1. Top genes by significance score
        plt.figure(figsize=(15, 8))
        top_30 = results_df.head(30)
        plt.barh(range(len(top_30)), top_30['significance_score'])
        plt.xlabel('Significance Score (Effect Size / Variability)')
        plt.ylabel('Gene Rank')
        plt.title('Top 30 Most Impactful Gene Perturbations')
        plt.yticks(range(len(top_30)), [f"Gene_{gene_id}" for gene_id in top_30['gene_id']])
        plt.gca().invert_yaxis()
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, 'top_genes_significance.png'), dpi=300, bbox_inches='tight')
        plt.show()
        
        # 2. Effect size distribution
        plt.figure(figsize=(12, 6))
        plt.hist(results_df['mean_shift_magnitude'], bins=50, alpha=0.7, edgecolor='black')
        plt.xlabel('Perturbation Effect Magnitude')
        plt.ylabel('Number of Genes')
        plt.title('Distribution of Perturbation Effect Magnitudes')
        plt.axvline(results_df['mean_shift_magnitude'].mean(), color='red', linestyle='--', 
                   label=f'Mean: {results_df["mean_shift_magnitude"].mean():.4f}')
        plt.legend()
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, 'effect_distribution.png'), dpi=300, bbox_inches='tight')
        plt.show()
        
        # 3. Scatter plot: Effect size vs consistency
        plt.figure(figsize=(10, 8))
        scatter = plt.scatter(results_df['mean_shift_magnitude'], 
                            1/(results_df['std_shift'] + 1e-8),
                            c=results_df['significance_score'], 
                            alpha=0.6, cmap='viridis')
        plt.xlabel('Effect Magnitude')
        plt.ylabel('Consistency (1/Std Dev)')
        plt.title('Gene Perturbation Effects: Magnitude vs Consistency')
        plt.colorbar(scatter, label='Significance Score')
        
        # Annotate top 10 genes
        top_10 = results_df.head(10)
        for idx, row in top_10.iterrows():
            plt.annotate(f"Gene_{row['gene_id']}", 
                        (row['mean_shift_magnitude'], 1/(row['std_shift'] + 1e-8)),
                        xytext=(5, 5), textcoords='offset points', fontsize=8)
        
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, 'effect_vs_consistency.png'), dpi=300, bbox_inches='tight')
        plt.show()
        
        # 4. Summary statistics
        print("\nSUMMARY STATISTICS:")
        print("="*40)
        print(f"Total genes analyzed: {len(results_df)}")
        print(f"Mean effect magnitude: {results_df['mean_shift_magnitude'].mean():.4f}")
        print(f"Median effect magnitude: {results_df['mean_shift_magnitude'].median():.4f}")
        print(f"Top 1% threshold: {results_df['significance_score'].quantile(0.99):.4f}")
        print(f"Top 5% threshold: {results_df['significance_score'].quantile(0.95):.4f}")
        
        # Identify highly impactful genes (top 5%)
        top_5_percent = results_df[results_df['significance_score'] >= results_df['significance_score'].quantile(0.95)]
        print(f"\nTop 5% most impactful genes ({len(top_5_percent)} genes):")
        for idx, row in top_5_percent.iterrows():
            print(f"  Gene ID {row['gene_id']}: Significance Score = {row['significance_score']:.4f}")
        
        print(f"\nAll visualizations and results saved in: {output_dir}")

except Exception as e:
    print(f"Error during perturbation analysis: {str(e)}")
    print("This might be due to data format issues or model compatibility.")
    print("Please check that the dataset format matches Geneformer expectations.")

# 8. Code template for your own data
print("\n" + "="*60)
print("TEMPLATE FOR YOUR OWN DATA")
print("="*60)

template_code = '''
# Template for running perturbation analysis on your own data

# 1. Prepare your data in Geneformer format
# Your data should be a .dataset file with the following structure:
# - Each sample should have "input_ids" (tokenized gene expression)
# - Metadata columns for filtering (e.g., "cell_type", "disease_state")

# 2. Load your dataset
your_dataset = load_dataset("path/to/your/dataset.dataset")

# 3. Customize perturbation parameters for your analysis
your_genes_to_perturb = [
    "ENSG00000xxxxxx",  # Replace with your genes of interest
    "ENSG00000yyyyyy",
    # Add more genes as needed
]

# 4. Initialize perturber with your specific parameters
your_isp = InSilicoPerturber(
    perturb_type="delete",  # or "overexpress", "inhibit"
    genes_to_perturb=your_genes_to_perturb,
    combos=0,  # Set to >0 for gene combinations
    model_type="Pretrained",
    emb_mode="cell",
    cell_emb_style="mean_pool",
    filter_data={"cell_type": ["your_cell_type"]},  # Customize filtering
    cell_states_to_model={"state_key": "your_condition", 
                         "start_state": "disease", 
                         "goal_state": "control"},  # Customize states
    max_ncells=1000,
    forward_batch_size=100,
    nproc=4
)

# 5. Run analysis on your data
your_results = your_isp.perturb_data(
    model_directory="ctheodoris/Geneformer",
    input_data_file=your_dataset["train"],
    output_directory="./your_perturbation_results",
    output_prefix="your_analysis"
)
'''

print(template_code)

print("="*60)
print("NEXT STEPS:")
print("1. Ensure your data is in the correct Geneformer format (.dataset file)")
print("2. Identify genes of interest for perturbation")
print("3. Customize the filter_data and cell_states_to_model parameters")
print("4. Run the analysis with your data")
print("5. Analyze results using the visualization code above")
