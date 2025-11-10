import pandas as pd

# Define a helper function
def add_all_genes_column(df):
    # exclude specific columns and compute mean across all remaining columns (axis=1)
    gene_cols = [c for c in df.columns if c not in ['umap_1', 'umap_2', 'index']]
    df['All_Genes'] = df[gene_cols].mean(axis=1)
    return df

# Load datasets
cortex = pd.read_parquet(r"E:\UIC_PHD\Matriviz_project\MatriViz-Data\cortex_v1.parquet")
heart = pd.read_parquet(r"E:\UIC_PHD\Matriviz_project\MatriViz-Data\heart_v1_ss.parquet")
kidney = pd.read_parquet(r"E:\UIC_PHD\Matriviz_project\MatriViz-Data\kidney_v1_ss.parquet")
pancrea = pd.read_parquet(r"E:\UIC_PHD\Matriviz_project\MatriViz-Data\pancrea_v1_ss.parquet")
pbmc = pd.read_parquet(r"E:\UIC_PHD\Matriviz_project\MatriViz-Data\pbmc_v1_ss.parquet")
tonsil = pd.read_parquet(r"E:\UIC_PHD\Matriviz_project\MatriViz-Data\tonsil_v1_ss.parquet")


test_centroid = pd.read_parquet(r"E:\UIC_PHD\Azimuth_11072025\PBMC\Azimuth_PBMC_centroid.parquet")
test_cortex = pd.read_parquet(r"E:\UIC_PHD\Azimuth_11072025\PBMC\Azimuth_PBMC_ADT_reference.parquet")

centroid = pd.read_parquet(r"E:\UIC_PHD\Matriviz_project\MatriViz-Data\cortex_centroid_v1.parquet")

# # Apply to each dataframe
# datasets = [cortex, heart, kidney, pancrea, pbmc, tonsil]
# datasets = [add_all_genes_column(df) for df in datasets]
#
# # Unpack back if you want to keep variable names
# cortex, heart, kidney, pancrea, pbmc, tonsil = datasets
