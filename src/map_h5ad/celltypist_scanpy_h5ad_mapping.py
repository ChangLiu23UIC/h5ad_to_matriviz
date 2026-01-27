import scanpy as sc
import os
import glob
import numpy as np
import pandas as pd
import celltypist
from scipy.sparse import csr_matrix
from tqdm import tqdm

# --- Path Configuration ---
ref_path = r"E:\UIC_PHD\Python_projects\parquet_read\src\map_h5ad\refs\knee_ref.h5ad"
input_dir = r"F:\12_30_2025_test_azimuth_annotation"
output_dir = r"E:\UIC_PHD\PythonProject\12-4-2025-Matriviz\raw\Mapped_knee_h5ad"
model_storage = os.path.join(os.path.dirname(ref_path), "knee_model.pkl")

SEED = 42
np.random.seed(SEED)

if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# --- 1. 准备 Reference 数据 ---
print("Loading and preprocessing reference data...")
ref_raw = sc.read_h5ad(ref_path)
ref_raw.var_names = ref_raw.var_names.str.split('.').str[0]

# --- CRITICAL FIX: Normalize reference for CellTypist ---
# CellTypist strictly expects log1p(count per 10k)
# If your reference already has decimals, check if it's already log1p
# If it is raw counts, run these:
sc.pp.normalize_total(ref_raw, target_sum=1e4)
sc.pp.log1p(ref_raw)

print("Training custom CellTypist model...")
# Now the matrix will satisfy the 'expect log1p normalized expression' check
knee_model = celltypist.train(
    ref_raw,
    labels='cell_type',
    genes=ref_raw.var_names.tolist(),
    check_normalization=True # This will now pass
)
knee_model.write(model_storage)

# The rest of your Ingest preprocessing (PCA/Neighbors) follows...
sc.pp.pca(ref_raw, random_state=SEED)
sc.pp.neighbors(ref_raw, n_neighbors=14, random_state=SEED)
sc.tl.umap(ref_raw, random_state=SEED)

# --- 2. Processing Loop ---
query_files = glob.glob(os.path.join(input_dir, "*.h5ad"))
pbar = tqdm(query_files, desc="Processing Samples", unit="sample")

for query_path in pbar:
    file_name = os.path.basename(query_path)
    pbar.set_description(f"Processing: {file_name}")
    save_path = os.path.join(output_dir, file_name)

    try:
        # A. Read Query
        adata = sc.read_h5ad(query_path)
        adata.var_names = adata.var_names.str.split('.').str[0]

        # B. Annotation via CellTypist (Uses Floats)
        # This provides the high-accuracy label you need
        predictions = celltypist.annotate(adata, model=model_storage, majority_voting=True)
        adata.obs['celltypist_label'] = predictions.predicted_labels['predicted_labels']
        adata.obs['celltypist_conf'] = predictions.probability_table.max(axis=1).values

        # C. UMAP Anchoring via Ingest
        # We handle alignment here similarly to your original logic but for mapping
        common_genes = adata.var_names.intersection(ref_raw.var_names)

        # Create a temporary alignment object
        n_obs = adata.n_obs
        new_X = csr_matrix((n_obs, ref_raw.n_vars))
        adata_map = sc.AnnData(X=new_X, obs=adata.obs.copy(), var=ref_raw.var.copy())

        if "unscaled" in adata.layers:
            source_data = adata[:, common_genes].layers["unscaled"]
            var_index_map = {gene: i for i, gene in enumerate(ref_raw.var_names)}
            target_indices = [var_index_map[gene] for gene in common_genes]

            temp_X = adata_map.X.tolil()
            temp_X[:, target_indices] = source_data
            adata_map.X = temp_X.tocsr()

            sc.pp.normalize_total(adata_map, target_sum=1e4)
            sc.pp.log1p(adata_map)

            # Execute Ingest for UMAP projection
            sc.tl.ingest(adata_map, ref_raw, obs='cell_type')

            # Sync mapping results back to original adata
            adata.obs['azimuth_label'] = adata_map.obs['cell_type']
            adata.obsm['X_umap_proj'] = adata_map.obsm['X_umap']
        else:
            tqdm.write(f"Warning: 'unscaled' layer not found in {file_name}. Ingest skipped.")

        # D. Save Result
        adata.write(save_path)

    except Exception as e:
        tqdm.write(f"Error processing {file_name}: {e}")

print("\nAll tasks completed! Random state fixed at 42.")