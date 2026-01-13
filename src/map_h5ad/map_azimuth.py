import scanpy as sc
import numpy as np
import umap
from sklearn.decomposition import PCA

# -----------------------
# 0) Load
# -----------------------
ref = sc.read_h5ad("reference.h5ad")
qry = sc.read_h5ad("query.h5ad")

# -----------------------
# 1) Match genes (critical)
# -----------------------
common = ref.var_names.intersection(qry.var_names)
ref = ref[:, common].copy()
qry = qry[:, common].copy()
qry = qry[:, ref.var_names].copy()  # same order

# -----------------------
# 2) Preprocess reference (choose ONE consistent recipe)
#    If ref.X is counts: normalize_total + log1p.
#    If ref.X already log-normalized: you can skip normalize/log.
# -----------------------
# (Counts-like recipe)
sc.pp.normalize_total(ref, target_sum=1e4)
sc.pp.log1p(ref)

# Pick HVGs on reference and reuse for query (helps stability)
sc.pp.highly_variable_genes(ref, n_top_genes=3000, flavor="seurat_v3")
hvg = ref.var["highly_variable"].values
ref = ref[:, hvg].copy()
qry = qry[:, hvg].copy()

# Scale reference and capture mean/std for applying to query
sc.pp.scale(ref, max_value=10)
ref_mean = ref.X.mean(axis=0)
ref_std  = ref.X.std(axis=0)
ref_std[ref_std == 0] = 1.0

# -----------------------
# 3) Fit PCA on reference, then UMAP on reference PCA
# -----------------------
n_pcs = 50
pca = PCA(n_components=n_pcs, svd_solver="randomized", random_state=0)
ref_pca = pca.fit_transform(ref.X)

umap_model = umap.UMAP(
    n_neighbors=30,
    min_dist=0.3,
    n_components=2,
    metric="euclidean",
    random_state=0
)
ref_umap = umap_model.fit_transform(ref_pca)

# Save ref UMAP into ref (optional)
ref.obsm["X_pca_ref"] = ref_pca
ref.obsm["X_umap_ref"] = ref_umap

# -----------------------
# 4) Preprocess query in the SAME way, then project + transform
# -----------------------
sc.pp.normalize_total(qry, target_sum=1e4)
sc.pp.log1p(qry)

# Apply reference scaling (don’t refit!)
# If qry.X is sparse, convert just once
Xq = qry.X
if not isinstance(Xq, np.ndarray):
    Xq = Xq.toarray()

Xq = (Xq - ref_mean) / ref_std
Xq = np.clip(Xq, -10, 10)

qry_pca = pca.transform(Xq)

# This is the key: anchor onto the reference UMAP
qry_umap = umap_model.transform(qry_pca)

qry.obsm["X_pca_ref"] = qry_pca
qry.obsm["X_umap_ref"] = qry_umap

qry.write_h5ad("query_anchored_to_ref_umap.h5ad")
print("Done. Wrote query_anchored_to_ref_umap.h5ad with obsm['X_umap_ref']")
