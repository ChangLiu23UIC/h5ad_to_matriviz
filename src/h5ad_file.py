import anndata as ad

# read the .h5ad file
# adata = ad.read_h5ad(r"C:\Users\Orangeeee\Downloads\expr.h5ad")
adata = ad.read_h5ad(r"C:\Users\Orangeeee\Downloads\0.h5ad")

# view basic info
print(adata)
