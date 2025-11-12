import anndata as ad
import h5py


def print_h5_structure(name, obj):
    if isinstance(obj, h5py.Dataset):
        print(f"📄 Dataset: {name}, shape={obj.shape}, dtype={obj.dtype}")
    elif isinstance(obj, h5py.Group):
        print(f"📂 Group: {name}")

path = r"C:\Users\Orangeeee\Downloads\expression_matrices.h5"
with h5py.File(path, "r") as f:
    print(f"=== File: {path} ===")
    f.visititems(print_h5_structure)


# read the .h5ad file
# adata = ad.read_h5ad(r"C:\Users\Orangeeee\Downloads\expr.h5ad")
