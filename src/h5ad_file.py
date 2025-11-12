import anndata as ad
import h5py
import json
import re


def print_h5_structure(name, obj):
    if isinstance(obj, h5py.Dataset):
        print(f"📄 Dataset: {name}, shape={obj.shape}, dtype={obj.dtype}")
    elif isinstance(obj, h5py.Group):
        print(f"📂 Group: {name}")



# read the .h5ad file
adata = ad.read_h5ad(r"C:\Users\Orangeeee\Downloads\secondary_analysis.h5ad")


genes = set(adata.var.index)

with open(r"E:\UIC_PHD\Python_projects\parquet_read\src\matrisome_list_ENSG.json", "r", encoding="utf-8") as f:
    data = json.load(f)

def get_leaf_values(data):
    """Recursively collect all most-nested (leaf) values from a dict/list structure."""
    values = set()

    if isinstance(data, dict):
        for v in data.values():
            values |= get_leaf_values(v)
    elif isinstance(data, list):
        for v in data:
            values |= get_leaf_values(v)
    else:
        # Base case: not a dict or list -> a leaf value
        values.add(data)

    return values


cleaned = {re.sub(r"\.\d+$", "", g) for g in genes}
genes_matrisome = get_leaf_values(data)

common = cleaned & genes_matrisome