import scanpy as sc

b = sc.read_h5ad(r"E:\UIC_PHD\PythonProject\1-28-2026-Matriviz_h5ad_mapping_Tabula\data\bladder_ref.h5ad")
k = sc.read_h5ad(r"E:\UIC_PHD\PythonProject\1-28-2026-Matriviz_h5ad_mapping_Tabula\data\knee_ref.h5ad")
li = sc.read_h5ad(r"E:\UIC_PHD\PythonProject\1-28-2026-Matriviz_h5ad_mapping_Tabula\data\large_intestine_ref.h5ad")
ln = sc.read_h5ad(r"E:\UIC_PHD\PythonProject\1-28-2026-Matriviz_h5ad_mapping_Tabula\data\lymph_node_ref.h5ad")
si = sc.read_h5ad(r"E:\UIC_PHD\PythonProject\1-28-2026-Matriviz_h5ad_mapping_Tabula\data\small_intestine_ref.h5ad")
s = sc.read_h5ad(r"E:\UIC_PHD\PythonProject\1-28-2026-Matriviz_h5ad_mapping_Tabula\data\spleen_ref.h5ad")

import numpy as np
import pandas as pd


def check_data_transformation(df):
    """
    Heuristically identifies if a dataframe/matrix is:
    1. Raw Counts (Integers)
    2. Log-transformed (Small positive floats)
    3. Z-score normalized (Centered around 0, includes negatives)
    """
    # Convert to numpy for faster stats if it's a DataFrame
    data = df.values if isinstance(df, pd.DataFrame) else df

    # 1. Check for Raw Counts
    # If all values are non-negative integers
    if np.all(data >= 0) and np.all(np.mod(data, 1) == 0):
        return "Raw Counts"

    # 2. Check for Z-score Normalized
    # If there are negative values and the mean is very close to 0
    if np.any(data < 0):
        mean_val = np.mean(data)
        if np.isclose(mean_val, 0, atol=1e-2):
            return "Z-score Normalized"
        return "Unknown (Contains negatives but not centered)"

    # 3. Check for Log-transformed
    # If all values are non-negative and the max value is relatively small (usually < 20 for scRNA-seq)
    if np.all(data >= 0):
        max_val = np.max(data)
        if max_val < 25:  # Typical threshold for log1p(CPM/TPM)
            return "Log-transformed"
        else:
            return "Normalized but not Log-transformed (e.g., CPM/TPM)"

    return "Unknown Format"