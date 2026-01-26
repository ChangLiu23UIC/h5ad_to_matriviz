import scanpy as sc
import pandas as pd
import os

# Define your input directory and output directory
input_dir = r"F:\12_30_2025_test_azimuth_annotation\bladder_mapped"
output_dir = r"F:\12_30_2025_test_azimuth_annotation\comparison_results"

# Create output directory if it doesn't exist
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# Iterate through all h5ad files
for filename in os.listdir(input_dir):
    if filename.endswith(".h5ad"):
        file_path = os.path.join(input_dir, filename)
        print(f"Processing: {filename}")

        try:
            # Read the file (backed='r' is faster if you only need .obs)
            adata = sc.read_h5ad(file_path, backed='r')

            # Check if both columns exist to avoid KeyErrors
            cols_to_extract = ['azimuth_label', 'azimuth_label_tabula']
            if all(col in adata.obs.columns for col in cols_to_extract):
                # Extract columns
                df_labels = adata.obs[cols_to_extract].copy()

                # Create output filename (e.g., HBM236.JPVT.769_compare.csv)
                base_name = os.path.splitext(filename)[0]
                output_file = os.path.join(output_dir, f"{base_name}_compare.csv")

                # Save to CSV
                df_labels.to_csv(output_file)
            else:
                print(f"Skipping {filename}: Required columns not found.")

        except Exception as e:
            print(f"Error processing {filename}: {e}")

print("Batch processing complete.")