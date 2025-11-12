# H5AD to MatriViz Converter

A Python tool for converting h5ad files (AnnData format) to MatriViz-compatible parquet format for visualization.

## Features

- Convert h5ad files to parquet format compatible with MatriViz
- Extract UMAP coordinates and gene expression data
- Generate category files and overview JSON
- GUI interface for easy conversion
- Support for h5seurat conversion (requires R)
- Version suffix matching for ENSG gene IDs
- All_Genes pseudo-gene for average expression visualization

## Installation

1. Install Python dependencies:
```bash
pip install -r ../requirements.txt
```

2. For h5seurat support, install R and SeuratDisk:
```R
install.packages("SeuratDisk")
```

## Usage

### GUI Interface

Run the GUI application:
```bash
python file_to_parquet.py
```

1. Select your input .h5ad file
2. Enter tissue name (used for file naming)
3. Select output directory
4. (Optional) Select category JSON template
5. Click "Convert"

## Output Files

After conversion, you'll get these files:

- `{tissue}_v1.parquet` - Main expression data with UMAP coordinates
- `{tissue}_centroid_v1.parquet` - Cell type centroids (if available)
- `{tissue}_v1_category.json` - Gene category definitions
- `{tissue}.json` - Overview configuration for MatriViz

## Requirements

- Python 3.8+
- h5ad file with UMAP coordinates in `adata.obsm['X_umap']`
- Gene expression matrix in `adata.X`

## File Structure

- `file_to_parquet.py` - Main GUI application
- `matrisome_list.json` - Gene category definitions (symbol names)
- `matrisome_list_ENSG.json` - Gene category definitions (ENSG IDs)

## License

This tool is provided for use with the MatriViz application.