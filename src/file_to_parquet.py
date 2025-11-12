import sys, os, json
import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import anndata as ad
from pathlib import Path
from PyQt6.QtWidgets import (
    QApplication, QWidget, QLabel, QPushButton, QLineEdit, QFileDialog,
    QVBoxLayout, QHBoxLayout, QMessageBox, QProgressBar
)
from PyQt6.QtCore import Qt, QThread, pyqtSignal
import re


# =========================
# Worker Thread
# =========================
class ConverterThread(QThread):
    progress = pyqtSignal(str)
    finished = pyqtSignal(bool, str)

    def __init__(self, input_file, tissue, output_dir, json_template=None):
        super().__init__()
        self.input_file = str(input_file)
        self.tissue = str(tissue).strip()
        self.output_dir = str(output_dir)
        self.json_template = json_template

    def run(self):
        try:
            self.progress.emit("Loading data...")
            os.makedirs(self.output_dir, exist_ok=True)

            # --- Load or convert input ---
            if self.input_file.endswith(".h5ad"):
                adata = ad.read_h5ad(self.input_file)
            elif self.input_file.endswith((".h5", ".h5seurat")):
                self.progress.emit("Converting .h5seurat → .h5ad ...")
                temp_h5ad = os.path.join(self.output_dir, "temp_converted.h5ad")
                import subprocess
                r_script = f"""
                suppressMessages(library(SeuratDisk));
                Convert("{self.input_file}", dest = "h5ad", overwrite = TRUE);
                """
                result = subprocess.run(["Rscript", "-e", r_script], capture_output=True, text=True)
                if result.returncode != 0:
                    raise RuntimeError(f"R conversion failed:\n{result.stderr}")
                adata = ad.read_h5ad(temp_h5ad)
            else:
                raise ValueError("Supported input: .h5ad or .h5/.h5seurat")

            # --- Expression matrix ---
            self.progress.emit("Extracting expression matrix...")
            expr = adata.X.toarray() if not isinstance(adata.X, np.ndarray) else adata.X
            expr_df = pd.DataFrame(expr, index=adata.obs_names, columns=adata.var_names)
            expr_df.columns = expr_df.columns.astype(str)

            # --- Add All_Genes pseudo-gene ---
            self.progress.emit("Calculating All_Genes average expression...")
            # Calculate mean expression across all genes for each cell
            all_genes_mean = expr_df.mean(axis=1)
            # Add as a new column to the expression dataframe
            expr_df["All_Genes"] = all_genes_mean

            # --- UMAP coordinates ---
            if "X_umap" not in adata.obsm:
                raise ValueError("No UMAP coordinates found in object.")
            umap_df = pd.DataFrame(
                adata.obsm["X_umap"], columns=["umap_1", "umap_2"], index=adata.obs_names
            )
            umap_df["index"] = umap_df.index.astype(str)

            # --- Metadata ---
            meta_df = adata.obs.copy()
            meta_df["index"] = meta_df.index
            celltype_cols = [c for c in meta_df.columns if "celltype" in c.lower()]
            label_col = celltype_cols[-1] if celltype_cols else None

            # --- Merge expression + UMAP ---
            self.progress.emit("Merging data...")
            merged = umap_df.join(expr_df, how="inner")
            # Remove any metadata columns that might interfere with gene selection
            merged = merged.loc[:, ~merged.columns.str.contains("celltype", case=False)]
            merged.columns = merged.columns.astype(str)

            # Ensure index column is first and properly formatted
            merged = merged.reset_index(drop=False)
            merged = merged.rename(columns={'index': 'original_index'})
            merged['index'] = merged['original_index'].astype(str)
            merged = merged.drop('original_index', axis=1)

            # Reorder columns to put index first
            cols = ['index', 'umap_1', 'umap_2'] + [col for col in merged.columns if col not in ['index', 'umap_1', 'umap_2']]
            merged = merged[cols]

            # =============================
            # 🔹 File naming
            # =============================
            base_name = str(self.tissue).lower().replace(" ", "_")
            expr_filename = f"{base_name}_v1.parquet"
            centroid_filename = f"{base_name}_centroid_v1.parquet"
            category_filename = f"{base_name}_v1_category.json"
            overview_filename = f"{base_name}.json"

            expr_path = Path(self.output_dir) / expr_filename
            centroid_path = Path(self.output_dir) / centroid_filename
            category_path = Path(self.output_dir) / category_filename
            overview_path = Path(self.output_dir) / overview_filename

            # --- Write expression parquet ---
            self.progress.emit("Writing expression parquet...")
            # Ensure index column is preserved as string type
            merged['index'] = merged['index'].astype(str)
            pq.write_table(pa.Table.from_pandas(merged), expr_path, compression="zstd")

            # --- Compute centroids ---
            self.progress.emit("Computing centroids...")
            if label_col:
                meta_sub = meta_df[[label_col, "index"]].merge(umap_df, on="index", how="left")
                centroids = (
                    meta_sub.groupby(label_col)[["umap_1", "umap_2"]]
                    .median()
                    .reset_index()
                    .rename(columns={label_col: "Type", "umap_1": "cen_x", "umap_2": "cen_y"})
                )
            else:
                centroids = pd.DataFrame(columns=["Type", "cen_x", "cen_y"])
            pq.write_table(pa.Table.from_pandas(centroids), centroid_path, compression="zstd")

            # --- Category JSON (flattened for MatriViz) ---
            self.progress.emit("Writing category JSON...")
            if self.json_template and os.path.exists(self.json_template):
                with open(self.json_template, "r", encoding="utf-8") as f:
                    template_json = json.load(f)
                # Get available genes from H5AD file for filtering
                available_genes = set(adata.var_names.astype(str))
                # Flatten nested structure for MatriViz compatibility and filter genes
                category_json = self._flatten_category_structure(template_json, available_genes)
            else:
                category_json = {
                    "ECM Glycoproteins": [],
                    "Secreted Factors": [],
                    "Collagens": [],
                    "ECM Regulators": [],
                    "Proteoglycans": [],
                    "ECM-affiliated Proteins": []
                }

            with open(category_path, "w", encoding="utf-8") as f:
                json.dump(category_json, f, indent=4)

            # --- Overview JSON ---
            self.progress.emit("Writing overview JSON...")
            overview_json = {
                "fileType": "matriviz",
                "version": "0.0.1",
                "category_name": base_name,
                "category_description": str(self.tissue),
                "parquet_file": expr_filename,
                "category_file": category_filename,
                "centroid_file": centroid_filename
            }

            with open(overview_path, "w", encoding="utf-8") as f:
                json.dump(overview_json, f, indent=4)

            self.finished.emit(True, f"✅ Conversion completed!\n\nOutput: {self.output_dir}")

        except Exception as e:
            import traceback
            tb = traceback.format_exc()
            self.finished.emit(False, f"❌ Error: {str(e)}\n\nTraceback:\n{tb}")

    def _flatten_category_structure(self, category_data, available_genes=None):
        """
        Flatten nested category structure to MatriViz-compatible format
        and filter genes based on what's available in the H5AD file

        Updated to:
        - Match ENSG IDs ignoring version suffix
        - Keep versioned ENSG IDs in the output
        - Use only secondary-level categories as keys

        Args:
            category_data: Nested category structure
            available_genes: Set of gene names available in the H5AD file (with version suffixes)

        Returns:
            dict: Flattened structure with secondary_category -> [versioned_gene_list]
        """
        flattened = {}

        # Create mapping from unversioned to versioned gene IDs
        versioned_gene_map = {}
        if available_genes:
            for versioned_gene in available_genes:
                # Remove version suffix to get base ENSG ID
                base_gene = re.sub(r"\.\d+$", "", versioned_gene)
                versioned_gene_map[base_gene] = versioned_gene

        def extract_categories(data, current_category=""):
            if isinstance(data, dict):
                for key, value in data.items():
                    # Track the current category level
                    if current_category == "":
                        # First level - primary category (e.g., "Matrisome-associated")
                        extract_categories(value, key)
                    elif current_category != "" and isinstance(value, dict):
                        # Second level - secondary category (e.g., "ECM-affiliated Proteins")
                        # This is the level we want to use as keys
                        extract_categories(value, key)
                    elif isinstance(value, list):
                        # Found a gene list - this is a subcategory
                        if available_genes:
                            # Filter genes to only include those available in H5AD
                            # Match by base ENSG ID (without version)
                            filtered_genes = []
                            for base_gene in value:
                                if base_gene in versioned_gene_map:
                                    # Use the versioned gene ID from H5AD
                                    filtered_genes.append(versioned_gene_map[base_gene])

                            if filtered_genes:  # Only add category if it has genes
                                # Use the secondary category as the key
                                if current_category not in flattened:
                                    flattened[current_category] = []
                                flattened[current_category].extend(filtered_genes)
                        else:
                            # If no available_genes provided, include all genes
                            if current_category not in flattened:
                                flattened[current_category] = []
                            flattened[current_category].extend(value)
                    else:
                        # Continue traversing nested structure
                        extract_categories(value, current_category)
            elif isinstance(data, list):
                # Direct gene list at current category level
                if current_category:
                    if available_genes:
                        # Filter genes to only include those available in H5AD
                        filtered_genes = []
                        for base_gene in data:
                            if base_gene in versioned_gene_map:
                                # Use the versioned gene ID from H5AD
                                filtered_genes.append(versioned_gene_map[base_gene])

                        if filtered_genes:  # Only add category if it has genes
                            if current_category not in flattened:
                                flattened[current_category] = []
                            flattened[current_category].extend(filtered_genes)
                    else:
                        # If no available_genes provided, include all genes
                        if current_category not in flattened:
                            flattened[current_category] = []
                        flattened[current_category].extend(data)

        extract_categories(category_data)

        # Remove duplicates from each category
        for category in flattened:
            flattened[category] = list(set(flattened[category]))

        return flattened


# =========================
# GUI
# =========================
class AzimuthApp(QWidget):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Azimuth → MatriViz Converter")
        self.setFixedSize(480, 380)

        layout = QVBoxLayout()

        # Input file
        self.input_label = QLabel("Input .h5ad file:")
        self.input_line = QLineEdit()
        self.input_btn = QPushButton("Browse")
        self.input_btn.clicked.connect(self.select_input)
        hl1 = QHBoxLayout()
        hl1.addWidget(self.input_line)
        hl1.addWidget(self.input_btn)

        # Tissue
        self.tissue_label = QLabel("Tissue name:")
        self.tissue_line = QLineEdit()

        # Output dir
        self.output_label = QLabel("Output directory:")
        self.output_line = QLineEdit()
        self.output_btn = QPushButton("Browse")
        self.output_btn.clicked.connect(self.select_output)
        hl2 = QHBoxLayout()
        hl2.addWidget(self.output_line)
        hl2.addWidget(self.output_btn)

        # JSON template
        self.json_label = QLabel("Category JSON template (optional):")
        self.json_line = QLineEdit()
        self.json_btn = QPushButton("Browse")
        self.json_btn.clicked.connect(self.select_json)
        hl3 = QHBoxLayout()
        hl3.addWidget(self.json_line)
        hl3.addWidget(self.json_btn)

        # Progress bar + convert button
        self.progress_bar = QProgressBar()
        self.progress_bar.setRange(0, 0)
        self.progress_bar.setVisible(False)
        self.convert_btn = QPushButton("Convert")
        self.convert_btn.clicked.connect(self.start_conversion)

        layout.addWidget(self.input_label)
        layout.addLayout(hl1)
        layout.addWidget(self.tissue_label)
        layout.addWidget(self.tissue_line)
        layout.addWidget(self.output_label)
        layout.addLayout(hl2)
        layout.addWidget(self.json_label)
        layout.addLayout(hl3)
        layout.addWidget(self.convert_btn)
        layout.addWidget(self.progress_bar)
        self.setLayout(layout)

    def select_input(self):
        file, _ = QFileDialog.getOpenFileName(self, "Select Input File", "", "H5AD/H5Seurat files (*.h5ad *.h5 *.h5seurat)")
        if file:
            self.input_line.setText(file)

    def select_output(self):
        folder = QFileDialog.getExistingDirectory(self, "Select Output Directory")
        if folder:
            self.output_line.setText(folder)

    def select_json(self):
        file, _ = QFileDialog.getOpenFileName(self, "Select Category JSON Template", "", "JSON files (*.json)")
        if file:
            self.json_line.setText(file)

    def start_conversion(self):
        input_file = self.input_line.text().strip()
        tissue = self.tissue_line.text().strip()
        output_dir = self.output_line.text().strip()
        json_template = self.json_line.text().strip() or None

        if not all([input_file, tissue, output_dir]):
            QMessageBox.warning(self, "Missing Info", "Please fill all required fields.")
            return

        self.convert_btn.setEnabled(False)
        self.progress_bar.setVisible(True)

        self.worker = ConverterThread(input_file, tissue, output_dir, json_template)
        self.worker.progress.connect(self.show_progress)
        self.worker.finished.connect(self.finish_conversion)
        self.worker.start()

    def show_progress(self, msg):
        self.progress_bar.setFormat(msg)

    def finish_conversion(self, success, msg):
        self.convert_btn.setEnabled(True)
        self.progress_bar.setVisible(False)
        QMessageBox.information(self, "Result", msg)


# =========================
# Run
# =========================
if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = AzimuthApp()
    window.show()
    sys.exit(app.exec())
