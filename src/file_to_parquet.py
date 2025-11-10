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


# =========================
# Worker Thread (Background Conversion)
# =========================
class ConverterThread(QThread):
    progress = pyqtSignal(str)
    finished = pyqtSignal(bool, str)

    def __init__(self, input_file, tissue, output_dir):
        super().__init__()
        self.input_file = input_file
        self.tissue = tissue
        self.output_dir = output_dir

    def run(self):
        try:
            self.progress.emit("Loading data...")
            os.makedirs(self.output_dir, exist_ok=True)
            if self.input_file.endswith(".h5ad"):
                adata = ad.read_h5ad(self.input_file)
            else:
                raise ValueError("Supported input: .h5ad")

            # Extract expression
            expr = adata.X.toarray() if not isinstance(adata.X, np.ndarray) else adata.X
            expr_df = pd.DataFrame(expr, index=adata.obs_names, columns=adata.var_names)

            # UMAP
            if "X_umap" not in adata.obsm:
                raise ValueError("No UMAP coordinates found in object.")
            umap_df = pd.DataFrame(
                adata.obsm["X_umap"], columns=["UMAP_1", "UMAP_2"], index=adata.obs_names
            )
            umap_df["index"] = umap_df.index

            # Metadata
            meta_df = adata.obs.copy()
            meta_df["index"] = meta_df.index
            celltype_cols = [c for c in meta_df.columns if "celltype" in c.lower()]
            label_col = celltype_cols[-1] if celltype_cols else None

            self.progress.emit("Merging data...")
            merged = umap_df.join(expr_df, how="inner")
            merged["All_Genes"] = merged.iloc[:, 2:].mean(axis=1, skipna=True)
            merged = merged.loc[:, ~merged.columns.str.contains("celltype", case=False)]

            # Expression parquet
            expr_filename = f"Azimuth_{self.tissue}_ADT_reference.parquet"
            expr_path = Path(self.output_dir) / expr_filename
            pq.write_table(pa.Table.from_pandas(merged.reset_index(drop=True)), expr_path, compression="zstd")

            self.progress.emit("Computing centroids...")
            if label_col:
                meta_sub = meta_df[[label_col, "index"]].merge(umap_df, on="index", how="left")
                centroids = (
                    meta_sub.groupby(label_col)[["UMAP_1", "UMAP_2"]]
                    .median()
                    .reset_index()
                    .rename(columns={label_col: "Type", "UMAP_1": "cen_x", "UMAP_2": "cen_y"})
                )
            else:
                centroids = pd.DataFrame(columns=["Type", "cen_x", "cen_y"])

            centroid_filename = f"Azimuth_{self.tissue}_centroid.parquet"
            centroid_path = Path(self.output_dir) / centroid_filename
            pq.write_table(pa.Table.from_pandas(centroids), centroid_path, compression="zstd")

            self.progress.emit("Writing metadata JSON...")
            meta_json = {
                "fileType": "matriviz",
                "version": "1.1.0",
                "category_name": self.tissue.lower(),
                "category_description": self.tissue,
                "parquet_file": expr_filename,
                "category_file": f"{self.tissue.lower()}_category.json",
                "centroid_file": centroid_filename
            }
            json_path = Path(self.output_dir) / f"{self.tissue.lower()}_matriviz_meta.json"
            with open(json_path, "w", encoding="utf-8") as f:
                json.dump(meta_json, f, indent=4)

            self.finished.emit(True, f"✅ Conversion completed!\n\nOutput: {self.output_dir}")
        except Exception as e:
            self.finished.emit(False, f"❌ Error: {str(e)}")


# =========================
# Main GUI Window
# =========================
class AzimuthApp(QWidget):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Azimuth → MatriViz Converter")
        self.setFixedSize(480, 300)

        layout = QVBoxLayout()

        # Input file
        self.input_label = QLabel("Input .h5ad file:")
        self.input_line = QLineEdit()
        self.input_btn = QPushButton("Browse")
        self.input_btn.clicked.connect(self.select_input)
        hl1 = QHBoxLayout()
        hl1.addWidget(self.input_line)
        hl1.addWidget(self.input_btn)

        # Tissue name
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

        # Progress bar
        self.progress_bar = QProgressBar()
        self.progress_bar.setRange(0, 0)
        self.progress_bar.setVisible(False)

        # Convert button
        self.convert_btn = QPushButton("Convert")
        self.convert_btn.clicked.connect(self.start_conversion)

        layout.addWidget(self.input_label)
        layout.addLayout(hl1)
        layout.addWidget(self.tissue_label)
        layout.addWidget(self.tissue_line)
        layout.addWidget(self.output_label)
        layout.addLayout(hl2)
        layout.addWidget(self.convert_btn)
        layout.addWidget(self.progress_bar)
        self.setLayout(layout)

    def select_input(self):
        file, _ = QFileDialog.getOpenFileName(self, "Select .h5ad File", "", "H5AD files (*.h5ad)")
        if file:
            self.input_line.setText(file)

    def select_output(self):
        folder = QFileDialog.getExistingDirectory(self, "Select Output Directory")
        if folder:
            self.output_line.setText(folder)

    def start_conversion(self):
        input_file = self.input_line.text().strip()
        tissue = self.tissue_line.text().strip()
        output_dir = self.output_line.text().strip()

        if not all([input_file, tissue, output_dir]):
            QMessageBox.warning(self, "Missing Info", "Please fill all fields.")
            return

        self.convert_btn.setEnabled(False)
        self.progress_bar.setVisible(True)

        self.worker = ConverterThread(input_file, tissue, output_dir)
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
# Run App
# =========================
if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = AzimuthApp()
    window.show()
    sys.exit(app.exec())
