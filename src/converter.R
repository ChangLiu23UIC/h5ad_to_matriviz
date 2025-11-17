library(tcltk)
library(Seurat)
library(dplyr)
library(arrow)
library(jsonlite)

# ================================
# 1. Flatten category structure for gene symbol matching
#    SYMBOL-ONLY, CASE-INSENSITIVE, grouped by L2 category
# ================================
flatten_category_structure_gene_names <- function(category_data, available_genes) {
  if (is.null(available_genes)) {
    stop("available_genes must not be NULL")
  }

  flattened <- list()
  available_upper <- toupper(available_genes)

  # Recursive traversal: track L1 / L2 / L3, but we group by L2
  walk_node <- function(node, L1 = "", L2 = "", L3 = "") {
    # Only lists exist in your JSON structure
    if (!is.list(node)) return(NULL)

    # Check if this node is a "gene list": list of scalar character values
    is_gene_list <- all(
      sapply(node, function(x) is.character(x) && length(x) == 1)
    )

    if (is_gene_list) {
      # ---- GENE LIST NODE ----
      genes_raw   <- unlist(node, use.names = FALSE)
      genes_upper <- toupper(genes_raw)

      # Case-insensitive matching to Seurat genes
      match_idx   <- available_upper %in% genes_upper
      matched     <- available_genes[match_idx]

      if (length(matched) > 0) {
        # Use L2 (secondary category) as the key, as in the debug script
        category_key <- L2
        if (!category_key %in% names(flattened)) {
          flattened[[category_key]] <<- character()
        }
        flattened[[category_key]] <<- unique(
          c(flattened[[category_key]], matched)
        )
      }
    } else {
      # ---- CATEGORY / SUBCATEGORY NODE ----
      if (length(node) == 0) return(NULL)

      for (nm in names(node)) {
        val <- node[[nm]]
        if (L1 == "") {
          # Top level: Matrisome-associated / Core matrisome
          walk_node(val, L1 = nm, L2 = "", L3 = "")
        } else if (L2 == "") {
          # Second level: ECM-affiliated Proteins / ECM Regulators / etc.
          walk_node(val, L1 = L1, L2 = nm, L3 = "")
        } else {
          # Third level: Syndecan / Galectin / Laminin / ...
          walk_node(val, L1 = L1, L2 = L2, L3 = nm)
        }
      }
    }
  }

  walk_node(category_data)
  return(flattened)
}

# ================================
# 2. Main converter
# ================================
convert_azimuth <- function(rds_path, tissue, output_dir, json_template = NULL) {
  cat("Loading reference...\n")
  ref <- readRDS(rds_path)

  # --- Ensure compatibility with new Seurat structure ---
  ref <- UpdateSeuratObject(ref)

  tissue_lower <- tolower(tissue)

  assays_available <- Assays(ref)
  if (!("integrated" %in% assays_available)) {
    stop("integrated assay not found in Seurat object.")
  }

  cat("Extracting integrated assay data...\n")
  expr_sparse <- GetAssayData(ref, assay = "integrated", layer = "data")
  # genes are rows (rownames(expr_sparse)), cells are columns

  # --- Select variable features to avoid huge memory load ---
  if (length(VariableFeatures(ref)) > 0) {
    genes_use <- VariableFeatures(ref)
    cat("Using variable features (", length(genes_use), " genes)...\n")
  } else {
    genes_use <- rownames(expr_sparse)[1:min(2000, nrow(expr_sparse))]
    cat("No variable features found — using first", length(genes_use), "genes.\n")
  }

  expr_sparse <- expr_sparse[genes_use, , drop = FALSE]
  available_genes <- rownames(expr_sparse)
  cat("Number of genes (after subsetting) used for matching:",
      length(available_genes), "\n")

  cat("Converting to data frame (subset)...\n")
  # Keep original orientation: genes as rows, cells as columns
  expr_df <- as.data.frame(as.matrix(expr_sparse))
  expr_df$gene <- rownames(expr_df)

  # --- Extract UMAP embeddings ---
  cat("Extracting UMAP...\n")
  umap_name <- if ("refUMAP" %in% names(ref@reductions)) "refUMAP" else "umap"
  umap <- as.data.frame(Embeddings(ref, umap_name))
  umap$index <- rownames(umap)
  colnames(umap) <- c("UMAP_1", "UMAP_2", "index")

  # --- Extract useful metadata ---
  meta_cols <- intersect(
    colnames(ref@meta.data),
    c("annotation.l1", "annotation.l2", "annotation.l3",
      "predicted.celltype.l1", "predicted.celltype.l2",
      "celltype.l1", "celltype.l2")
  )

  meta <- ref@meta.data
  meta$index <- rownames(meta)
  if (length(meta_cols) > 0) {
    meta <- meta[, c("index", meta_cols), drop = FALSE]
  }

  # --- Merge data (cells as rows, genes as columns) ---
  cat("Merging data (cells as rows, genes as columns)...\n")
  expr_t <- as.data.frame(t(expr_df[, !colnames(expr_df) %in% "gene"]))
  expr_t$index <- rownames(expr_t)

  merged <- dplyr::left_join(umap, expr_t, by = "index")
  if (length(meta_cols) > 0) merged <- dplyr::left_join(merged, meta, by = "index")

  # --- Add mean gene expression summary ---
  gene_cols <- setdiff(colnames(merged), c("index", "UMAP_1", "UMAP_2", meta_cols))
  merged$All_Genes <- rowMeans(merged[, gene_cols], na.rm = TRUE)
  merged <- merged[, !grepl("celltype", colnames(merged), ignore.case = TRUE)]

  # --- File naming ---
  expr_file     <- file.path(output_dir, paste0(tissue_lower, "_rna_v1.parquet"))
  centroid_file <- file.path(output_dir, paste0(tissue_lower, "_rna_centroid_v1.0.parquet"))
  category_file <- file.path(output_dir, paste0(tissue_lower, "_rna_v1_category.json"))
  meta_file     <- file.path(output_dir, paste0(tissue_lower, "_rna.json"))

  # --- Save expression parquet ---
  cat("Writing expression parquet (cells × genes + UMAP + metadata)...\n")
  write_parquet(merged, expr_file, compression = "zstd", version = "1.0")

  # --- Compute centroids by cell type ---
  label_col <- if ("annotation.l1" %in% meta_cols) {
    "annotation.l1"
  } else if ("annotation.l2" %in% meta_cols) {
    "annotation.l2"
  } else if ("annotation.l3" %in% meta_cols) {
    "annotation.l3"
  } else if ("celltype.l2" %in% meta_cols) {
    "celltype.l2"
  } else if ("predicted.celltype.l2" %in% meta_cols) {
    "predicted.celltype.l2"
  } else if (length(meta_cols) > 0) {
    meta_cols[1]
  } else {
    NULL
  }

  cat("Computing centroids...\n")
  if (!is.null(label_col)) {
    centroids <- ref@meta.data %>%
      dplyr::select(all_of(label_col)) %>%
      mutate(index = rownames(ref@meta.data)) %>%
      dplyr::left_join(umap, by = "index") %>%
      dplyr::group_by(.data[[label_col]]) %>%
      dplyr::summarise(
        cen_x = mean(UMAP_1, na.rm = TRUE),
        cen_y = mean(UMAP_2, na.rm = TRUE)
      ) %>%
      dplyr::arrange(.data[[label_col]]) %>%
      dplyr::rename(Type = .data[[label_col]])
  } else {
    centroids <- data.frame(
      Type = character(),
      cen_x = numeric(),
      cen_y = numeric()
    )
    cat("No cell type columns found - creating empty centroids file\n")
  }

  write_parquet(centroids, centroid_file, compression = "snappy", version = "1.0")

  # ==========================
  # Category JSON processing
  # ==========================
  cat("Processing category JSON (gene symbols only, case-insensitive)...\n")

  category_json <- list()

  if (!is.null(json_template) && file.exists(json_template)) {
    template_json <- fromJSON(json_template, simplifyVector = FALSE)

    # Use the SAME logic as the working debug inspector
    matrisome_categories <- flatten_category_structure_gene_names(
      category_data   = template_json,
      available_genes = available_genes
    )

    # Print summary to console so you can confirm it matches debug_summary
    if (length(matrisome_categories) > 0) {
      cat("Matched categories and gene counts:\n")
      for (nm in names(matrisome_categories)) {
        cat("  -", nm, ":", length(matrisome_categories[[nm]]), "genes\n")
      }
    } else {
      cat("No matrisome genes matched; will fall back to empty skeleton.\n")
    }

    # Add to category_json
    for (category_name in names(matrisome_categories)) {
      genes_vec <- matrisome_categories[[category_name]]
      if (length(genes_vec) > 0) {
        category_json[[category_name]] <- genes_vec
      }
    }
  }

  # If still nothing, use default empty structure
  if (length(category_json) == 0) {
    category_json <- list(
      "ECM Glycoproteins"      = list(),
      "Secreted Factors"       = list(),
      "Collagens"              = list(),
      "ECM Regulators"         = list(),
      "Proteoglycans"          = list(),
      "ECM-affiliated Proteins"= list()
    )
  }

  write_json(category_json, category_file, pretty = TRUE, auto_unbox = TRUE)

  # --- Metadata JSON ---
  meta_json <- list(
    fileType            = "matriviz",
    version             = "1.1.0",
    category_name       = tissue_lower,
    category_description= paste(tissue, "(RNA reference)"),
    parquet_file        = basename(expr_file),
    category_file       = basename(category_file),
    centroid_file       = basename(centroid_file)
  )
  write_json(meta_json, meta_file, pretty = TRUE, auto_unbox = TRUE)

  tkmessageBox(
    title   = "Success",
    message = paste("✅ RNA Conversion complete!\n\nSaved to:", output_dir)
  )
}

# =============================
# GUI
# =============================
win <- tktoplevel()
tkwm.title(win, "Azimuth → MatriViz Converter (RDS with Categories)")
tkwm.resizable(win, 0, 0)

lbl_input <- tklabel(win, text = "Input .Rds file:")
ent_input <- tkentry(win, width = 50)
btn_input <- tkbutton(win, text = "Browse", command = function() {
  path <- tclvalue(tkgetOpenFile(filetypes = "{{RDS Files} {.Rds}}"))
  if (path != "") tkinsert(ent_input, "end", path)
})

lbl_tissue <- tklabel(win, text = "Tissue name:")
ent_tissue <- tkentry(win, width = 50)

lbl_json <- tklabel(win, text = "Category JSON template (optional):")
ent_json <- tkentry(win, width = 50)
btn_json <- tkbutton(win, text = "Browse", command = function() {
  path <- tclvalue(tkgetOpenFile(filetypes = "{{JSON Files} {.json}}"))
  if (path != "") tkinsert(ent_json, "end", path)
})

lbl_out <- tklabel(win, text = "Output directory:")
ent_out <- tkentry(win, width = 50)
btn_out <- tkbutton(win, text = "Browse", command = function() {
  dir <- tclvalue(tkchooseDirectory())
  if (dir != "") tkinsert(ent_out, "end", dir)
})

btn_convert <- tkbutton(win, text = "Convert", command = function() {
  rds_path      <- tclvalue(tkget(ent_input))
  tissue        <- tclvalue(tkget(ent_tissue))
  json_template <- tclvalue(tkget(ent_json))
  out_dir       <- tclvalue(tkget(ent_out))
  if (rds_path == "" || tissue == "" || out_dir == "") {
    tkmessageBox(title = "Error", message = "Please fill in all required fields.")
  } else {
    if (json_template == "") json_template <- NULL
    convert_azimuth(rds_path, tissue, out_dir, json_template)
  }
})

tkgrid(lbl_input, sticky = "w", padx = 10, pady = 5)
tkgrid(ent_input, btn_input, padx = 10, pady = 2)
tkgrid(lbl_tissue, sticky = "w", padx = 10, pady = 5)
tkgrid(ent_tissue, padx = 10, pady = 2)
tkgrid(lbl_json, sticky = "w", padx = 10, pady = 5)
tkgrid(ent_json, btn_json, padx = 10, pady = 2)
tkgrid(lbl_out, sticky = "w", padx = 10, pady = 5)
tkgrid(ent_out, btn_out, padx = 10, pady = 2)
tkgrid(btn_convert, columnspan = 2, pady = 15)
