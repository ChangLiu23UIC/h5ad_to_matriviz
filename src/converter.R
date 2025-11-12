library(tcltk)
library(Seurat)
library(dplyr)
library(arrow)
library(jsonlite)

# Function to flatten category structure for gene name matching
flatten_category_structure_gene_names <- function(category_data, available_genes = NULL) {
  flattened <- list()

  extract_categories <- function(data, current_category = "") {
    if (is.list(data)) {
      for (key in names(data)) {
        value <- data[[key]]
        # Track the current category level
        if (current_category == "") {
          # First level - primary category (e.g., "Matrisome-associated")
          extract_categories(value, key)
        } else if (current_category != "" && is.list(value)) {
          # Second level - secondary category (e.g., "ECM-affiliated Proteins")
          # This is the level we want to use as keys
          extract_categories(value, key)
        } else if (is.vector(value) && !is.list(value)) {
          # Found a gene list - this is a subcategory
          if (!is.null(available_genes)) {
            # Filter genes to only include those available in Seurat object
            # Match by exact gene name
            filtered_genes <- intersect(value, available_genes)

            if (length(filtered_genes) > 0) {  # Only add category if it has genes
              # Use the secondary category as the key
              if (!(current_category %in% names(flattened))) {
                flattened[[current_category]] <- list()
              }
              flattened[[current_category]] <- unique(c(flattened[[current_category]], filtered_genes))
            }
          } else {
            # If no available_genes provided, include all genes
            if (!(current_category %in% names(flattened))) {
              flattened[[current_category]] <- list()
            }
            flattened[[current_category]] <- unique(c(flattened[[current_category]], value))
          }
        } else {
          # Continue traversing nested structure
          extract_categories(value, current_category)
        }
      }
    } else if (is.vector(data) && !is.list(data)) {
      # Direct gene list at current category level
      if (current_category != "") {
        if (!is.null(available_genes)) {
          # Filter genes to only include those available in Seurat object
          filtered_genes <- intersect(data, available_genes)

          if (length(filtered_genes) > 0) {  # Only add category if it has genes
            if (!(current_category %in% names(flattened))) {
              flattened[[current_category]] <- list()
            }
            flattened[[current_category]] <- unique(c(flattened[[current_category]], filtered_genes))
          }
        } else {
          # If no available_genes provided, include all genes
          if (!(current_category %in% names(flattened))) {
            flattened[[current_category]] <- list()
          }
          flattened[[current_category]] <- unique(c(flattened[[current_category]], data))
        }
      }
    }
  }

  extract_categories(category_data)

  # Convert list structure to proper JSON format
  result <- list()
  for (category in names(flattened)) {
    result[[category]] <- flattened[[category]]
  }

  return(result)
}

# Function to flatten category structure for MatriViz compatibility
flatten_category_structure <- function(category_data, available_genes = NULL) {
  flattened <- list()

  # Create mapping from unversioned to versioned gene IDs
  versioned_gene_map <- list()
  if (!is.null(available_genes)) {
    for (versioned_gene in available_genes) {
      # Remove version suffix to get base ENSG ID
      base_gene <- gsub("\\.\\d+$", "", versioned_gene)
      versioned_gene_map[[base_gene]] <- versioned_gene
    }
  }

  extract_categories <- function(data, current_category = "") {
    if (is.list(data)) {
      for (key in names(data)) {
        value <- data[[key]]
        # Track the current category level
        if (current_category == "") {
          # First level - primary category (e.g., "Matrisome-associated")
          extract_categories(value, key)
        } else if (current_category != "" && is.list(value)) {
          # Second level - secondary category (e.g., "ECM-affiliated Proteins")
          # This is the level we want to use as keys
          extract_categories(value, key)
        } else if (is.vector(value) && !is.list(value)) {
          # Found a gene list - this is a subcategory
          if (!is.null(available_genes)) {
            # Filter genes to only include those available in Seurat object
            # Match by base ENSG ID (without version)
            filtered_genes <- c()
            for (base_gene in value) {
              if (base_gene %in% names(versioned_gene_map)) {
                # Use the versioned gene ID from Seurat object
                filtered_genes <- c(filtered_genes, versioned_gene_map[[base_gene]])
              }
            }

            if (length(filtered_genes) > 0) {  # Only add category if it has genes
              # Use the secondary category as the key
              if (!(current_category %in% names(flattened))) {
                flattened[[current_category]] <- list()
              }
              flattened[[current_category]] <- unique(c(flattened[[current_category]], filtered_genes))
            }
          } else {
            # If no available_genes provided, include all genes
            if (!(current_category %in% names(flattened))) {
              flattened[[current_category]] <- list()
            }
            flattened[[current_category]] <- unique(c(flattened[[current_category]], value))
          }
        } else {
          # Continue traversing nested structure
          extract_categories(value, current_category)
        }
      }
    } else if (is.vector(data) && !is.list(data)) {
      # Direct gene list at current category level
      if (current_category != "") {
        if (!is.null(available_genes)) {
          # Filter genes to only include those available in Seurat object
          filtered_genes <- c()
          for (base_gene in data) {
            if (base_gene %in% names(versioned_gene_map)) {
              # Use the versioned gene ID from Seurat object
              filtered_genes <- c(filtered_genes, versioned_gene_map[[base_gene]])
            }
          }

          if (length(filtered_genes) > 0) {  # Only add category if it has genes
            if (!(current_category %in% names(flattened))) {
              flattened[[current_category]] <- list()
            }
            flattened[[current_category]] <- unique(c(flattened[[current_category]], filtered_genes))
          }
        } else {
          # If no available_genes provided, include all genes
          if (!(current_category %in% names(flattened))) {
            flattened[[current_category]] <- list()
          }
          flattened[[current_category]] <- unique(c(flattened[[current_category]], data))
        }
      }
    }
  }

  extract_categories(category_data)

  # Convert list structure to proper JSON format
  result <- list()
  for (category in names(flattened)) {
    result[[category]] <- flattened[[category]]
  }

  return(result)
}

convert_azimuth <- function(rds_path, tissue, output_dir, json_template = NULL) {
  cat("Loading reference...\n")
  ref <- readRDS(rds_path)
  
  # --- Ensure compatibility with new Seurat structure ---
  ref <- UpdateSeuratObject(ref)
  
  tissue_lower <- tolower(tissue)
  
  assays_available <- Assays(ref)
  if (!("integrated" %in% assays_available)) stop("integrated assay not found in Seurat object.")
  
  cat("Extracting integrated assay data...\n")
  expr_sparse <- GetAssayData(ref, assay = "integrated", layer = "data")

  # --- Select variable features to avoid huge memory load ---
  if (length(VariableFeatures(ref)) > 0) {
    genes_use <- VariableFeatures(ref)
    cat("Using variable features (", length(genes_use), " genes)...\n")
  } else {
    genes_use <- rownames(expr_sparse)[1:min(2000, nrow(expr_sparse))]
    cat("No variable features found — using first", length(genes_use), "genes.\n")
  }

  expr_sparse <- expr_sparse[genes_use, , drop = FALSE]

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
  meta_cols <- intersect(colnames(ref@meta.data),
                         c("annotation.l1", "annotation.l2", "annotation.l3",
                           "predicted.celltype.l1", "predicted.celltype.l2",
                           "celltype.l1", "celltype.l2"))
  
  meta <- ref@meta.data
  meta$index <- rownames(meta)
  if (length(meta_cols) > 0) {
    meta <- meta[, c("index", meta_cols), drop = FALSE]
  }
  
  # --- Merge data ---
  cat("Merging data...\n")
  # Transpose expression data to get cells as rows for merging with UMAP
  expr_t <- as.data.frame(t(expr_df[, !colnames(expr_df) %in% "gene"]))
  expr_t$index <- rownames(expr_t)

  merged <- left_join(umap, expr_t, by = "index")
  if (length(meta_cols) > 0) merged <- left_join(merged, meta, by = "index")
  
  # --- Add mean gene expression summary ---
  gene_cols <- setdiff(colnames(merged), c("index", "UMAP_1", "UMAP_2", meta_cols))
  merged$All_Genes <- rowMeans(merged[, gene_cols], na.rm = TRUE)
  merged <- merged[, !grepl("celltype", colnames(merged), ignore.case = TRUE)]
  
  # --- File naming ---
  expr_file <- file.path(output_dir, paste0(tissue_lower, "_rna_v1.parquet"))
  centroid_file <- file.path(output_dir, paste0(tissue_lower, "_rna_centroid_v1.0.parquet"))
  category_file <- file.path(output_dir, paste0(tissue_lower, "_rna_v1_category.json"))
  meta_file <- file.path(output_dir, paste0(tissue_lower, "_rna.json"))
  
  # --- Save expression parquet ---
  cat("Writing expression parquet...\n")
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
    # Compute centroids for annotation labels (these will be used for category positioning)
    centroids <- ref@meta.data %>%
      dplyr::select(all_of(label_col)) %>%
      mutate(index = rownames(ref@meta.data)) %>%
      left_join(umap, by = "index") %>%
      group_by(.data[[label_col]]) %>%
      summarise(
        cen_x = mean(UMAP_1, na.rm = TRUE),
        cen_y = mean(UMAP_2, na.rm = TRUE)
      ) %>%
      arrange(.data[[label_col]]) %>%
      rename(Type = .data[[label_col]])
  } else {
    # Create empty centroids dataframe if no cell type columns found
    centroids <- data.frame(
      Type = character(),
      cen_x = numeric(),
      cen_y = numeric()
    )
    cat("No cell type columns found - creating empty centroids file\n")
  }
  
  write_parquet(centroids, centroid_file, compression = "zstd", version = "1.0")
  
  # --- Category JSON processing ---
  cat("Processing category JSON...\n")

  # Get available genes from integrated assay (genes are rownames in sparse matrix)
  available_genes <- rownames(expr_sparse)

  # Initialize category JSON structure
  category_json <- list()

  # Process matrisome gene categories if template provided
  if (!is.null(json_template) && file.exists(json_template)) {
    # Load JSON template
    template_json <- fromJSON(json_template, simplifyVector = FALSE)

    # Flatten matrisome category structure for gene name matching
    matrisome_categories <- flatten_category_structure_gene_names(template_json, available_genes)

    # Add matrisome categories to the main category JSON (genes only, no centroids)
    for (category_name in names(matrisome_categories)) {
      if (length(matrisome_categories[[category_name]]) > 0) {
        category_json[[category_name]] <- matrisome_categories[[category_name]]
      }
    }
  }

  # If no categories were created, use default structure
  if (length(category_json) == 0) {
    category_json <- list(
      "ECM Glycoproteins" = list(),
      "Secreted Factors" = list(),
      "Collagens" = list(),
      "ECM Regulators" = list(),
      "Proteoglycans" = list(),
      "ECM-affiliated Proteins" = list()
    )
  }

  write_json(category_json, category_file, pretty = TRUE, auto_unbox = TRUE)
  
  # --- Metadata JSON ---
  meta_json <- list(
    fileType = "matriviz",
    version = "1.1.0",
    category_name = tissue_lower,
    category_description = paste(tissue, "(RNA reference)"),
    parquet_file = basename(expr_file),
    category_file = basename(category_file),
    centroid_file = basename(centroid_file)
  )
  write_json(meta_json, meta_file, pretty = TRUE, auto_unbox = TRUE)
  
  tkmessageBox(title = "Success",
               message = paste("✅ RNA Conversion complete!\n\nSaved to:", output_dir))
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
  rds_path <- tclvalue(tkget(ent_input))
  tissue <- tclvalue(tkget(ent_tissue))
  json_template <- tclvalue(tkget(ent_json))
  out_dir <- tclvalue(tkget(ent_out))
  if (rds_path == "" || tissue == "" || out_dir == "") {
    tkmessageBox(title = "Error", message = "Please fill in all required fields.")
  } else {
    # Convert empty string to NULL for optional parameter
    if (json_template == "") {
      json_template <- NULL
    }
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
