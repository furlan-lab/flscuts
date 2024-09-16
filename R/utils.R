#' @title Convert a GMT File to a List of Gene Sets
#'
#' @description
#' This function reads a Gene Matrix Transposed (GMT) file and converts it into a list where each element corresponds to a gene set. The names of the list elements are the gene set names, and the contents are vectors of genes belonging to each set.
#'
#' @param gmt_file A character string specifying the path to the GMT file.
#' @return A named list where each element is a vector of gene symbols for a gene set.
#' @details
#' The GMT file format is commonly used for storing gene sets, such as those provided by MSigDB. Each line in a GMT file corresponds to a gene set and is formatted as follows:
#' \enumerate{
#'   \item Gene set name
#'   \item Description (can be empty)
#'   \item Gene symbols (tab-separated)
#' }
#' This function reads the GMT file line by line, splits each line by tabs, and extracts the gene set name and the list of genes.
#'
#' @importFrom utils read.table
#' @export
#'
#' @examples
#' \dontrun{
#' # Convert a GMT file to a list of gene sets
#' gene_sets <- GMT_to_list("path/to/your/file.gmt")
#' }
GMT_to_list <- function(gmt_file) {
  # Check if the file exists
  if (!file.exists(gmt_file)) {
    stop("The specified GMT file does not exist.")
  }

  # Read the GMT file
  gmt_lines <- readLines(gmt_file)

  # Initialize an empty list to store gene sets
  gene_sets <- list()

  # Process each line
  for (line in gmt_lines) {
    # Split the line by tabs
    elements <- unlist(strsplit(line, "\t"))

    # Ensure that the line has at least three elements
    if (length(elements) < 3) {
      warning("A line in the GMT file does not have the expected format and will be skipped.")
      next
    }

    # Extract gene set name and genes
    gene_set_name <- elements[1]
    genes <- elements[3:length(elements)]

    # Add to the list
    gene_sets[[gene_set_name]] <- genes
  }

  return(gene_sets)
}


#' @title Convert Monocle3 Cell Data Set to Seurat Object
#'
#' @description
#' Converts a Monocle3 `cell_data_set` object into a Seurat object. Optionally, dimensionality reduction embeddings can be transferred.
#'
#' @param cds A Monocle3 `cell_data_set` object to convert.
#' @param assay_name A character string specifying the assay name for the Seurat object. Default is \code{"RNA"}.
#' @param monocle_reduction The name of the reduction in Monocle3 to transfer (e.g., \code{"UMAP"}). Use \code{NULL} to skip transferring reductions. Default is \code{"UMAP"}.
#' @param seurat_reduction The name to assign to the reduction in the Seurat object (e.g., \code{"umap"}). Use \code{NULL} to skip transferring reductions. Default is \code{"umap"}.
#' @param row.names A character string specifying the column in \code{rowData(cds)} to use as feature names. Default is \code{"gene_short_name"}.
#' @param handle_duplicate_features How to handle duplicate feature names. Options are \code{"make_unique"} or \code{"remove"}. Default is \code{"make_unique"}.
#' @param normalize Logical, whether to run \code{NormalizeData()} on the Seurat object. Default is \code{TRUE}.
#'
#' @return A Seurat object.
#'
#' @details
#' This function converts a Monocle3 `cell_data_set` object into a Seurat object by extracting the counts data, cell metadata, and feature metadata. It can also transfer dimensionality reduction embeddings from the Monocle3 object to the Seurat object if specified.
#'
#' @importFrom Seurat CreateSeuratObject CreateDimReducObject NormalizeData
#' @importFrom SingleCellExperiment reducedDims
#' @importFrom SummarizedExperiment assay
#' @export
#'
#' @examples
#' \dontrun{
#' # Assuming 'cds' is a Monocle3 cell_data_set
#' seurat_obj <- monocle3_to_seurat(
#'   cds,
#'   assay_name = "RNA",
#'   monocle_reduction = "UMAP",
#'   seurat_reduction = "umap"
#' )
#' }
monocle3_to_seurat <- function(cds,
                               assay_name = "RNA",
                               monocle_reduction = "UMAP",
                               seurat_reduction = "umap",
                               row.names = "gene_short_name",
                               handle_duplicate_features = c("make_unique", "remove"),
                               normalize = TRUE) {
  # Input validation
  if (!inherits(cds, "cell_data_set")) {
    stop("Input must be a Monocle3 cell_data_set object.")
  }

  handle_duplicate_features <- match.arg(handle_duplicate_features)

  # Extract counts data
  counts <- SummarizedExperiment::assay(cds)

  # Extract feature names
  feature_names <- SummarizedExperiment::rowData(cds)[[row.names]]
  if (is.null(feature_names)) {
    stop(paste("Row data column", row.names, "not found in the cell_data_set object."))
  }

  # Set feature names
  rownames(counts) <- feature_names

  # Handle duplicate feature names
  if (handle_duplicate_features == "remove") {
    duplicate_features <- duplicated(rownames(counts))
    if (any(duplicate_features)) {
      counts <- counts[!duplicate_features, ]
    }
  } else if (handle_duplicate_features == "make_unique") {
    rownames(counts) <- make.unique(rownames(counts))
  }

  # Extract cell metadata
  cell_metadata <- as.data.frame(SummarizedExperiment::colData(cds))

  # Create the Seurat object
  seurat_obj <- Seurat::CreateSeuratObject(
    counts = counts,
    assay = assay_name,
    meta.data = cell_metadata
  )

  # Transfer reduction if specified
  if (!is.null(monocle_reduction) && !is.null(seurat_reduction)) {
    if (monocle_reduction %in% names(SingleCellExperiment::reducedDims(cds))) {
      embeddings <- SingleCellExperiment::reducedDims(cds)[[monocle_reduction]]
      colnames(embeddings) <- paste0(seurat_reduction, "_", 1:ncol(embeddings))
      seurat_obj[[seurat_reduction]] <- Seurat::CreateDimReducObject(
        embeddings = embeddings,
        key = paste0(seurat_reduction, "_"),
        assay = assay_name
      )
    } else {
      warning(paste("Reduction", monocle_reduction, "not found in the cell_data_set object."))
    }
  }

  # Normalize data if requested
  if (normalize) {
    seurat_obj <- Seurat::NormalizeData(seurat_obj)
  }

  return(seurat_obj)
}


#' @title Convert Seurat Object to Monocle3 Cell Data Set
#'
#' @description
#' Converts a Seurat object into a Monocle3 `cell_data_set` object. Optionally, dimensionality reduction embeddings can be transferred.
#'
#' @param seurat_obj A Seurat object to convert.
#' @param assay_name A character string specifying the assay to use from the Seurat object. Default is \code{"RNA"}.
#' @param seurat_reduction The name of the reduction in Seurat to transfer (e.g., \code{"umap"}). Use \code{NULL} to skip transferring reductions. Default is \code{"umap"}.
#' @param monocle_reduction The name to assign to the reduction in the Monocle3 object (e.g., \code{"UMAP"}). Use \code{NULL} to skip transferring reductions. Default is \code{"UMAP"}.
#'
#' @return A Monocle3 `cell_data_set` object.
#'
#' @details
#' This function converts a Seurat object into a Monocle3 `cell_data_set` object by extracting the counts data, cell metadata, and gene metadata. It can also transfer dimensionality reduction embeddings from the Seurat object to the Monocle3 object if specified.
#'
#' @importFrom Seurat GetAssayData Embeddings
#' @importFrom monocle3 new_cell_data_set
#' @importFrom SingleCellExperiment reducedDims<-
#' @export
#'
#' @examples
#' \dontrun{
#' # Assuming 'seurat_obj' is a Seurat object
#' cds <- seurat_to_monocle3(
#'   seurat_obj,
#'   assay_name = "RNA",
#'   seurat_reduction = "umap",
#'   monocle_reduction = "UMAP"
#' )
#' }
seurat_to_monocle3 <- function(seurat_obj,
                               assay_name = "RNA",
                               seurat_reduction = "umap",
                               monocle_reduction = "UMAP") {
  # Input validation
  if (!inherits(seurat_obj, "Seurat")) {
    stop("Input must be a Seurat object.")
  }
  if (!assay_name %in% names(seurat_obj@assays)) {
    stop(paste("Assay", assay_name, "not found in the Seurat object."))
  }

  # Extract counts data
  counts <- Seurat::GetAssayData(seurat_obj, assay = assay_name, slot = "counts")

  # Extract cell metadata
  cell_metadata <- seurat_obj@meta.data

  # Create gene metadata
  gene_metadata <- data.frame(
    gene_id = rownames(counts),
    gene_short_name = rownames(counts),
    row.names = rownames(counts),
    stringsAsFactors = FALSE
  )

  # Create the cell_data_set object
  cds <- monocle3::new_cell_data_set(
    counts,
    cell_metadata = cell_metadata,
    gene_metadata = gene_metadata
  )

  # Transfer reduction if specified
  if (!is.null(seurat_reduction) && !is.null(monocle_reduction)) {
    if (seurat_reduction %in% names(seurat_obj@reductions)) {
      embeddings <- Seurat::Embeddings(seurat_obj, reduction = seurat_reduction)
      SingleCellExperiment::reducedDims(cds)[[monocle_reduction]] <- embeddings
    } else {
      warning(paste("Reduction", seurat_reduction, "not found in the Seurat object."))
    }
  }

  return(cds)
}

#' @export
#' @title Add Souporcell Clustering Output to a Seurat Object
#' @description This function adds an assay to a Seurat object based on Souporcell clustering results.
#' The assay contains log-transformed and normalized cluster probabilities from the Souporcell `clusters.tsv` file.
#' Principal Component Analysis (PCA) is performed on these probabilities, and the resulting components are added as a dimensionality reduction object.
#' The Souporcell 'assignment' is added to the Seurat object's metadata under the specified label.
#' Optionally, assignments can be renamed to be 1-indexed and multiplets collapsed into a single category.
#'
#' @param seurat_obj A Seurat object.
#' @param souporcell_file Path to the Souporcell `clusters.tsv` file.
#' @param prefix Optional prefix to prepend to cell barcodes.
#' @param assay_name Name of the assay to add to the Seurat object. Default is "GENO".
#' @param key Key for the dimensionality reduction. Default is "gpca_".
#' @param meta_data_label Name of the metadata column to store Souporcell assignments. Default is "geno".
#' @param rd_label Name of the dimensionality reduction object to store PCA results. Default is "gpca".
#' @param rename_assignments Logical indicating whether to rename Souporcell assignments to be 1-indexed and collapse multiplets. Default is TRUE.
#' @return A Seurat object with the added assay, dimensionality reduction, and updated metadata.
#' @importFrom Seurat CreateAssayObject CreateDimReducObject
#' @examples
#' # Assuming 'seurat_obj' is your Seurat object and 'clusters.tsv' is your Souporcell output file:
#' seurat_obj <- add_souporcell_seurat(seurat_obj, "path/to/clusters.tsv")
add_souporcell_seurat <- function(
    seurat_obj,
    souporcell_file,
    prefix = NULL,
    assay_name = "GENO",
    key = "gpca_",
    meta_data_label = "geno",
    rd_label = "gpca",
    rename_assignments = TRUE
) {
  # Check if the file exists
  if (length(souporcell_file) > 1) {
    stop("Only supports one file addition at a time.")
  }
  if (!file.exists(souporcell_file)) {
    stop("Souporcell file does not exist!")
  }

  # Check that seurat_obj is a Seurat object
  if (!inherits(seurat_obj, "Seurat")) {
    stop("Input object is not a Seurat object.")
  }

  # Read Souporcell data
  souporcell_data <- read.table(souporcell_file, header = TRUE, stringsAsFactors = FALSE)

  # Check for required columns
  required_cols <- c("barcode", "status", "assignment", "log_prob_singleton", "log_prob_doublet")
  if (!all(required_cols %in% colnames(souporcell_data))) {
    stop(sprintf(
      "Souporcell file %s does not have the required columns: %s",
      souporcell_file, paste(required_cols, collapse = ", ")
    ))
  }

  # Set row names
  if (is.null(prefix)) {
    rownames(souporcell_data) <- souporcell_data$barcode
  } else {
    rownames(souporcell_data) <- paste0(prefix[1], souporcell_data$barcode)
  }

  # Check that all cells in seurat_obj are in souporcell_data
  missing_cells <- setdiff(Cells(seurat_obj), rownames(souporcell_data))
  if (length(missing_cells) > 0) {
    stop("Not all cells in the Seurat object are found in the Souporcell data.")
  }

  # Subset souporcell_data to only include cells in seurat_obj
  souporcell_data <- souporcell_data[Cells(seurat_obj), ]

  # Extract cluster probabilities
  cluster_cols <- grep("^cluster", colnames(souporcell_data), value = TRUE)
  if (length(cluster_cols) == 0) {
    stop("No cluster columns found in Souporcell data.")
  }
  cluster_probs <- as.matrix(souporcell_data[, cluster_cols])

  # Normalize and log-transform cluster probabilities
  col_means <- colMeans(cluster_probs)
  cluster_probs_normalized <- sweep(cluster_probs, 2, col_means, "/")

  # Handle zeros and negative values before log transformation
  epsilon <- .Machine$double.eps
  cluster_probs_normalized[cluster_probs_normalized <= 0] <- epsilon
  cluster_probs_log <- log10(cluster_probs_normalized)

  # Perform PCA
  pca_results <- stats::prcomp(cluster_probs_log, center = TRUE, scale. = TRUE)

  # Prepare PCA scores
  pca_scores <- pca_results$x
  rownames(pca_scores) <- souporcell_data$barcode
  colnames(pca_scores) <- paste0("PC", 1:ncol(pca_scores))

  # Prepare assignments
  assignments <- souporcell_data$assignment
  if (rename_assignments) {
    fixed_assignments <- fix_assignment(assignments)
  } else {
    fixed_assignments <- assignments
  }

  # Create assay
  assay_data <- t(cluster_probs_log)
  seurat_obj[[assay_name]] <- Seurat::CreateAssayObject(data = assay_data)
  DefaultAssay(seurat_obj) <- assay_name

  # Create dimensionality reduction object
  seurat_obj[[rd_label]] <- Seurat::CreateDimReducObject(
    embeddings = pca_scores,
    key = key,
    assay = assay_name
  )

  # Update metadata
  if (!(meta_data_label %in% colnames(seurat_obj@meta.data))) {
    seurat_obj@meta.data[[meta_data_label]] <- NA
  }
  seurat_obj@meta.data[Cells(seurat_obj), meta_data_label] <- fixed_assignments

  # Return updated Seurat object
  return(seurat_obj)
}

#' Helper function for add_souporcell_seurat
fix_assignment <- function(vector) {
  vector[grepl("\\/", vector)] <- "Multiplet"
  vector[vector != "Multiplet"] <- as.numeric(vector[vector != "Multiplet"]) + 1
  vector
}

#'
#' #' @title Lighten or Darken a Color
#' #' @description Adjusts the brightness of a hex color by a specified amount.
#' #' @param col A character string representing a hex color (e.g., "#FF0000").
#' #' @param amt Integer value between -255 and 255 indicating the amount to lighten or darken the color.
#' #' @return A hex color string representing the adjusted color.
#' #' @export
#' lighten_darken_color <- function(col, amt) {
#'   if (substring(col, 1, 1) == "#") {
#'     col <- substring(col, 2)
#'   }
#'   num <- as.integer(paste0("0x", col))
#'   r <- (num >> 16) & 0xFF
#'   g <- (num >> 8) & 0xFF
#'   b <- num & 0xFF
#'
#'   adjust <- function(x) {
#'     x <- x + amt
#'     x <- max(min(255, x), 0)
#'     return(x)
#'   }
#'
#'   r <- adjust(r)
#'   g <- adjust(g)
#'   b <- adjust(b)
#'
#'   new_col <- sprintf("#%02X%02X%02X", r, g, b)
#'   return(new_col)
#' }


#' @title Generate a Color Palette
#' @description Creates a color palette of specified length, optionally scrambled.
#' @param n Integer specifying the number of colors to generate.
#' @param scramble Logical indicating whether to randomly shuffle the colors. Default is FALSE.
#' @return A character vector of hex color codes.
#' @export
sfc <- function(n, scramble = FALSE) {
  if (!is.numeric(n) || n <= 0 || n %% 1 != 0) {
    stop("Please input a positive integer for 'n'.")
  }
  base_colors <- c(
    "#16482A", "#1C7C40", "#45AC49", "#69BC9A", "#FBD43F",
    "#E77A2B", "#DC3F32", "#932528", "#50191E", "#96C4D9",
    "#2394C4", "#4575AD", "#8681B0", "#6C5492", "#8C4A8D",
    "#9E2563", "#492C74", "#E9E52F", "#F8C566", "#D85191"
  )
  palette_func <- colorRampPalette(base_colors)
  colors <- palette_func(n)
  if (scramble) {
    colors <- sample(colors)
  }
  return(colors)
}



