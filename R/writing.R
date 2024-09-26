
#' @title Write a compressed MatrixMarket file
#'
#' @description This function writes a sparse matrix in the MatrixMarket format to a compressed `.gz` file.
#' The function handles both real and integer matrix types.
#'
#' @param x A sparse matrix (typically a \code{dgCMatrix} or \code{ngCMatrix} object).
#' @param file A character string specifying the output file name, which will be compressed into `.gz` format.
#'
#' @details
#' This function writes the matrix in the MatrixMarket coordinate format.
#' It first writes the header indicating the matrix type and size, and then appends the matrix data.
#' If the matrix is an `ngCMatrix`, it is treated as an integer matrix, otherwise as a real matrix.
#' The function compresses the output into a `.gz` file.
#'
#' @importFrom data.table fwrite
#' @importFrom Matrix summary
#' @keywords internal
#' @return This function does not return a value. It writes a file as a side effect.
#' @export
#' @examples
#' \dontrun{
#' library(Matrix)
#' m <- Matrix(c(0, 1, 0, 2), 2, 2, sparse = TRUE)
#' writeMMgz(m, "matrix.mtx.gz")
#' }
writeMMgz <- function(x, file) {
  mtype <- "real"
  if (is(x, "ngCMatrix")) {
    mtype <- "integer"
  }
  writeLines(
    c(
      sprintf("%%%%MatrixMarket matrix coordinate %s general", mtype),
      sprintf("%s %s %s", x@Dim[1], x@Dim[2], length(x@x))
    ),
    gzfile(file)
  )
  fwrite(
    x = summary(x),
    file = file,
    append = TRUE,
    sep = " ",
    row.names = FALSE,
    col.names = FALSE
  )
}




#' @title Export Seurat Object Data to 10X-Style Format with Optional Reductions
#'
#' @description
#' This function exports data from a Seurat object into a 10X Genomics-style format. The output includes files for the expression matrix, feature (gene) information, barcodes, metadata, dimensional reductions (e.g., UMAP, PCA), and variable features. These files are written in a compressed format where applicable.
#'
#' @param seurat_obj A Seurat object containing the data to be exported.
#' @param assay A character string indicating which assay to use from the Seurat object. Default is \code{"RNA"}.
#' @param dir A character string specifying the directory where the output files will be saved. The directory must already exist.
#' @param include_reductions Logical, whether to include cell embeddings from reductions (e.g., UMAP, PCA) in the output. Default is \code{TRUE}.
#'
#' @details
#' The function creates several files in a subdirectory called \code{"export"} within the specified directory:
#' \itemize{
#'   \item \code{matrix.mtx.gz}: A compressed Matrix Market file containing the assay data (expression matrix).
#'   \item \code{features.tsv.gz}: A tab-separated file with feature (gene) information, including gene IDs and names.
#'   \item \code{barcodes.tsv.gz}: A tab-separated file with cell barcodes.
#'   \item \code{metadata.csv}: A CSV file containing metadata from the Seurat object.
#'   \item \code{<reduction>_embeddings.tsv.gz}: A compressed file with cell embeddings for each reduction (e.g., UMAP, PCA), if \code{include_reductions} is set to \code{TRUE}.
#'   \item \code{variable_features.tsv.gz}: A compressed file listing the variable features.
#' }
#'
#' If reductions (like UMAP or PCA) are present in the Seurat object and \code{include_reductions} is \code{TRUE}, the cell embeddings from each reduction will be written to separate files in the format \code{<reduction>_embeddings.tsv.gz}.
#' If no reductions are found and \code{include_reductions} is \code{TRUE}, the function will issue a warning but will still generate the other files.
#' The function will create the \code{"export"} subdirectory within the specified directory if it doesn't exist.
#'
#' @importFrom Seurat GetAssayData Cells VariableFeatures Assays
#' @importFrom utils write.csv write.table
#' @importFrom Matrix writeMM
#' @importFrom methods is
#' @export
#' @return The function does not return a value. It writes several files as a side effect.
#'
#' @examples
#' \dontrun{
#' library(Seurat)
#' seurat_obj <- CreateSeuratObject(counts = matrix(rnorm(100), 10, 10))
#' export_seurat_to_10x(seurat_obj, assay = "RNA", dir = "output_directory")
#'
#' # Export Seurat object with reductions
#' export_seurat_to_10x(seurat_obj, assay = "RNA", dir = "output_directory", include_reductions = TRUE)
#'
#' # Export Seurat object without reductions
#' export_seurat_to_10x(seurat_obj, assay = "RNA", dir = "output_directory", include_reductions = FALSE)
#' }
export_seurat_to_10x <- function(seurat_obj, assay = "RNA", dir, include_reductions = TRUE) {
  # Check if the directory exists
  if (!dir.exists(dir)) {
    stop("The specified directory does not exist. Please provide a valid directory.")
  }

  # Check if the assay exists in the Seurat object
  if (!assay %in% Seurat::Assays(seurat_obj)) {
    stop(paste("Assay", assay, "not found in the Seurat object."))
  }

  # Create the export subdirectory if it doesn't exist
  out_dir <- file.path(dir, "export")
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }

  # Get assay data
  assay_data <- Seurat::GetAssayData(seurat_obj, assay = assay, slot = "counts")

  # Write matrix.mtx.gz file (expression matrix)
  message("Writing expression matrix...")
  matrix_file <- file.path(out_dir, "matrix.mtx.gz")
  gz1 <- gzfile(matrix_file, "w")
  Matrix::writeMM(assay_data, file = gz1)
  close(gz1)

  # Write features.tsv.gz file (feature/gene information)
  message("Writing feature information...")
  features <- data.frame(
    gene_id = rownames(assay_data),
    gene_name = rownames(assay_data),
    feature_type = rep("Gene Expression", nrow(assay_data)),
    stringsAsFactors = FALSE
  )
  features_file <- file.path(out_dir, "features.tsv.gz")
  gz2 <- gzfile(features_file, "w")
  write.table(features, gz2, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  close(gz2)

  # Write barcodes.tsv.gz file (cell barcodes)
  message("Writing cell barcodes...")
  barcodes <- Seurat::Cells(seurat_obj)
  barcodes_file <- file.path(out_dir, "barcodes.tsv.gz")
  gz3 <- gzfile(barcodes_file, "w")
  write.table(barcodes, gz3, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  close(gz3)

  # Write metadata.csv
  message("Writing metadata...")
  metadata_file <- file.path(out_dir, "metadata.csv")
  write.csv(seurat_obj@meta.data, file = metadata_file, quote = TRUE, row.names = TRUE)

  # Write variable features to variable_features.tsv.gz
  message("Writing variable features...")
  variable_features <- Seurat::VariableFeatures(seurat_obj)
  if (length(variable_features) > 0) {
    var_features_file <- file.path(out_dir, "variable_features.tsv.gz")
    gz4 <- gzfile(var_features_file, "w")
    write.table(variable_features, gz4, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
    close(gz4)
  } else {
    warning("No variable features found in the Seurat object.")
  }

  # Write reductions if requested
  if (include_reductions) {
    reductions <- names(seurat_obj@reductions)
    if (length(reductions) == 0) {
      warning("No reductions found in the Seurat object.")
    } else {
      for (red in reductions) {
        message(paste("Writing reduction:", red))
        # Retrieve the reduction embeddings
        embeddings <- seurat_obj@reductions[[red]]@cell.embeddings
        # Open gzipped reduction file
        red_file <- file.path(out_dir, paste0(red, "_embeddings.tsv.gz"))
        gz_red <- gzfile(red_file, "w")
        # Write the reduction embeddings to the file
        write.table(embeddings, gz_red, sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)
        close(gz_red)
      }
    }
  }

  message("Export completed.")
}


