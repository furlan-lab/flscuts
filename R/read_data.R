#' @title Read Cell Ranger H5 File into Monocle3 cell_data_set
#'
#' @description
#' Reads a Cell Ranger output H5 file and converts it into a \code{cell_data_set} object for use with Monocle3.
#'
#' @param h5_file Path to the Cell Ranger H5 file.
#'
#' @return A \code{cell_data_set} object containing the expression data and associated metadata.
#'
#' @details
#' This function reads the H5 file produced by Cell Ranger and constructs a \code{cell_data_set} object, including expression matrices, cell metadata, and gene metadata.
#'
#' @importFrom hdf5r H5File
#' @importFrom Matrix sparseMatrix
#' @importFrom monocle3 new_cell_data_set
#' @export
#'
#' @examples
#' \dontrun{
#' cds <- read_cds_cellranger_h5_file("path/to/filtered_feature_bc_matrix.h5")
#' }
read_cds_cellranger_h5_file <- function(h5_file) {
  if (!file.exists(h5_file)) {
    stop(paste0("File ", h5_file, " not found"))
  }

  # Open the H5 file
  h5 <- hdf5r::H5File$new(h5_file, mode = "r")

  # Read data from H5 file
  barcodes <- h5[["matrix/barcodes"]][]
  gene_ids <- h5[["matrix/features/id"]][]
  gene_names <- h5[["matrix/features/name"]][]
  feature_type <- h5[["matrix/features/feature_type"]][]
  data <- h5[["matrix/data"]][]
  indices <- h5[["matrix/indices"]][] + 1  # Convert to 1-based indexing
  indptr <- h5[["matrix/indptr"]][]
  shape <- h5[["matrix/shape"]][]
  h5$close_all()

  if (length(unique(feature_type)) > 1) {
    warning("Multiple feature types found in H5 file")
  }

  # Construct the sparse matrix
  expression_matrix <- Matrix::sparseMatrix(
    i = indices,
    p = indptr,
    x = data,
    dims = shape,
    index1 = TRUE
  )

  # Create cell metadata DataFrame
  cell_metadata <- data.frame(
    barcode = barcodes,
    row.names = barcodes,
    stringsAsFactors = FALSE
  )

  # Create gene metadata DataFrame
  gene_metadata <- data.frame(
    id = gene_ids,
    gene_short_name = gene_names,
    feature_type = feature_type,
    row.names = gene_ids,
    stringsAsFactors = FALSE
  )

  # Assign row and column names to the expression matrix
  rownames(expression_matrix) <- row.names(gene_metadata)
  colnames(expression_matrix) <- row.names(cell_metadata)

  # Create the cell_data_set
  cds <- monocle3::new_cell_data_set(
    expression_data = expression_matrix,
    cell_metadata = cell_metadata,
    gene_metadata = gene_metadata
  )

  return(cds)
}


#' @title Read STARsolo Output Files into Monocle3 cell_data_set
#'
#' @description
#' Reads STARsolo output files (matrix.mtx, barcodes.tsv, features.tsv) and converts them into a \code{cell_data_set} object for use with Monocle3.
#'
#' @param folder Path to the folder containing STARsolo output files.
#'
#' @return A \code{cell_data_set} object containing the expression data and associated metadata.
#'
#' @details
#' This function reads the STARsolo output files and constructs a \code{cell_data_set} object, including expression matrices, cell metadata, and gene metadata. The function supports both compressed (.gz) and uncompressed files.
#'
#' @importFrom Matrix readMM
#' @importFrom data.table fread
#' @importFrom monocle3 new_cell_data_set
#' @export
#'
#' @examples
#' \dontrun{
#' cds <- read_cds_starsolo_file("path/to/STARsolo/output")
#' }
read_cds_starsolo_file <- function(folder) {
  if (!dir.exists(folder)) {
    stop(paste0("Folder ", folder, " not found"))
  }

  # Define possible file names
  possible_files <- list(
    files = c("features.tsv", "barcodes.tsv", "matrix.mtx"),
    genes_files = c("genes.tsv", "barcodes.tsv", "matrix.mtx"),
    files_gz = c("features.tsv.gz", "barcodes.tsv.gz", "matrix.mtx.gz"),
    genes_files_gz = c("genes.tsv.gz", "barcodes.tsv.gz", "matrix.mtx.gz")
  )

  # Check which set of files exists
  found <- sapply(possible_files, function(file_set) {
    all(file.exists(file.path(folder, file_set)))
  })

  if (!any(found)) {
    stop(paste0("Required files not found in: ", folder))
  }

  # Use the first set of files that exists
  files_to_use <- possible_files[[which(found)[1]]]
  file_paths <- file.path(folder, files_to_use)

  # Read barcodes
  barcodes <- data.table::fread(file_paths[2], header = FALSE)$V1

  # Read features (genes)
  features <- data.table::fread(file_paths[1], header = FALSE)
  gene_ids <- features$V1
  gene_names <- if (ncol(features) >= 2) features$V2 else features$V1

  # Read expression matrix
  expression_matrix <- Matrix::readMM(file_paths[3])

  # Create cell metadata DataFrame
  cell_metadata <- data.frame(
    barcode = barcodes,
    row.names = barcodes,
    stringsAsFactors = FALSE
  )

  # Create gene metadata DataFrame
  gene_metadata <- data.frame(
    id = gene_ids,
    gene_short_name = gene_names,
    row.names = gene_ids,
    stringsAsFactors = FALSE
  )

  # Assign row and column names to the expression matrix
  rownames(expression_matrix) <- row.names(gene_metadata)
  colnames(expression_matrix) <- row.names(cell_metadata)

  # Create the cell_data_set
  cds <- monocle3::new_cell_data_set(
    expression_data = expression_matrix,
    cell_metadata = cell_metadata,
    gene_metadata = gene_metadata
  )

  return(cds)
}


#' @title Find Cell Ranger Data Files
#'
#' @description
#' Identifies the correct path for Cell Ranger output files (matrix.mtx, barcodes.tsv, genes.tsv or features.tsv), possibly in compressed or uncompressed forms.
#'
#' @param folder Path to the Cell Ranger output folder.
#'
#' @return The path to the folder containing the data files. If no files are found, returns \code{NULL}.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' data_folder <- files3_prep("/path/to/cellranger/output")
#' }
files3_prep <- function(folder) {
  # Define possible filenames
  files_to_check <- list(
    uncompressed = c("matrix.mtx", "genes.tsv", "features.tsv", "barcodes.tsv"),
    compressed = c("matrix.mtx.gz", "genes.tsv.gz", "features.tsv.gz", "barcodes.tsv.gz")
  )

  # Check in the main folder and subfolders
  paths_to_check <- c(
    folder,
    file.path(folder, "outs"),
    file.path(folder, "filtered_feature_bc_matrix"),
    file.path(folder, "raw_feature_bc_matrix"),
    file.path(folder, "outs", "filtered_feature_bc_matrix"),
    file.path(folder, "outs", "raw_feature_bc_matrix")
  )

  for (path in paths_to_check) {
    if (dir.exists(path)) {
      for (file_set in files_to_check) {
        files_exist <- sapply(file_set, function(f) {
          file.exists(file.path(path, f))
        })
        if (any(files_exist)) {
          return(path)
        }
      }
    }
  }

  warning(paste0("Required files not found in: ", folder))
  return(NULL)
}


