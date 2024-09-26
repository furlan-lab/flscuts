#' Project Data into UMAP Embedding of Single-Cell Object (Compatible with Monocle3 and Seurat)
#'
#' This function projects data from one single-cell object (or a similar structure) into the UMAP embedding of another object. The embedding must be precomputed using methods such as LSI or PCA. This allows for co-embedding of new data into an existing embedding space.
#'
#' @details The function supports projections using models built from UMAP embeddings, which are ideal for comparing single-cell data, and can also simulate pseudo-single cells for bulk data. It is adapted from: *ArchR: An integrative and scalable software package for single-cell chromatin accessibility analysis* by Granja et al. (2020).
#'
#' @param projector A `cell_data_set` or Seurat object containing a reduced dimension matrix (e.g., LSI) and the model used to generate a UMAP embedding.
#' @param projectee A `SummarizedExperiment`, `cell_data_set`, or Seurat object that will be projected into the UMAP space defined by the projector.
#' @param make_pseudo_single_cells Logical, whether to simulate pseudo-single cells from the `projectee` data, useful for bulk data. Default is `FALSE`.
#' @param ncells_coembedding Numeric, the number of cells to use from the `projector` for co-embedding with the pseudo-single cells. Default is 5000, or the total number of cells if less than 5000.
#' @param reduced_dim Character, specifies the reduced dimension method (e.g., "LSI" or "PCA"). Default is "LSI".
#' @param embedding Character, specifies the embedding type to be used (currently only "UMAP" supported). Default is "UMAP".
#' @param n Integer, the number of subsampled "pseudo single cells" per bulk sample, relevant when `make_pseudo_single_cells` is `TRUE`. Default is 250.
#' @param verbose Logical, whether to print verbose output during execution. Default is `TRUE`.
#' @param threads Integer, the number of threads for parallel execution. Default is 6.
#' @param seed Integer, seed for reproducibility. Default is 2020.
#' @param force Logical, whether to continue even if the overlap ratio of features between the `projector` and `projectee` is below a certain threshold. Default is `FALSE`.
#' @param scale Logical, whether to scale the projected data after projection. Default is `FALSE`.
#' @param features Character vector, specifying whether the projection should be based on "annotation-based" or "range-based" features. Default is `c("annotation-based", "range-based")`.
#'
#' @return A `SimpleList` object containing the UMAP coordinates of the projected data (both original single-cell and simulated bulk data) and the reduced dimension matrix.
#' @export
project_data <- function(
    projector = NULL,
    projectee = NULL,
    ncells_coembedding = 5000,
    scale = FALSE,
    reduced_dim = "LSI",
    embedding = "UMAP",
    make_pseudo_single_cells = FALSE,
    features = c("annotation-based", "range-based"),
    n = 250,
    verbose = TRUE,
    threads = 6,
    seed = 2020,
    force = FALSE
){
  features <- match.arg(features)
  # Check object types
  if (methods::is(projector, "Seurat")) {
    object_type_projector <- "seurat"
  } else if (methods::is(projector, "cell_data_set")) {
    object_type_projector <- "monocle3"
  } else {
    stop("The projector must be a Seurat object or a Monocle3 cell_data_set object.")
  }

  if (methods::is(projectee, "Seurat") || methods::is(projectee, "cell_data_set") || methods::is(projectee, "SummarizedExperiment")) {
    object_type_projectee <- "seurat_or_sce"
  } else {
    stop("The projectee must be a Seurat object, a Monocle3 cell_data_set object, or a SummarizedExperiment object.")
  }

  # Extract reduced dimensions and embeddings from the projector
  if (object_type_projector == "monocle3") {
    if (!reduced_dim %in% names(SingleCellExperiment::reducedDims(projector))) {
      stop(paste("Reduced dimension", reduced_dim, "not found in projector"))
    }
    rD <- SingleCellExperiment::reducedDims(projector)[[reduced_dim]]
    if (!embedding %in% names(SingleCellExperiment::reducedDims(projector))) {
      stop(paste("Embedding", embedding, "not found in projector"))
    }
    sc_embedding <- SingleCellExperiment::reducedDims(projector)[[embedding]]
    rownames(sc_embedding) <- colnames(projector)
    # Retrieve features and LSI model
    if (is.null(projector@int_metadata$LSI_model)) {
      stop("LSI model not found in projector object.")
    }
    lsi_model <- projector@int_metadata$LSI_model
  } else if (object_type_projector == "seurat") {
    if (!reduced_dim %in% names(projector@reductions)) {
      stop(paste("Reduced dimension", reduced_dim, "not found in projector"))
    }
    rD <- projector@reductions[[reduced_dim]]@cell.embeddings
    if (!embedding %in% names(projector@reductions)) {
      stop(paste("Embedding", embedding, "not found in projector"))
    }
    sc_embedding <- projector@reductions[[embedding]]@cell.embeddings
    rownames(sc_embedding) <- colnames(projector)
    # Retrieve features and LSI model
    if (is.null(projector@reductions[[reduced_dim]]@misc$svd)) {
      stop("LSI model not found in projector object.")
    }
    lsi_model <- projector@reductions[[reduced_dim]]@misc
  }

  embedding_num_dim <- lsi_model$num_dim

  # Extract features used in the reduced dimension
  if (features == "annotation-based") {
      features_used <- lsi_model$features
    } else if (features == "range-based") {
      if (object_type_projector == "monocle3") {
        features_used <- lsi_model$granges
      } else if (object_type_projector == "seurat") {
        stop("Range-based features not supported for Seurat objects.")
      }
  }

  # Extract shared data between projector and projectee
  shared_rd <- extract_data(features_used, projectee)
  # Calculate overlap ratio
  overlap_ratio <- shared_rd$overlap
  message(paste0("Overlap Ratio of Reduced Dims Features = ", round(overlap_ratio, 3)))

  if (overlap_ratio < 0.25) {
    if (force) {
      warning("Less than 25% of the features are present in the projectee data set! Continuing since force = TRUE!")
    } else {
      stop("Less than 25% of the features are present in the projectee data set! Set force = TRUE to continue!")
    }
  }

  # Simulate single cells and project using original LSI/SVD model
  if (make_pseudo_single_cells) {
    depthN <- round(sum(lsi_model$row_sums) / nrow(rD))
    nRep <- 5
    n2 <- ceiling(n / nRep)
    ratios <- c(2, 1.5, 1, 0.5, 0.25) # Range of ratios of number of fragments

    if (verbose) message(paste0("Simulating ", (n * dim(shared_rd$mat)[2]), " single cells"))
    projRD <- pbmcapply::pbmclapply(seq_len(ncol(shared_rd$mat)), function(x){
      counts <- shared_rd$mat[, x]
      counts <- rep(seq_along(counts), counts)
      simMat <- lapply(seq_len(nRep), function(y){
        ratio <- ratios[y]
        simMat <- matrix(sample(x = counts, size = ceiling(ratio * depthN) * n2, replace = TRUE), ncol = n2)
        simMat <- Matrix::summary(as(simMat, "dgCMatrix"))[, -1, drop = FALSE]
        simMat[, 1] <- simMat[, 1] + (y - 1) * n2
        simMat
      }) %>% Reduce("rbind", .)
      simMat <- Matrix::sparseMatrix(i = simMat[, 2], j = simMat[, 1], x = rep(1, nrow(simMat)), dims = c(nrow(shared_rd$mat), n2 * nRep))
      projRD <- as.matrix(projectLSI(simMat, LSI = lsi_model, verbose = verbose))
      rownames(projRD) <- paste0(colnames(shared_rd$mat)[x], "#", seq_len(nrow(projRD)))
      projRD
    }, mc.cores = threads) %>% Reduce("rbind", .)

    # Deal with NaN
    if (any(is.nan(projRD))) {
      projRD[is.nan(projRD)] <- 0
      warning("NaN calculated during single cell generation")
    }

    if (scale) {
      projRD <- scale_dims(projRD)
    }
  } else {
    projRD <- as.matrix(projectLSI(shared_rd$mat, LSI = lsi_model, verbose = verbose))
  }

  # Check LSI and Embedding SVD columns
  if (embedding_num_dim != ncol(projRD)) {
    stop("Error inconsistency found with matching LSI dimensions to those used in embedding")
  }

  # Get Previous UMAP Model
  if (object_type_projector == "monocle3") {
    umap_model_file <- projector@int_metadata$UMAP_model_file
  } else if (object_type_projector == "seurat") {
    if (!is.null(projector@reductions[[embedding]]@misc$model)) {
      umap_model <- projector@reductions[[embedding]]@misc$model
    } else {
      stop("UMAP model not found in projector object. Please run RunUMAP with return.model = TRUE.")
    }
  }

  # Subsample cells for co-embedding
  idx <- sort(sample(seq_len(nrow(rD)), min(nrow(rD), ncells_coembedding)))
  rD_ss <- rD[idx, , drop = FALSE]

  # Project UMAP
  if (verbose & make_pseudo_single_cells) message(paste0("Projecting simulated cells onto manifold"))
  if (verbose & !make_pseudo_single_cells) message(paste0("Projecting projectee cells onto manifold"))
  set.seed(seed)

  if (object_type_projector == "monocle3") {
    # Load UMAP model
    umap_model <- load_umap_model(umap_model_file, embedding_num_dim)
  }
  # Combine the existing reduced dimensions with the projected ones
  combined_rD <- rbind(rD_ss, projRD)

  # Project into UMAP space
  simUMAP <- uwot::umap_transform(
    X = combined_rD,
    model = umap_model,
    verbose = verbose,
    n_threads = threads
  )
  rownames(simUMAP) <- c(rownames(rD_ss), rownames(projRD))

  # Check correlation of subsampled cells
  c1 <- cor(simUMAP[rownames(rD_ss), 1], sc_embedding[rownames(rD_ss), 1])
  c2 <- cor(simUMAP[rownames(rD_ss), 2], sc_embedding[rownames(rD_ss), 2])
  if (min(c1, c2) < 0.8) {
    message(paste0("Warning: projection correlation is less than 0.8 (R = ", round(min(c1, c2), 4), ").\nThese results may not be accurate because of the lack of heterogeneity in the single-cell data."))
  }

  dfUMAP <- sc_embedding
  colnames(dfUMAP) <- c("UMAP1", "UMAP2")
  colnames(simUMAP) <- c("UMAP1", "UMAP2")
  dfUMAP <- DataFrame(dfUMAP)
  dfUMAP$Type <- Rle("single_cell", lengths = nrow(dfUMAP))

  simUMAP <- DataFrame(simUMAP[rownames(projRD), , drop = FALSE])
  simUMAP$Type <- Rle(stringr::str_split(rownames(simUMAP), pattern = "#", simplify = TRUE)[, 1])

  out <- SimpleList(
    projectedUMAP = simUMAP,
    singleCellUMAP = dfUMAP,
    projectedReducedDims = projRD
  )
  return(out)
}


#' Extract Data Using Specified Features (Compatible with Seurat Objects)
#'
#' This function extracts data from a `SummarizedExperiment`-like object or a Seurat object using features specified in the `query`. Features can be identified by either name (annotation-based) or genomic ranges (range-based).
#'
#' @param query A vector of feature names (for annotation-based search) or a `GRanges` object (for range-based search).
#' @param subject A `SummarizedExperiment` or Seurat object from which to extract the data.
#' @param annotation Character, specifies the target feature column in `subject` to search for features. Default is "row.names".
#' @param duplicate_hits Character, specifies how to resolve multiple hits in range-based searches. Options are: "max.mean", "max.var", "max.disp", "min.mean", "min.var", "min.disp". Default is "max.disp".
#' @param fill_output Logical, whether to return the `subject` object filled with zeros for features in `query` that are not found in `subject`. Default is `TRUE`.
#' @param ignore_strand Logical, whether to ignore strand information in range-based searches. Default is `TRUE`.
#' @param verbose Logical, whether to print verbose output. Default is `TRUE`.
#'
#' @return A list containing:
#' - `mat`: A matrix of the extracted data.
#' - `overlap`: The ratio of features found in the `subject` compared to the `query`.
#' - `not_found`: Features not found in the `subject`.
#' @export

extract_data <- function(query,
                         subject,
                         annotation = "row.names",
                         verbose = TRUE,
                         fill_output = TRUE,
                         duplicate_hits = "max.disp",
                         ignore_strand = TRUE) {
  ##################################################
  # Check Inputs
  ##################################################
  # Ensure required packages are loaded
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("Package 'Matrix' is required but not installed.")
  }
  if (!requireNamespace("GenomicRanges", quietly = TRUE)) {
    stop("Package 'GenomicRanges' is required but not installed.")
  }
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Package 'data.table' is required but not installed.")
  }
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package 'dplyr' is required but not installed.")
  }
  if (!requireNamespace("matrixStats", quietly = TRUE)) {
    stop("Package 'matrixStats' is required but not installed.")
  }
  if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
    stop("Package 'SummarizedExperiment' is required but not installed.")
  }
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Package 'Seurat' is required but not installed.")
  }

  # Check if subject is a SummarizedExperiment or Seurat object
  if (methods::is(subject, "SummarizedExperiment")) {
    subject_type <- "SummarizedExperiment"
  } else if (methods::is(subject, "Seurat")) {
    subject_type <- "Seurat"
  } else {
    stop("The 'subject' must be a SummarizedExperiment or Seurat object.")
  }

  # Check if query is a character vector or GRanges
  if (is.character(query)) {
    query_type <- "character"
  } else if (methods::is(query, "GRanges")) {
    query_type <- "GRanges"
  } else {
    stop("The 'query' must be a character vector or a GRanges object.")
  }

  ##################################################
  # Parse duplicate_hits behavior
  ##################################################
  # Parse the function and variable to use for resolving duplicate hits
  duplicate_parts <- strsplit(duplicate_hits, "\\.")[[1]]
  if (length(duplicate_parts) != 2) {
    stop("Invalid 'duplicate_hits' format. Should be in the form 'function.variable', e.g., 'max.mean'.")
  }
  FUN <- get(duplicate_parts[1])
  Var1 <- duplicate_parts[2]

  ##################################################
  # Extract Data Based on Query Type
  ##################################################
  if (query_type == "GRanges") {
    # Range-based search
    if (subject_type == "SummarizedExperiment") {
      # Ensure subject has rowRanges
      if (is.null(SummarizedExperiment::rowRanges(subject))) {
        stop("The 'subject' does not have rowRanges for range-based querying.")
      }
      # Get the overlaps
      se_sub <- GenomicRanges::subsetByOverlaps(subject, query, ignore.strand = ignore_strand, type = "any")
      if (nrow(se_sub) == 0) {
        stop("No overlap between query and subject found.")
      }
      # Get the data matrix
      dat <- SummarizedExperiment::assay(se_sub)
      # Find overlaps and resolve duplicates
      hits <- data.table::as.data.table(GenomicRanges::findOverlaps(query, SummarizedExperiment::rowRanges(se_sub), ignore.strand = ignore_strand))
    } else if (subject_type == "Seurat") {
      # Seurat object does not support range-based querying unless feature ranges are available
      stop("Range-based querying is not supported for Seurat objects without feature genomic ranges.")
    }

    # Continue processing overlaps
    fs <- names(query) <- paste("f", seq_along(query), sep = "_")
    # Calculate statistics for duplicate resolution
    hits$var <- matrixStats::rowVars(dat)[hits$subjectHits]
    hits$mean <- Matrix::rowMeans(dat)[hits$subjectHits]
    hits$disp <- sqrt(hits$var) / hits$mean * 100
    # Resolve duplicates
    best_hits <- hits %>%
      dplyr::group_by(queryHits) %>%
      dplyr::filter(get(Var1) == FUN(get(Var1)))
    if (nrow(best_hits) == 0) {
      stop("No overlap between query and subject found after resolving duplicates.")
    }
    mat <- dat[best_hits$subjectHits, ]
    rownames(mat) <- names(query)[best_hits$queryHits]
    if (fill_output) {
      mat <- safe_subset(mat, subsetRows = fs)
    }
    overlap <- nrow(best_hits) / length(fs)
    not_found <- query[!best_hits$queryHits %in% seq_along(fs)]
  } else {
    # Annotation-based search
    if (annotation == "row.names") {
      fs <- query
      if (subject_type == "SummarizedExperiment") {
        # Use rownames from SummarizedExperiment
        subject_features <- rownames(subject)
        if (is.null(subject_features)) {
          stop("The 'subject' does not have row names.")
        }
        fidx <- which(subject_features %in% fs)
        if (length(fidx) == 0) {
          stop("No overlap between query and subject found.")
        }
        dat <- SummarizedExperiment::assay(subject)
        mat <- dat[fidx, , drop = FALSE]
        rownames(mat) <- subject_features[fidx]
      } else if (subject_type == "Seurat") {
        # Use features from Seurat object
        default_assay <- Seurat::DefaultAssay(subject)
        dat <- Seurat::GetAssayData(subject, assay = default_assay, slot = "data")
        subject_features <- rownames(dat)
        fidx <- which(subject_features %in% fs)
        if (length(fidx) == 0) {
          stop("No overlap between query and subject found.")
        }
        mat <- dat[fidx, , drop = FALSE]
        rownames(mat) <- subject_features[fidx]
      }
      if (fill_output) {
        mat <- safe_subset(mat, subsetRows = fs)
      }
      overlap <- length(fidx) / length(fs)
      not_found <- fs[!fs %in% subject_features]
    } else {
      stop("Annotation-based search with custom annotations is not implemented.")
    }
  }

  if (verbose) {
    message(paste("Extracted data with overlap ratio:", round(overlap, 3)))
    if (length(not_found) > 0) {
      message(paste("Number of features not found:", length(not_found)))
    }
  }

  return(list(mat = mat, overlap = overlap, not_found = not_found))
}





#' Check Input Validity
#'
#' This is an internal helper function that checks the validity of input parameters.
#' @param input The input object to validate.
#' @param name The name of the parameter being checked.
#' @param valid A vector of valid types/classes for the input.
#' @export
#' @keywords internal

check_input<-function (input = NULL, name = NULL, valid = NULL)
{
  valid <- unique(valid)
  if (is.character(valid)) {
    valid <- tolower(valid)
  }
  else {
    stop("Validator must be a character!")
  }
  if (!is.character(name)) {
    stop("name must be a character!")
  }
  if ("null" %in% tolower(valid)) {
    valid <- c("null", valid[which(tolower(valid) != "null")])
  }
  av <- FALSE
  for (i in seq_along(valid)) {
    vi <- valid[i]
    if (vi == "integer" | vi == "wholenumber") {
      if (all(is.numeric(input))) {
        cv <- min(abs(c(input%%1, input%%1 - 1))) < .Machine$double.eps^0.5
      }
      else {
        cv <- FALSE
      }
    }
    else if (vi == "null") {
      cv <- is.null(input)
    }
    else if (vi == "bool" | vi == "boolean" | vi == "logical") {
      cv <- is.logical(input)
    }
    else if (vi == "numeric") {
      cv <- is.numeric(input)
    }
    else if (vi == "vector") {
      cv <- is.vector(input)
    }
    else if (vi == "matrix") {
      cv <- is.matrix(input)
    }
    else if (vi == "sparsematrix") {
      cv <- is(input, "dgCMatrix")
    }
    else if (vi == "character") {
      cv <- is.character(input)
    }
    else if (vi == "factor") {
      cv <- is.factor(input)
    }
    else if (vi == "cell_data_set") {
      cv <- is(input, "cell_data_set")
    }
    else if (vi == "rlecharacter") {
      cv1 <- is(input, "Rle")
      if (cv1) {
        cv <- is(input@values, "factor") || is(input@values,
                                               "character")
      }
      else {
        cv <- FALSE
      }
    }
    else if (vi == "palette") {
      cv <- all(.isColor(input))
    }
    else if (vi == "timestamp") {
      cv <- is(input, "POSIXct")
    }
    else if (vi == "dataframe" | vi == "data.frame" | vi ==
             "df") {
      cv1 <- is.data.frame(input)
      cv2 <- is(input, "DataFrame")
      cv <- any(cv1, cv2)
    }
    else if (vi == "fileexists") {
      cv <- all(file.exists(input))
    }
    else if (vi == "direxists") {
      cv <- all(dir.exists(input))
    }
    else if (vi == "granges" | vi == "gr") {
      cv <- is(input, "GRanges")
    }
    else if (vi == "grangeslist" | vi == "grlist") {
      cv <- .isGRList(input)
    }
    else if (vi == "list" | vi == "simplelist") {
      cv1 <- is.list(input)
      cv2 <- is(input, "SimpleList")
      cv <- any(cv1, cv2)
    }
    else if (vi == "bsgenome") {
      cv1 <- is(input, "BSgenome")
      cv2 <- tryCatch({
        library(input)
        eval(parse(text = input))
      }, error = function(e) {
        FALSE
      })
      cv <- any(cv1, cv2)
    }
    else if (vi == "se" | vi == "summarizedexperiment") {
      cv <- is(input, "SummarizedExperiment")
    }
    else if (vi == "seurat" | vi == "seuratobject") {
      cv <- is(input, "Seurat")
    }
    else if (vi == "txdb") {
      cv <- is(input, "TxDb")
    }
    else if (vi == "orgdb") {
      cv <- is(input, "OrgDb")
    }
    else if (vi == "bsgenome") {
      cv <- is(input, "BSgenome")
    }
    else if (vi == "parallelparam") {
      cv <- is(input, "BatchtoolsParam")
    }
    else {
      stop("Validator is not currently supported")
    }
    if (cv) {
      av <- TRUE
      break
    }
  }
  if (av) {
    return(invisible(TRUE))
  }
  else {
    stop("Input value for '", name, "' is not a ", paste(valid,
                                                         collapse = ","), ", (", name, " = ", class(input),
         ") please supply valid input!")
  }
}

#' Subset Matrix Safely
#'
#' Safely subsets a matrix, filling in missing rows or columns with zeros if necessary.
#' @param mat A matrix to subset.
#' @param subsetRows Character vector specifying row names to subset.
#' @param subsetCols Character vector specifying column names to subset.
#' @return A subsetted matrix, potentially filled with zeros.
#' @export
#' @keywords internal
safe_subset<-function (mat = NULL, subsetRows = NULL, subsetCols = NULL)
{
  if (!is.null(subsetRows)) {
    idxNotIn <- which(!subsetRows %in% rownames(mat))
    if (length(idxNotIn) > 0) {
      subsetNamesNotIn <- subsetRows[idxNotIn]
      matNotIn <- Matrix::sparseMatrix(i = 1, j = 1, x = 0,
                                       dims = c(length(idxNotIn), ncol = ncol(mat)))
      dim(matNotIn)
      dim(mat)
      rownames(matNotIn) <- subsetNamesNotIn
      mat <- rbind(mat, matNotIn)
    }
    mat <- mat[subsetRows, ]
  }
  if (!is.null(subsetCols)) {
    idxNotIn <- which(subsetCols %ni% colnames(mat))
    if (length(idxNotIn) > 0) {
      subsetNamesNotIn <- subsetCols[idxNotIn]
      matNotIn <- Matrix::sparseMatrix(i = 1, j = 1, x = 0,
                                       dims = c(nrow(mat), ncol = length(idxNotIn)))
      colnames(matNotIn) <- subsetNamesNotIn
      mat <- cbind(mat, matNotIn)
    }
    mat <- mat[, subsetCols]
  }
  mat
}

#' get_assay helper
#'
#' @description
#' Adapted from: Jeffrey M. Granja, M. Ryan Corces, Sarah E. Pierce, S. Tansu Bagdatli, Hani Choudhry, Howard Y. Chang, William J. Greenleaf
#' doi: https://doi.org/10.1101/2020.04.28.066498
#'
#' @export
get_assay<-function (se = NULL, assayName = NULL)
{
  .assayNames <- function(se) {
    names(SummarizedExperiment::assays(se))
  }
  if (is.null(assayName)) {
    o <- SummarizedExperiment::assay(se)
  }
  else if (assayName %in% .assayNames(se)) {
    o <- SummarizedExperiment::assays(se)[[assayName]]
  }
  else {
    stop(sprintf("assayName '%s' is not in assayNames of se : %s",
                 assayName, paste(.assayNames(se), collapse = ", ")))
  }
  return(o)
}

#' @export
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}


#' @keywords internal
#' @export
projectLSI<-function (mat = NULL, LSI = NULL, returnModel = FALSE, verbose = FALSE, seed=2020)
{
  out2 <- tryCatch({
    require(Matrix)
    set.seed(seed)
    if(verbose) message(sprintf("Projecting LSI, Input Matrix = %s GB",
                                round(object.size(mat)/10^9, 3)))
    if(verbose) message("Subsetting by Non-Zero features in inital Matrix")
    #mat <- mat[LSI$idx, ] I don't see what this does after safeSubset
    if (LSI$binarize) {
      if(verbose) message("Binarizing Matrix")
      mat@x[mat@x > 0] <- 1
    }
    if(verbose) message("Computing Term Frequency")
    colSm <- Matrix::colSums(mat)
    if (any(colSm == 0)) {
      exclude <- which(colSm == 0)
      mat <- mat[, -exclude]
      colSm <- colSm[-exclude]
    }
    mat@x <- mat@x/rep.int(colSm, Matrix::diff(mat@p))
    if (LSI$LSI_method == 1) {
      #Adapted from Cusanovich et al.

      if(verbose) message("Computing Inverse Document Frequency")
      idf   <- as(log(1 + nrow(LSI$svd$v) / LSI$row_sums), "sparseVector")
      if(verbose) message("Computing TF-IDF Matrix")
      mat <- as(Matrix::Diagonal(x = as.vector(idf)), "sparseMatrix") %*%
        mat
    }
    else if (LSI$LSI_method == 2) {
      #Adapted from Stuart et al.
      if(verbose) message("Computing Inverse Document Frequency")
      idf   <- as( nrow(LSI$svd$v) / LSI$row_sums, "sparseVector")
      if(verbose) message("Computing TF-IDF Matrix")
      mat <- as(Matrix::Diagonal(x = as.vector(idf)), "sparseMatrix") %*%
        mat
      mat@x <- log(mat@x * LSI$scale_to + 1)
    }else if (LSI$LSI_method == 3) {
      mat@x <- log(mat@x + 1)
      if(verbose) message("Computing Inverse Document Frequency")
      idf <- as(log(1 + nrow(LSI$svd$v) /LSI$row_sums), "sparseVector")
      if(verbose) message("Computing TF-IDF Matrix")
      mat <- as(Matrix::Diagonal(x = as.vector(idf)), "sparseMatrix") %*%
        mat
    }else {
      stop("LSIMethod unrecognized please select valid method!")
    }
    gc()
    idxNA <- Matrix::which(is.na(mat), arr.ind = TRUE)
    if (length(idxNA) > 0) {
      if(verbose) message((sprintf("Zeroing %s NA elements", length(idxNA))))
      mat[idxNA] <- 0
    }
    if(verbose) message("Calculating V Matrix")
    V <- Matrix::t(mat) %*% LSI$svd$u %*% Matrix::diag(1/LSI$svd$d)
    if(verbose) message("Computing Projected Coordinates")
    svdDiag <- matrix(0, nrow = LSI$num_dim, ncol = LSI$num_dim)
    diag(svdDiag) <- LSI$svd$d
    matSVD <- Matrix::t(svdDiag %*% Matrix::t(V))
    matSVD <- as.matrix(matSVD)
    rownames(matSVD) <- colnames(mat)
    colnames(matSVD) <- paste0("LSI", seq_len(ncol(matSVD)))
    if (returnModel) {
      if(verbose) message("Calculating Re-Projected Matrix")
      X <- LSI$svd$u %*% diag(LSI$svd$d) %*% t(V)
      out <- list(matSVD = matSVD, V = V, X = X)
    }
    else {
      out <- matSVD
    }
    out
  }, error = function(e) {
    errorList <- list(mat = mat, colSm = if (exists("colSm",
                                                    inherits = FALSE)) colSm else "Error with colSm!",
                      idf = if (exists("idf", inherits = FALSE)) idf else "Error with idf!",
                      V = if (exists("V", inherits = FALSE)) V else "Error with V!",
                      matSVD = if (exists("matSVD", inherits = FALSE)) matSVD else "Error with matSVD!")
    errorLog(e, fn = "projectLSI", info = "", errorList = errorList)
  })
  out2
}

#' @keywords internal
#' @export
scale_dims<-function(x, scale_max = NULL){
  if(!is.null(scale_max)){
    row_z_scores(m=x, min=-scale_max, max = scale_max, limit = TRUE)
  }else{
    row_z_scores(m=x)
  }
}

#' @keywords internal
#' @export
row_z_scores<-function(m = NULL, min = -2, max = 2, limit = FALSE){
  z <- sweep(m - Matrix::rowMeans(m), 1, matrixStats::rowSds(m),`/`)
  if(limit){
    z[z > max] <- max
    z[z < min] <- min
  }
  return(z)
}

#' @keywords internal
#' @export
errorLog<-function(
    e = NULL,
    fn = NULL,
    info = NULL,
    errorList = NULL,
    throwError = TRUE
){

  header <- "************************************************************"
  if(!is.null(errorList)){
    tryCatch({
      saveRDS(errorList, "Save-Error.rds")
      message("Saving a list of errors to Save-Error.rds")
    }, error = function(e){
      message("Error recording errorList")
    })
  }
  print(e)
  cat(sprintf("\n%s\n\n", header))
  if(throwError) stop("Exiting See Error Above")
}


#' Helper function for summing sparse matrix groups
#' @references Granja, J. M.et al. (2019). Single-cell multiomic analysis identifies regulatory programs in mixed-phenotype
#' acute leukemia. Nature Biotechnology, 37(12), 1458–1465.
#' @export
#' @keywords internal
#' @keywords internal
sparseRowVariances <- function (m){
  rM <- Matrix::rowMeans(m)
  rV <- computeSparseRowVariances(m@i + 1, m@x, rM, ncol(m))
  return(rV)
}
#
#
# #Sparse Variances Rcpp
# Rcpp::sourceCpp(code='
#   #include <Rcpp.h>
#   using namespace Rcpp;
#   using namespace std;
#   // [[Rcpp::export]]
#   Rcpp::NumericVector computeSparseRowVariances(IntegerVector j, NumericVector val, NumericVector rm, int n) {
#     const int nv = j.size();
#     const int nm = rm.size();
#     Rcpp::NumericVector rv(nm);
#     Rcpp::NumericVector rit(nm);
#     int current;
#     // Calculate RowVars Initial
#     for (int i = 0; i < nv; ++i) {
#       current = j(i) - 1;
#       rv(current) = rv(current) + (val(i) - rm(current)) * (val(i) - rm(current));
#       rit(current) = rit(current) + 1;
#     }
#     // Calculate Remainder Variance
#     for (int i = 0; i < nm; ++i) {
#       rv(i) = rv(i) + (n - rit(i))*rm(i)*rm(i);
#     }
#     rv = rv / (n - 1);
#     return(rv);
#   }'
# )

#New save_umap
#' @keywords internal
#' @export
save_umap_model <- function(model, file){
  if (!is.null(model$nn_index$ann)) {
    save_umap_model_new(model = model, file = file)
  }
  else {
    save_umap_model_depracated(model = model, file = file) #backwards to previous version
  }
}

#' @keywords internal
#' @export
save_umap_model_new<-function(model, file){
  file <- file.path(normalizePath(dirname(file)), basename(file))
  wd <- getwd()
  mod_dir <- tempfile(pattern = "dir")
  dir.create(mod_dir)
  uwot_dir <- file.path(mod_dir, "uwot")
  dir.create(uwot_dir)
  model_tmpfname <- file.path(uwot_dir, "model")
  saveRDS(model, file = model_tmpfname)
  metrics <- names(model$metric)
  n_metrics <- length(metrics)
  for (i in seq_len(n_metrics)) {
    nn_tmpfname <- file.path(uwot_dir, paste0("nn", i))
    if (n_metrics == 1) {
      model$nn_index$ann$save(nn_tmpfname)
      model$nn_index$ann$unload()
      model$nn_index$ann$load(nn_tmpfname)
    }
    else {
      model$nn_index[[i]]$ann$save(nn_tmpfname)
      model$nn_index[[i]]$ann$unload()
      model$nn_index[[i]]$ann$load(nn_tmpfname)
    }
  }
  setwd(mod_dir)
  system2("tar", "-cvf uwot.tar uwot", stdout = NULL, stderr = NULL)
  o <- file_rename("uwot.tar", file)
  setwd(wd)
  if (file.exists(mod_dir)) {
    unlink(mod_dir, recursive = TRUE)
  }
  return(o)
}


#' @keywords internal
#' @export
save_umap_model_depracated<-function(model, file){
  file <- file.path(normalizePath(dirname(file)), basename(file))
  wd <- getwd()
  mod_dir <- tempfile(pattern = "dir")
  dir.create(mod_dir)
  uwot_dir <- file.path(mod_dir, "uwot")
  dir.create(uwot_dir)
  model_tmpfname <- file.path(uwot_dir, "model")
  saveRDS(model, file = model_tmpfname)
  metrics <- names(model$metric)
  n_metrics <- length(metrics)
  for (i in seq_len(n_metrics)) {
    nn_tmpfname <- file.path(uwot_dir, paste0("nn", i))
    if (n_metrics == 1) {
      model$nn_index$save(nn_tmpfname)
      model$nn_index$unload()
      model$nn_index$load(nn_tmpfname)
    }
    else {
      model$nn_index[[i]]$save(nn_tmpfname)
      model$nn_index[[i]]$unload()
      model$nn_index[[i]]$load(nn_tmpfname)
    }
  }
  setwd(mod_dir)
  system2("tar", "-cvf uwot.tar uwot", stdout = NULL, stderr = NULL)
  o <- file_rename("uwot.tar", file)
  setwd(wd)
  if (file.exists(mod_dir)) {
    unlink(mod_dir, recursive = TRUE)
  }
  return(o)
}

#' @keywords internal
#' @export
file_rename<-function(from = NULL, to = NULL){

  if(!file.exists(from)){
    stop("Input file does not exist!")
  }

  tryCatch({

    o <- .suppressAll(file.rename(from, to))

    if(!o){
      stop("retry with mv")
    }

  }, error = function(x){

    tryCatch({

      system(paste0("mv '", from, "' '", to, "'"))

      return(to)

    }, error = function(y){

      stop("File Moving/Renaming Failed!")

    })

  })

}


#New load_umap_model
#' @export
#' @keywords internal
#'
load_umap_model <- function(file, num_dim = NULL){
  tryCatch({
    load_umap_model_new(file = file, num_dim = num_dim)
  }, error = function(x){
    load_umap_model_depracated(file = file, num_dim = num_dim)
  })
}

#' @export
#' @keywords internal
load_umap_model_new<-function(file, num_dim = NULL){
  model <- NULL
  tryCatch({
    mod_dir <- tempfile(pattern = "dir")
    dir.create(mod_dir)
    utils::untar(file, exdir = mod_dir)
    model_fname <- file.path(mod_dir, "uwot/model")
    if (!file.exists(model_fname)) {
      stop("Can't find model in ", file)
    }
    model <- readRDS(file = model_fname)
    metrics <- names(model$metric)
    n_metrics <- length(metrics)
    for (i in seq_len(n_metrics)){
      nn_fname <- file.path(mod_dir, paste0("uwot/nn", i))
      if (!file.exists(nn_fname)) {
        stop("Can't find nearest neighbor index ", nn_fname, " in ", file)
      }
      metric <- metrics[[i]]
      if(length(model$metric[[i]]) == 0){
        if(!is.null(num_dim)){
          num_dim2 <- num_dim
        }else{
          num_dim2 <- length(model$metric[[i]])
        }
      }
      if(!is.null(num_dim)){
        num_dim2 <- num_dim
      }
      ann <- uwot:::create_ann(metric, ndim = num_dim2)
      ann$load(nn_fname)
      if (n_metrics == 1) {
        model$nn_index$ann <- ann
      }else{
        model$nn_index[[i]]$ann <- ann
      }
    }
  }, finally = {
    if (file.exists(mod_dir)) {
      unlink(mod_dir, recursive = TRUE)
    }
  })
  model
}

#' @keywords internal
#' @export
load_umap_model_depracated<-function(file, num_dim = NULL){
  model <- NULL
  tryCatch({
    mod_dir <- tempfile(pattern = "dir")
    dir.create(mod_dir)
    utils::untar(file, exdir = mod_dir)
    model_fname <- file.path(mod_dir, "uwot/model")
    if (!file.exists(model_fname)) {
      stop("Can't find model in ", file)
    }
    model <- readRDS(file = model_fname)
    metrics <- names(model$metric)
    n_metrics <- length(metrics)
    for (i in seq_len(n_metrics)){
      nn_fname <- file.path(mod_dir, paste0("uwot/nn", i))
      if (!file.exists(nn_fname)) {
        stop("Can't find nearest neighbor index ", nn_fname, " in ", file)
      }
      metric <- metrics[[i]]
      if(length(model$metric[[i]]) == 0){
        if(!is.null(num_dim)){
          num_dim2 <- num_dim
        }else{
          num_dim2 <- length(model$metric[[i]])
        }
      }
      if(!is.null(num_dim)){
        num_dim2 <- num_dim
      }
      ann <- uwot:::create_ann(metric, ndim = num_dim2)
      ann$load(nn_fname)
      if (n_metrics == 1) {
        model$nn_index <- ann
      }else{
        model$nn_index[[i]] <- ann
      }
    }
  }, finally = {
    if (file.exists(mod_dir)) {
      unlink(mod_dir, recursive = TRUE)
    }
  })
  model
}

#' Liftover projection
#'
#' @description This function transfers labels from a projector to a projectee
#' @param projection projection object created by project_data function
#' @param projector cell_data_set object with a reduced dimension matrix (currently LSI supported) as specified in
#' reduced_dim argument and a model used to create a low dimensional embedding
#' @param projectee a SummarizedExperiment type object (cell_data_set currently supported) to be projected using
#' the models contained in the projector
#' @param projector_col column name of projector to color cells by using nearest neighbor
#' @param projectee_col column name of projectee to color cells by
#' @param colors color palatte of projectee_col cells
#' @export
#' @keywords internal
liftover_projection<-function(projection, projector, projectee, projector_col=NULL){
  closest<-RANN::nn2(data = data.frame(UMAP1=projection[[2]]$UMAP1, UMAP2=projection[[2]]$UMAP2), query = data.frame(UMAP1=projection[[1]]$UMAP1, UMAP2=projection[[1]]$UMAP2), k = 1)
  projectee$liftover<-projector[[projector_col]][closest$nn.idx]
  projectee
}


#' Plot Projection of Projected Data onto Projector's UMAP Embedding (Compatible with Monocle3 and Seurat)
#'
#' @description This function visualizes the projection of data from one single-cell object onto the UMAP embedding of another object. It allows for coloring the projected data points based on metadata from either the projector or projectee object.
#'
#' @param projection A list or `SimpleList` containing the UMAP coordinates of the projected data. Typically, this is the output from the `project_data` function.
#' @param projector The original single-cell object (Monocle3 `cell_data_set` or Seurat object) used to generate the UMAP embedding.
#' @param projectee The single-cell object (Monocle3 `cell_data_set` or Seurat object) that was projected into the UMAP space.
#' @param projector_col Optional character string specifying the column name in the projector's metadata to use for coloring the data points.
#' @param projectee_col Optional character string specifying the column name in the projectee's metadata to use for coloring the data points.
#'
#' @return A `ggplot` object visualizing the projection.
#' @export
#'
#' @examples
#' \dontrun{
#' # Assuming 'projection' is the output from project_data()
#' plot_projection(
#'   projection = projection,
#'   projector = projector_object,
#'   projectee = projectee_object,
#'   projector_col = "cell_type",
#'   projectee_col = "sample_id"
#' )
#' }
plot_projection <- function(projection, projector, projectee, projector_col = NULL, projectee_col = NULL) {
  # Ensure 'projection' is a list
  projection <- as.list(projection)

  # Find nearest neighbors in UMAP space
  closest <- RANN::nn2(
    data = data.frame(UMAP1 = projection[[2]]$UMAP1, UMAP2 = projection[[2]]$UMAP2),
    query = data.frame(UMAP1 = projection[[1]]$UMAP1, UMAP2 = projection[[1]]$UMAP2),
    k = 1
  )

  # Assign labels to the projector data
  projection[[2]]$neighbor <- "00_PROJECTOR"

  # Assign labels to the projectee data
  if (!is.null(projectee_col)) {
    # Get labels from projectee
    if (methods::is(projectee, "Seurat")) {
      labels <- projectee@meta.data[[projectee_col]]
    } else if (methods::is(projectee, "cell_data_set")) {
      labels <- colData(projectee)[[projectee_col]]
    } else {
      stop("Projectee must be a Seurat object or a Monocle3 cell_data_set object.")
    }
    # Assign labels to projection[[1]]
    projection[[1]]$neighbor <- labels
  } else if (!is.null(projector_col)) {
    # Get labels from projector's nearest neighbors
    if (methods::is(projector, "Seurat")) {
      labels <- projector@meta.data[[projector_col]][closest$nn.idx]
    } else if (methods::is(projector, "cell_data_set")) {
      labels <- colData(projector)[[projector_col]][closest$nn.idx]
    } else {
      stop("Projector must be a Seurat object or a Monocle3 cell_data_set object.")
    }
    # Assign labels to projection[[1]]
    projection[[1]]$neighbor <- labels
  } else {
    # Default label if no metadata column is provided
    projection[[1]]$neighbor <- "projectee"
  }

  # Combine the projections
  combined_data <- as.data.frame(do.call(rbind, projection[2:1]))

  # Plot using ggplot2
  g <- ggplot(combined_data, aes(x = UMAP1, y = UMAP2, color = neighbor)) +
    geom_point(size = 0.3) +
    theme_minimal() +
    labs(title = "Projection onto UMAP Embedding", color = "Group")

  return(g)
}


