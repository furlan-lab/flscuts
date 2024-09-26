#' @title Iterative Latent Semantic Indexing (LSI) for Single-Cell Data (Compatible with Monocle3 and Seurat) with optional UMAP
#'
#' @description
#' This function performs iterative LSI on single-cell data to minimize batch effects and accentuate cell type differences. It can accept either a Monocle3 `cell_data_set` object or a Seurat object as input. The function iteratively:
#' 1. Applies TF-IDF transformation and Singular Value Decomposition (SVD) to normalize the data.
#' 2. Clusters the normalized data using Leiden clustering in high-dimensional space.
#' 3. Identifies over-represented features in the resulting clusters using a simple counting method.
#'
#' These steps are repeated, using features identified in step 3 to subset the normalization matrix in step 1, and the process is repeated for a specified number of iterations. This method is inspired by Granja et al. (2019) and aims to enhance the separation of cell types while reducing batch effects.
#'
#' @param object A Monocle3 `cell_data_set` object or a Seurat object.
#' @param num_dim Integer specifying the number of principal components to use in downstream analysis. Default is 25.
#' @param starting_features Optional character vector of starting features (e.g., genes or peaks) to use in the first iteration.
#' @param resolution Numeric vector specifying the resolution parameters for Leiden clustering at each iteration. The number of iterations is determined by the length of this vector.
#' @param num_features Integer or numeric vector specifying the number of features to use for dimensionality reduction at each iteration. If a single integer is provided, it is used for all iterations. Default is 3000.
#' @param exclude_features Optional character vector of features (rownames of the data) to exclude from analysis.
#' @param do_tf_idf Logical indicating whether to perform TF-IDF transformation. Default is \code{TRUE}.
#' @param binarize Logical indicating whether to binarize the data prior to TF-IDF transformation. Default is \code{FALSE}.
#' @param scale Logical indicating whether to scale the data to \code{scale_to}. Default is \code{TRUE}.
#' @param log_transform Logical indicating whether to log-transform the data after scaling. Default is \code{TRUE}.
#' @param scale_to Numeric value specifying the scaling factor if \code{scale} is \code{TRUE}. Default is 10000.
#' @param leiden_k Integer specifying the number of nearest neighbors (k) for Leiden clustering. Default is 20.
#' @param leiden_weight Logical indicating whether to use edge weights in Leiden clustering. Default is \code{FALSE}.
#' @param leiden_iter Integer specifying the number of iterations for Leiden clustering. Default is 1.
#' @param random_seed Integer specifying the random seed for reproducibility. Default is 2020.
#' @param verbose Logical indicating whether to display progress messages. Default is \code{FALSE}.
#' @param run_umap Logical indicating whether to run UMAP.
#' @param return_object Logical indicating whether to return the updated input object with LSI reduction and clustering results. Default is \code{TRUE}.
#' @param ... Additional arguments passed to lower-level functions.
#'
#' @return If \code{return_object} is \code{TRUE}, returns the updated input object (Monocle3 `cell_data_set` or Seurat object) with LSI reduction and clustering results added. If \code{FALSE}, returns a list with elements:
#' \describe{
#'   \item{\code{lsi_embeddings}}{The final LSI embeddings.}
#'   \item{\code{clusters}}{The clustering assignments.}
#'   \item{\code{iterations}}{A list containing intermediate results from each iteration.}
#' }
#'
#' @details
#' The function performs iterative LSI as described in Granja et al. (2019), adapting methods from Cusanovich et al. (2018). It is suitable for processing single-cell ATAC-seq or RNA-seq data to identify meaningful clusters and reduce batch effects.
#'
#' @references
#' Granja, J. M., et al. (2019). Single-cell multiomic analysis identifies regulatory programs in mixed-phenotype acute leukemia. \emph{Nature Biotechnology}, 37(12), 1458–1465.
#'
#' Cusanovich, D. A., et al. (2018). The cis-regulatory dynamics of embryonic development at single-cell resolution. \emph{Nature}, 555(7697), 538–542.
#'
#' @importFrom methods is
#' @importFrom monocle3 cluster_cells find_partitions clusters
#' @importFrom SingleCellExperiment reducedDims<- rowData
#' @importFrom Seurat GetAssayData CreateDimReducObject FindNeighbors FindClusters Idents
#' @importFrom Matrix colSums rowSums rowMeans
#' @importFrom irlba irlba
#' @importFrom edgeR cpm
#' @importFrom matrixStats rowVars
#' @export
#'
#' @examples
#' \dontrun{
#' # For a Monocle3 cell_data_set object:
#' cds <- iterative_LSI(
#'   object = cds,
#'   num_dim = 30,
#'   resolution = c(1e-4, 3e-4, 5e-4)
#' )
#'
#' # For a Seurat object:
#' seurat_obj <- iterative_LSI(
#'   object = seurat_obj,
#'   num_dim = 30,
#'   resolution = c(0.2, 0.5, 0.8)
#' )
#' }
iterative_LSI <- function (object, num_dim = 25, starting_features = NULL, resolution = c(1e-04,
                                                                                          3e-04, 5e-04), do_tf_idf = TRUE, num_features = c(3000, 3000,
                                                                                                                                            3000), exclude_features = NULL, binarize = FALSE, scale = TRUE,
                           log_transform = TRUE, LSI_method = 1, partition_qval = 0.05,
                           seed = 2020, scale_to = 10000, leiden_k = 20, leiden_weight = FALSE,
                           leiden_iter = 1, verbose = FALSE, return_iterations = FALSE, run_umap = FALS, ...)
{
  # Check object type
  if (is(object, "Seurat")) {
    object_type <- "seurat"
  } else if (is(object, "cell_data_set")) {
    object_type <- "monocle3"
  } else {
    stop("The object must be a Seurat object or a Monocle3 cell_data_set object.")
  }

  # Handle starting_features and num_features length
  if (!is.null(starting_features)) {
    if (length(num_features) != length(resolution)) {
      num_features <- c(length(starting_features), num_features)
    }
  }
  if (length(num_features) != length(resolution)) {
    message("Numbers of elements for resolution and num_features do not match. Will use num_features[1]...")
    num_features <- rep(num_features, length(resolution))
  }

  # Get the expression matrix
  if (!is.null(exclude_features)) {
    if (object_type == "seurat") {
      mat <- GetAssayData(object)
      mat <- mat[!rownames(mat) %in% exclude_features, ]
    } else {
      mat <- assay(object)
      mat <- mat[!rownames(mat) %in% exclude_features, ]
    }
  } else {
    if (object_type == "seurat") {
      mat <- GetAssayData(object)
    } else {
      mat <- assay(object)
    }
  }

  original_features <- rownames(mat)
  set.seed(seed)

  if (binarize) {
    message("Binarizing...")
    mat@x[mat@x > 0] <- 1
  }

  outlist <- list()

  if (scale) {
    matNorm <- t(t(mat)/Matrix::colSums(mat)) * scale_to
  } else {
    matNorm <- mat
  }

  if (log_transform) {
    matNorm@x <- log2(matNorm@x + 1)
  }

  message("Performing LSI/SVD for iteration 1....")

  if (!is.null(starting_features)) {
    if (!all(starting_features %in% rownames(mat))) {
      stop("Not all starting features found in data")
    }
    f_idx <- which(rownames(mat) %in% starting_features)
  } else {
    # Compute variances and select features
    if (requireNamespace("matrixStats", quietly = TRUE)) {
      feature_vars <- matrixStats::rowVars(as.matrix(matNorm))
    } else {
      stop("Package 'matrixStats' is required for variance calculation.")
    }
    f_idx <- head(order(feature_vars, decreasing = TRUE), num_features[1])
  }

  if (do_tf_idf) {
    tf <- tf_idf_transform(mat[f_idx, ], method = LSI_method)
    row_sums <- Matrix::rowSums(mat[f_idx, ])
    tf@x[is.na(tf@x)] <- 0
  } else {
    tf <- mat[f_idx, ]
    row_sums <- Matrix::rowSums(mat[f_idx, ])
  }

  svd_list <- svd_lsi(tf, num_dim, mat_only = FALSE)

  # Perform clustering
  if (object_type == "monocle3") {
    # For monocle3, use monocle3:::leiden_clustering
    cluster_result <- monocle3:::leiden_clustering(data = svd_list$matSVD,
                                                   pd = colData(object), k = leiden_k, weight = leiden_weight,
                                                   num_iter = leiden_iter, resolution_parameter = resolution[1],
                                                   random_seed = seed, verbose = verbose, nn_control = list("method"="nn2"), ...)
    clusters <- factor(igraph::membership(cluster_result$optim_res))
  } else if (object_type == "seurat") {
    # For Seurat, use FindNeighbors and FindClusters
    object@reductions[["lsi"]] <- Seurat::CreateDimReducObject(embeddings = svd_list$matSVD,
                                                    key = "LSI_", assay = Seurat::DefaultAssay(object))
    object <- Seurat::FindNeighbors(object, reduction = "lsi", dims = 1:num_dim, k.param = leiden_k)
    object <- Seurat::FindClusters(object, resolution = resolution[1], algorithm = 1, random.seed = seed, verbose = verbose)
    clusters <- Seurat::Idents(object)
  }

  # Proceed with the rest of the function, adapting as needed
  # For example, calculate clusterMat, etc.
  clusterMat <- edgeR::cpm(groupSums(mat, clusters, sparse = TRUE),
                           log = TRUE, prior.count = 3)

  if (length(resolution) == 1) {
    # Final iteration
    if (object_type == "monocle3") {
      SingleCellExperiment::reducedDims(object)[["LSI"]] <- svd_list$matSVD
      # Save other components as needed
      # Store clusters
      if (length(unique(clusters)) > 1) {
        cluster_graph_res <- monocle3:::compute_partitions(cluster_result$g,
                                                           cluster_result$optim_res, partition_qval, verbose)
        partitions <- igraph::components(cluster_graph_res$cluster_g)$membership[cluster_result$optim_res$membership]
        partitions <- as.factor(partitions)
      } else {
        partitions <- rep(1, nrow(colData(object)))
      }
      names(partitions) <- row.names(colData(object))
      object@clusters[["LSI"]] <- list(cluster_result = cluster_result,
                                       partitions = partitions, clusters = clusters)
      if (run_umap){
        object <- run_umap(object, ...)
      }
    } else if (object_type == "seurat") {
      # Clusters are already stored in object@meta.data
      object@reductions[["lsi"]]@misc <- list(svd=svd_list$svd, features=original_features[f_idx],
                                          row_sums = row_sums, seed=seed, binarize=binarize,
                                          scale_to=scale_to, num_dim=num_dim, resolution=resolution,
                                          granges=NULL, LSI_method=LSI_method, outliers=NULL)
      if (run_umap){
        object <- run_umap(object, ...)
      }

    }

    if (return_iterations) {
      outlist[["iteration_1"]] <- list(matSVD = svd_list$matSVD,
                                       features = original_features[f_idx], clusters = clusters)
      return(list(object = object, iterationlist = outlist))
    } else {
      return(object)
    }
  }

  # For multiple iterations
  for (iteration in 2:length(resolution)) {
    message("Performing LSI/SVD for iteration ", iteration, "....")
    f_idx <- head(order(matrixStats::rowVars(clusterMat), decreasing = TRUE),
                  num_features[iteration])
    if (do_tf_idf) {
      tf <- tf_idf_transform(mat[f_idx, ], method = LSI_method)
      tf@x[is.na(tf@x)] <- 0
      row_sums <- Matrix::rowSums(mat[f_idx, ])
    } else {
      tf <- mat[f_idx, ]
      row_sums <- Matrix::rowSums(mat[f_idx, ])
    }

    svd_list <- svd_lsi(tf, num_dim, mat_only = FALSE)

    # Clustering
    if (object_type == "monocle3") {
      cluster_result <- monocle3:::leiden_clustering(data = svd_list$matSVD,
                                                     pd = colData(object), k = leiden_k, weight = leiden_weight,
                                                     num_iter = leiden_iter, resolution_parameter = resolution[iteration],
                                                     random_seed = seed, verbose = verbose, nn_control = list("method"="nn2"), ...)
      clusters <- factor(igraph::membership(cluster_result$optim_res))
    } else if (object_type == "seurat") {
      object@reductions[["lsi"]] <- Seurat::CreateDimReducObject(embeddings = svd_list$matSVD,
                                                      key = "LSI_", assay = Seurat::DefaultAssay(object))
      object <- Seurat::FindNeighbors(object, reduction = "lsi", dims = 1:num_dim, k.param = leiden_k)
      object <- Seurat::FindClusters(object, resolution = resolution[iteration], algorithm = 1, random.seed = seed, verbose = verbose)
      clusters <- Seurat::Idents(object)
    }

    clusterMat <- edgeR::cpm(groupSums(mat, clusters, sparse = TRUE),
                             log = TRUE, prior.count = 3)

    if (iteration == length(resolution)) {
      # Final iteration
      if (object_type == "monocle3") {
        SingleCellExperiment::reducedDims(object)[["LSI"]] <- svd_list$matSVD
        # Store clusters and partitions
        if (length(unique(clusters)) > 1) {
          cluster_graph_res <- monocle3:::compute_partitions(cluster_result$g,
                                                             cluster_result$optim_res, partition_qval, verbose)
          partitions <- igraph::components(cluster_graph_res$cluster_g)$membership[cluster_result$optim_res$membership]
          partitions <- as.factor(partitions)
        } else {
          partitions <- rep(1, nrow(colData(object)))
        }
        names(partitions) <- row.names(colData(object))
        object@clusters[["LSI"]] <- list(cluster_result = cluster_result,
                                         partitions = partitions, clusters = clusters)
      } else if (object_type == "seurat") {
        # Clusters are already stored
      }
      if (return_iterations) {
        it_count <- paste0("iteration_", iteration)
        outlist[[it_count]] <- list(matSVD = svd_list$matSVD,
                                    features = original_features[f_idx], clusters = clusters)
        return(list(object = object, iterationlist = outlist))
      } else {
        return(object)
      }
    } else {
      if (return_iterations) {
        it_count <- paste0("iteration_", iteration)
        outlist[[it_count]] <- list(matSVD = svd_list$matSVD,
                                    features = original_features[f_idx], clusters = clusters)
      }
      next
    }
  }
}

#' @keywords internal
#' @importFrom uwot umap
#' @export
#'
run_umap <- function(object, ...) {
  # Check object type
  if (is(object, "Seurat")) {
    object_type <- "seurat"
  } else if (is(object, "cell_data_set")) {
    object_type <- "monocle3"
  } else {
    stop("The object must be a Seurat object or a Monocle3 cell_data_set object.")
  }

  # Default UMAP parameters
  default_params <- list(
    n_neighbors = 30L,
    n_components = 2L,
    metric = "cosine",
    n_epochs = NULL,
    learning_rate = 1,
    min_dist = 0.3,
    spread = 1,
    set_op_mix_ratio = 1,
    local_connectivity = 1L,
    repulsion_strength = 1,
    negative_sample_rate = 5,
    verbose = TRUE,
    ret_model = TRUE,
    ret_nn = TRUE
  )

  # Capture additional arguments from ...
  user_params <- list(...)

  # Merge user-provided parameters with the default ones
  umap_params <- modifyList(default_params, user_params)

  if (object_type == "monocle3") {
    ##TODO
    print("Under construction for Monocle3 objects.")
    return(0)
  }

  if (object_type == "seurat") {
    umap_params$X <- object@reductions$lsi@cell.embeddings
    # Run UMAP with the final set of parameters
    umap_res <- do.call(umap, umap_params)
    #umap_res = umap(umap_params$X, ret_model = TRUE,  ret_extra = "model", verbose = umap_params$verbose)
    # Rename UMAP dimensions
    colnames(umap_res$embedding) <- c("UMAP_1", "UMAP_2")
    # Store the UMAP embedding in the Seurat object
    object@reductions[["umap"]] <- Seurat::CreateDimReducObject(embeddings = umap_res$embedding, key = "UMAP_", assay = DefaultAssay(object))
    object@reductions[["umap"]]@misc$model <- umap_res
    # Optionally return the modified object
    return(object)
  }
}


#' Cluster LSI (Compatible with Monocle3 and Seurat)
#'
#' @description This function extracts clustering from the last iteration of LSI (see \code{iterative_LSI}) in a single-cell experiment. It uses Leiden clustering and computes partitions. Compatible with both Monocle3 `cell_data_set` and Seurat objects.
#'
#' @param object The `cell_data_set` or Seurat object upon which to perform this operation.
#' @param k Integer number of nearest neighbors to use when creating the k nearest neighbor graph for Leiden clustering. Default is 20.
#' @param weight A logical argument to determine whether or not to use Jaccard coefficients for two nearest neighbors (based on the overlapping of their kNN) as the weight used for Louvain clustering. Default is FALSE.
#' @param num_iter Integer number of iterations used for Leiden clustering. Default is 1.
#' @param resolution Parameter that controls the resolution of clustering. If NULL (Default), the parameter is determined automatically.
#' @param random_seed The seed used by the random number generator. Default is 2020.
#' @param verbose A logical flag to determine whether or not to print the run details.
#' @param partition_qval Numeric, the q-value cutoff to determine when to partition. Default is 0.05.
#' @references Granja, J. M.et al. (2019). Single-cell multiomic analysis identifies regulatory programs in mixed-phenotype acute leukemia. Nature Biotechnology, 37(12), 1458–1465.
#' @export
#' @keywords internal
cluster_LSI <- function(object,
                        k = 20,
                        weight = FALSE,
                        num_iter = 1,
                        resolution = NULL,
                        random_seed = 2020,
                        verbose = TRUE,
                        partition_qval = 0.05) {

  # Check object type
  if (is(object, "Seurat")) {
    object_type <- "seurat"
  } else if (is(object, "cell_data_set")) {
    object_type <- "monocle3"
  } else {
    stop("The object must be a Seurat object or a Monocle3 cell_data_set object.")
  }

  if (object_type == "monocle3") {
    lsi_embeddings <- SingleCellExperiment::reducedDims(object)[["LSI"]]
    cluster_result <- monocle3:::leiden_clustering(data = lsi_embeddings,
                                                   pd = colData(object), k = k, weight = weight,
                                                   num_iter = num_iter, nn_control = list("method" = "nn2"),
                                                   resolution_parameter = resolution,
                                                   random_seed = random_seed, verbose = verbose)
    cluster_graph_res <- monocle3:::compute_partitions(cluster_result$g,
                                                       cluster_result$optim_res, partition_qval, verbose = verbose)
    partitions <- as.factor(igraph::components(cluster_graph_res$cluster_g)$membership[cluster_result$optim_res$membership])
    clusters <- factor(igraph::membership(cluster_result$optim_res))
    object@clusters[["UMAP"]] <- list(cluster_result = cluster_result,
                                      partitions = partitions, clusters = clusters)
  } else if (object_type == "seurat") {
    object <- Seurat::FindNeighbors(object, reduction = "lsi", dims = 1:ncol(object[["lsi"]]@cell.embeddings), k.param = k)
    object <- Seurat::FindClusters(object, resolution = resolution, algorithm = 4, random.seed = random_seed, verbose = verbose)
  }

  return(object)
}




#' Performs TF-IDF transformation on a cell_data_set
#'
#' @description Just like it sounds.
#'
#' @param cds_list Input cell_data_set object or sparse matrix.
#' @importFrom Matrix rowSums
#' @importFrom Matrix colSums
#' @importFrom Matrix Diagonal
#' @importFrom Matrix t
#' @export
#' @keywords internal
tf_idf_transform <- function(input, method=1, verbose=T){
  if(class(input)=="cell_data_set"){
    mat<-exprs(input)
  }else{
    mat<-input
  }
  rn <- rownames(mat)
  row_sums<-rowSums(mat)
  nz<-which(row_sums>0)
  mat <- mat[nz,]
  rn <- rn[nz]
  row_sums <- row_sums[nz]
  col_sums <- colSums(mat)

  #column normalize
  mat <-Matrix::t(Matrix::t(mat)/col_sums)


  if (method == 1) {
    #Adapted from Casanovich et al.
    if(verbose) message("Computing Inverse Document Frequency")
    idf   <- as(log(1 + ncol(mat) / row_sums), "sparseVector")
    if(verbose) message("Computing TF-IDF Matrix")
    mat <- as(Diagonal(x = as.vector(idf)), "sparseMatrix") %*%
      mat
  }
  else if (method == 2) {
    #Adapted from Stuart et al.
    if(verbose) message("Computing Inverse Document Frequency")
    idf   <- as( ncol(mat) / row_sums, "sparseVector")
    if(verbose) message("Computing TF-IDF Matrix")
    mat <- as(Diagonal(x = as.vector(idf)), "sparseMatrix") %*%
      mat
    mat@x <- log(mat@x * scale_to + 1)
  }else if (method == 3) {
    mat@x <- log(mat@x + 1)
    if(verbose) message("Computing Inverse Document Frequency")
    idf <- as(log(1 + ncol(mat) /row_sums), "sparseVector")
    if(verbose) message("Computing TF-IDF Matrix")
    mat <- as(Diagonal(x = as.vector(idf)), "sparseMatrix") %*%
      mat
  }else {
    stop("LSIMethod unrecognized please select valid method!")
  }
  rownames(mat) <- rn
  if(class(input)=="cell_data_set"){
    input@assays$data$counts<-mat
    return(input)
  }else{
    return(mat)
  }
}


#' @title Sparse Row Variances
#' @description Computes the variances of rows in a sparse matrix efficiently.
#' @param m A sparse matrix of class \code{dgCMatrix}.
#' @return A numeric vector of row variances.
#' @export
#' @keywords internal
sparseRowVariances <- function(m) {
  row_means <- Matrix::rowMeans(m)
  row_means_sq <- row_means^2
  row_vars <- Matrix::rowMeans(m^2) - row_means_sq
  return(row_vars)
}

#' @title Group Sums for Sparse Matrices
#' @description Sums the columns of a matrix grouped by a factor, optimized for sparse matrices.
#' @param mat A numeric matrix or sparse matrix.
#' @param groups A factor or character vector indicating group membership for each column of \code{mat}.
#' @param sparse Logical indicating whether the input matrix is sparse. Default is \code{FALSE}.
#' @return A matrix with rows corresponding to the rows of \code{mat} and columns corresponding to the unique levels of \code{groups}.
#' @export
#' @keywords internal
groupSums <- function(mat, groups = NULL, sparse = FALSE) {
  stopifnot(!is.null(groups))
  stopifnot(length(groups) == ncol(mat))
  group_levels <- unique(groups)
  group_matrices <- lapply(group_levels, function(group) {
    cols <- which(groups == group)
    if (sparse) {
      Matrix::rowSums(mat[, cols, drop = FALSE])
    } else {
      rowSums(mat[, cols, drop = FALSE])
    }
  })
  result <- do.call(cbind, group_matrices)
  colnames(result) <- group_levels
  return(result)
}


#' @title Perform SVD for Latent Semantic Indexing
#' @description Computes Singular Value Decomposition (SVD) on a term-frequency matrix for dimensionality reduction.
#' @param sp_mat A term-frequency sparse matrix (genes x cells).
#' @param num_dim The number of dimensions to retain.
#' @param mat_only Logical indicating whether to return only the transformed matrix. Default is TRUE.
#' @return If `mat_only` is TRUE, returns the transformed matrix; otherwise, returns a list containing the matrix and SVD components.
#' @importFrom irlba irlba
#' @export
#' @keywords internal
svd_lsi <- function (sp_mat, num_dim, mat_only = T)
{
  svd <- irlba::irlba(sp_mat, num_dim, num_dim)
  svdDiag <- matrix(0, nrow = num_dim, ncol = num_dim)
  diag(svdDiag) <- svd$d
  matSVD <- t(svdDiag %*% t(svd$v))
  rownames(matSVD) <- colnames(sp_mat)
  colnames(matSVD) <- seq_len(ncol(matSVD))
  if (mat_only) {
    return(matSVD)
  }
  else {
    return(list(matSVD = matSVD, svd = svd))
  }
}




#' @title Find Partitions in a Seurat Object Similar to Monocle3
#'
#' @description
#' Finds partitions (low-resolution clusters) in a Seurat object using Louvain or Leiden clustering, similar to the \code{find_partition} function in Monocle3. This function computes partitions by clustering the cells and then grouping clusters into partitions based on connectivity.
#'
#' @param obj A Seurat object.
#' @param method A character string specifying the clustering algorithm to use: \code{"louvain"} or \code{"leiden"}. Default is \code{"louvain"}.
#' @param k Integer specifying the number of nearest neighbors (k) to use when constructing the k-NN graph. Default is 20.
#' @param reduction A character string specifying the dimensionality reduction to use (e.g., \code{"pca"}, \code{"umap"}). Default is \code{"umap"}.
#' @param dims A vector of integers specifying the dimensions to use from the specified reduction. If \code{NULL}, all dimensions are used. Default is \code{NULL}.
#' @param weight Logical indicating whether to use edge weights in clustering. Default is \code{FALSE}.
#' @param num_iter Integer specifying the number of iterations for the Leiden algorithm. Default is 1.
#' @param resolution_parameter Numeric value specifying the resolution parameter for clustering. Default is \code{NULL}, which uses the Seurat default.
#' @param random_seed Integer specifying the random seed for reproducibility. Default is 2020.
#' @param verbose Logical indicating whether to display messages during computation. Default is \code{TRUE}.
#' @param partition_q_value Not used in this implementation. Included for compatibility.
#'
#' @return A Seurat object with partitions added to the metadata (\code{obj$partitions}).
#'
#' @details
#' This function performs clustering using the specified method and then computes partitions by constructing a cluster graph and finding connected components. Cells belonging to clusters within the same connected component are assigned to the same partition.
#'
#' @importFrom Seurat FindNeighbors FindClusters Embeddings Reductions Idents DefaultAssay
#' @importFrom igraph as.igraph set_vertex_attr as_data_frame graph_from_data_frame components V
#' @export
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' library(Seurat)
#' # Assuming 'seurat_obj' is a Seurat object with PCA and UMAP computed
#' seurat_obj <- find_partitions(seurat_obj, method = "leiden", k = 20, reduction = "umap", dims = 1:10)
#' # The partitions can be accessed via seurat_obj$partitions
#' }
find_partitions <- function(obj, method = "louvain", k = 20, reduction = "umap", dims = NULL, weight = FALSE,
                            num_iter = 1, resolution_parameter = NULL, random_seed = 2020,
                            verbose = TRUE, partition_q_value = 0.05) {
  # Ensure necessary packages are available
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Seurat package is required but not installed.")
  }
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("igraph package is required but not installed.")
  }

  # Check that obj is a Seurat object
  if (!inherits(obj, "Seurat")) {
    stop("Input object must be a Seurat object.")
  }

  # Check that reduction exists
  if (!reduction %in% Seurat::Reductions(obj)) {
    stop(paste0("Reduction '", reduction, "' not found in the Seurat object."))
  }

  # Determine dimensions to use
  if (is.null(dims)) {
    dims <- 1:ncol(Seurat::Embeddings(obj, reduction))
  }

  # Compute the nearest neighbor graph
  obj <- Seurat::FindNeighbors(obj, reduction = reduction, dims = dims, k.param = k, verbose = verbose)

  # Perform clustering
  if (method == "leiden") {
    obj <- Seurat::FindClusters(
      obj,
      algorithm = 4,
      resolution = ifelse(is.null(resolution_parameter), 0.8, resolution_parameter),
      random.seed = random_seed,
      verbose = verbose
    )
  } else if (method == "louvain") {
    obj <- Seurat::FindClusters(
      obj,
      algorithm = 1,
      resolution = ifelse(is.null(resolution_parameter), 0.8, resolution_parameter),
      random.seed = random_seed,
      verbose = verbose
    )
  } else {
    stop("Method must be either 'leiden' or 'louvain'.")
  }

  # Get the shared nearest neighbor (SNN) graph
  graph_name <- paste0(Seurat::DefaultAssay(obj), "_snn")
  snn_graph <- obj@graphs[[graph_name]]

  # Convert the SNN graph to an igraph object
  snn_igraph <- igraph::as.igraph(snn_graph)

  # Get the clusters
  clusters <- as.character(Seurat::Idents(obj))

  # Add cluster membership as a vertex attribute
  snn_igraph <- igraph::set_vertex_attr(snn_igraph, "cluster", value = clusters)

  # Build cluster graph
  # For each edge, get the clusters of the two nodes
  edges <- igraph::as_data_frame(snn_igraph, what = "edges")
  vertices <- igraph::as_data_frame(snn_igraph, what = "vertices")

  edges$cluster_from <- vertices$cluster[match(edges$from, vertices$name)]
  edges$cluster_to <- vertices$cluster[match(edges$to, vertices$name)]

  # Keep only edges between different clusters
  cluster_edges <- edges[edges$cluster_from != edges$cluster_to, ]

  # Create an edge list between clusters
  cluster_edge_list <- unique(cluster_edges[, c("cluster_from", "cluster_to")])

  # Build cluster graph
  unique_clusters <- unique(c(cluster_edge_list$cluster_from, cluster_edge_list$cluster_to))
  cluster_vertices <- data.frame(name = unique_clusters, stringsAsFactors = FALSE)
  cluster_graph <- igraph::graph_from_data_frame(cluster_edge_list, directed = FALSE, vertices = cluster_vertices)

  # Compute connected components
  components <- igraph::components(cluster_graph)
  cluster_partitions <- components$membership
  names(cluster_partitions) <- igraph::V(cluster_graph)$name

  # Map cells to partitions
  cell_partitions <- cluster_partitions[clusters]
  partitions <- as.factor(cell_partitions)

  # Add partitions to metadata
  obj$partitions <- partitions

  message(paste0("Found ", length(unique(partitions)), " partitions!"))

  return(obj)
}

