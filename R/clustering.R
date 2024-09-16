#' @title Iterative Latent Semantic Indexing (LSI) for Single-Cell Data
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
iterative_LSI <- function(object,
                          num_dim = 25,
                          starting_features = NULL,
                          resolution = c(1e-4, 3e-4, 5e-4),
                          num_features = 3000,
                          exclude_features = NULL,
                          do_tf_idf = TRUE,
                          binarize = FALSE,
                          scale = TRUE,
                          log_transform = TRUE,
                          scale_to = 10000,
                          leiden_k = 20,
                          leiden_weight = FALSE,
                          leiden_iter = 1,
                          random_seed = 2020,
                          verbose = FALSE,
                          return_object = TRUE,
                          ...) {
  # Set random seed for reproducibility
  set.seed(random_seed)

  # Determine the class of the input object
  if (methods::is(object, "cell_data_set")) {
    # Monocle3 cell_data_set
    object_type <- "monocle"
  } else if (methods::is(object, "Seurat")) {
    # Seurat object
    object_type <- "seurat"
  } else {
    stop("The input object must be either a Monocle3 cell_data_set or a Seurat object.")
  }

  # Ensure num_features is a vector matching the length of resolution
  if (length(num_features) == 1) {
    num_features <- rep(num_features, length(resolution))
  } else if (length(num_features) != length(resolution)) {
    stop("Length of num_features must be 1 or match the length of resolution.")
  }

  # Extract expression matrix based on object type
  if (object_type == "monocle") {
    mat <- SummarizedExperiment::assay(object)
  } else if (object_type == "seurat") {
    mat <- Seurat::GetAssayData(object, assay = DefaultAssay(object), slot = "counts")
  }

  # Exclude specified features
  if (!is.null(exclude_features)) {
    mat <- mat[!rownames(mat) %in% exclude_features, ]
  }

  # Binarize data if specified
  if (binarize) {
    if (verbose) {
      message("Binarizing data...")
    }
    mat@x[mat@x > 0] <- 1
  }

  # Initialize variables
  original_features <- rownames(mat)
  iterations_list <- list()

  # Scale data if specified
  if (scale) {
    mat_norm <- t(t(mat) / Matrix::colSums(mat)) * scale_to
  } else {
    mat_norm <- mat
  }

  # Log-transform data if specified
  if (log_transform) {
    mat_norm@x <- log2(mat_norm@x + 1)
  }

  # Start iterations
  for (iteration in seq_along(resolution)) {
    if (verbose) {
      message("Performing iteration ", iteration, "...")
    }

    # Select features
    if (iteration == 1 && !is.null(starting_features)) {
      if (!all(starting_features %in% rownames(mat_norm))) {
        stop("Not all starting_features found in the data.")
      }
      feature_indices <- which(rownames(mat_norm) %in% starting_features)
    } else if (iteration == 1) {
      feature_variances <- sparseRowVariances(mat_norm)
      feature_indices <- head(order(feature_variances, decreasing = TRUE), num_features[iteration])
    } else {
      feature_variances <- matrixStats::rowVars(cluster_mat)
      feature_indices <- head(order(feature_variances, decreasing = TRUE), num_features[iteration])
    }

    # Subset matrix to selected features
    mat_subset <- mat_norm[feature_indices, ]

    # Perform TF-IDF transformation if specified
    if (do_tf_idf) {
      tfidf_mat <- tf_idf_transform(mat_subset)
      tfidf_mat@x[is.na(tfidf_mat@x)] <- 0
    } else {
      tfidf_mat <- mat_subset
    }

    # Perform SVD using irlba
    svd_result <- irlba::irlba(tfidf_mat, nv = num_dim)
    lsi_embeddings <- svd_result$u %*% diag(svd_result$d)
    rownames(lsi_embeddings) <- colnames(mat)
    colnames(lsi_embeddings) <- paste0("LSI_", 1:ncol(lsi_embeddings))

    # Clustering based on object type
    if (object_type == "monocle") {
      # Store LSI embeddings in reducedDims
      SingleCellExperiment::reducedDims(object)[["LSI"]] <- lsi_embeddings

      # Perform clustering
      object <- monocle3::cluster_cells(
        object,
        reduction_method = "LSI",
        k = leiden_k,
        weight = leiden_weight,
        num_iter = leiden_iter,
        resolution = resolution[iteration],
        random_seed = random_seed,
        verbose = verbose,
        ...
      )

      clusters <- monocle3::clusters(object)

      # Compute cluster matrix
      cluster_mat <- groupSums(mat, groups = clusters, sparse = TRUE)
      cluster_mat <- edgeR::cpm(cluster_mat, log = TRUE, prior.count = 3)
    } else if (object_type == "seurat") {
      # Store LSI embeddings in Seurat object
      object[["lsi"]] <- Seurat::CreateDimReducObject(embeddings = lsi_embeddings, key = "LSI_", assay = DefaultAssay(object))

      # Perform clustering
      object <- Seurat::FindNeighbors(
        object,
        reduction = "lsi",
        dims = 1:num_dim,
        k.param = leiden_k,
        verbose = verbose
      )
      object <- Seurat::FindClusters(
        object,
        resolution = resolution[iteration],
        algorithm = 4,  # Algorithm 4 is Leiden in Seurat
        random.seed = random_seed,
        verbose = verbose
      )

      clusters <- Seurat::Idents(object)

      # Compute cluster matrix
      cluster_mat <- groupSums(mat, groups = clusters, sparse = TRUE)
      cluster_mat <- edgeR::cpm(cluster_mat, log = TRUE, prior.count = 3)
    }

    # Store intermediate results if requested
    iteration_name <- paste0("iteration_", iteration)
    iterations_list[[iteration_name]] <- list(
      lsi_embeddings = lsi_embeddings,
      features = original_features[feature_indices],
      clusters = clusters
    )

    # Update mat_norm for next iteration (if needed)
    # In this implementation, mat_norm remains the same unless custom updates are needed
  }

  if (return_object) {
    return(object)
  } else {
    return(list(
      lsi_embeddings = lsi_embeddings,
      clusters = clusters,
      iterations = iterations_list
    ))
  }
}

# Helper Functions

#' @title TF-IDF Transformation
#' @description Performs TF-IDF transformation on a matrix.
#' @param mat A sparse matrix.
#' @return A sparse matrix after TF-IDF transformation.
#' @export
tf_idf_transform <- function(mat) {
  # Compute term frequency (TF)
  tf <- Matrix::t(Matrix::t(mat) / Matrix::colSums(mat))
  # Compute inverse document frequency (IDF)
  idf <- log(1 + ncol(mat) / Matrix::rowSums(mat))
  # Compute TF-IDF
  tf_idf <- Matrix::Diagonal(x = as.vector(idf)) %*% tf
  return(tf_idf)
}

#' @title Sparse Row Variances
#' @description Computes the variances of rows in a sparse matrix efficiently.
#' @param m A sparse matrix of class \code{dgCMatrix}.
#' @return A numeric vector of row variances.
#' @export
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

