#' @title Plot Grouped Geneset Scores
#'
#' @description
#' Plots geneset scores for each single cell, grouped by a specified variable. The plot includes a violin plot, a boxplot, and optionally jittered points of the raw data.
#'
#' Geneset scores are calculated per cell based on a set of genes. When using method \code{"totals"}, the sum of the size-factor corrected, log-normalized gene expression for the given set of genes is computed. When using method \code{"corrected"}, the geneset scores are corrected using 100 random genes with similar expression levels.
#'
#' @param cds A \code{cell_data_set} object from the \code{monocle3} package.
#' @param marker_set A character vector of gene identifiers corresponding to the genes in the geneset.
#' @param name A character string representing the name of the geneset (used as a label).
#' @param by A character string specifying the cell metadata column in \code{colData(cds)} to group cells by in the plot.
#' @param fData_col A character string specifying the feature data column in \code{rowData(cds)} that contains the gene identifiers in \code{marker_set}. Default is \code{"gene_short_name"}.
#' @param scale A character string specifying how violins should be scaled, passed to \code{\link[ggplot2]{geom_violin}}. Default is \code{"width"}.
#' @param facet A character string specifying an optional cell metadata column in \code{colData(cds)} to facet the plot by. Default is \code{NULL}.
#' @param adjust Numeric value for bandwidth adjustment in \code{geom_violin}. Default is \code{1.4}.
#' @param size Numeric value for point size in \code{geom_jitter}. Default is \code{0.05}.
#' @param alpha Numeric value for point transparency in \code{geom_jitter}. Default is \code{0.1}.
#' @param method A character string specifying the method to calculate geneset scores. Either \code{"totals"} or \code{"corrected"}. Default is \code{"totals"}.
#' @param overlay_violinandbox Logical indicating whether to overlay violin and box plots. Default is \code{TRUE}.
#' @param box_width Numeric value specifying the width of the boxplot. Default is \code{0.3}.
#' @param rotate_x Logical indicating whether to rotate x-axis labels by 90 degrees. Default is \code{TRUE}.
#' @param jitter Logical indicating whether to add jittered points to the plot. Default is \code{TRUE}.
#' @param return_values Logical indicating whether to return the plot and the scores data frame. If \code{FALSE}, only the plot is returned. Default is \code{FALSE}.
#'
#' @return If \code{return_values} is \code{TRUE}, a list with elements:
#' \describe{
#'   \item{\code{plot}}{The ggplot object.}
#'   \item{\code{scores}}{A data frame containing the geneset scores and grouping variables.}
#' }
#' If \code{return_values} is \code{FALSE}, only the ggplot object is returned.
#'
#' @import ggplot2
#' @importFrom monocle3 colData
#' @importFrom methods is
#' @export
#' @keywords internal
#' @examples
#' # Assuming `cds` is a cell_data_set object, and `genes` is a vector of gene identifiers
#' # plot_grouped_geneset(cds, marker_set = genes, name = "MyGeneset", by = "cell_type")
#'
#' @references
#' Puram, S. V. et al. (2017). Single-Cell Transcriptomic Analysis of Primary and Metastatic Tumor Ecosystems in Head and Neck Cancer. \emph{Cell}, 171(7), 1611–1624.e24.
plot_grouped_geneset <- function(cds,
                                 marker_set,
                                 name,
                                 by,
                                 fData_col = "gene_short_name",
                                 scale = "width",
                                 facet = NULL,
                                 adjust = 1.4,
                                 size = 0.05,
                                 alpha = 0.1,
                                 method = c("totals", "corrected"),
                                 overlay_violinandbox = TRUE,
                                 box_width = 0.3,
                                 rotate_x = TRUE,
                                 jitter = TRUE,
                                 return_values = FALSE) {
  method <- match.arg(method)

  if (!methods::is(cds, "cell_data_set")) {
    stop("cds must be a cell_data_set object from the monocle3 package.")
  }

  if (!(by %in% colnames(colData(cds)))) {
    stop("The 'by' parameter must be a column name in colData(cds).")
  }

  if (!is.null(facet) && !(facet %in% colnames(colData(cds)))) {
    stop("The 'facet' parameter must be a column name in colData(cds).")
  }

  if (method == "totals") {
    colData(cds)[[name]] <- estimate_score(cds, marker_set, fData_col = fData_col)
  } else if (method == "corrected") {
    colData(cds)[[name]] <- estimate_corrected_score(cds, marker_set, fData_col = fData_col)
  }

  scores <- data.frame(
    geneset_score = colData(cds)[[name]],
    group = as.factor(colData(cds)[[by]])
  )

  if (!is.null(facet)) {
    scores$facet_var <- colData(cds)[[facet]]
  }

  g <- ggplot(scores, aes(x = group, y = geneset_score, fill = group))

  if (jitter) {
    g <- g + geom_jitter(size = size, alpha = alpha)
  }

  if (overlay_violinandbox) {
    g <- g +
      geom_violin(scale = scale, adjust = adjust) +
      geom_boxplot(width = box_width, fill = "white", outlier.size = 0)
  }

  if (!is.null(facet)) {
    g <- g + facet_wrap(~facet_var, scales = "free")
  }

  if (rotate_x) {
    g <- g + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  }

  if (return_values) {
    return(list(plot = g, scores = scores))
  } else {
    return(g)
  }
}


#' @title Plot Geneset Scores on Cell Embeddings
#'
#' @description
#' Plots geneset scores for each single cell on a dimensionality reduction embedding (e.g., UMAP or t-SNE). The geneset scores are calculated per cell based on a set of genes.
#'
#' When using method \code{"totals"}, the sum of the size-factor corrected, log-normalized gene expression for the given set of genes is computed. When using method \code{"corrected"}, the geneset scores are corrected using 100 random genes with similar expression levels.
#'
#' @param cds A \code{cell_data_set} object from the \code{monocle3} package.
#' @param marker_set A character vector of gene identifiers corresponding to the genes in the geneset.
#' @param name A character string representing the name of the geneset (used as a label).
#' @param fData_col A character string specifying the feature data column in \code{rowData(cds)} that contains the gene identifiers in \code{marker_set}. Default is \code{"gene_short_name"}.
#' @param method A character string specifying the method to calculate geneset scores. Either \code{"totals"} or \code{"corrected"}. Default is \code{"totals"}.
#' @param reduction_method A character string specifying the dimensionality reduction method to use for plotting. Default is \code{"UMAP"}.
#' @param cell_size Numeric value specifying the size of the points in the plot. Default is \code{0.5}.
#'
#' @return A ggplot object representing the cells colored by geneset scores.
#'
#' @importFrom monocle3 colData plot_cells
#' @importFrom methods is
#' @export
#' @keywords internal
#' @examples
#' # Assuming `cds` is a cell_data_set object, and `genes` is a vector of gene identifiers
#' # plot_geneset(cds, marker_set = genes, name = "MyGeneset")
#'
#' @references
#' Puram, S. V. et al. (2017). Single-Cell Transcriptomic Analysis of Primary and Metastatic Tumor Ecosystems in Head and Neck Cancer. \emph{Cell}, 171(7), 1611–1624.e24.
plot_geneset <- function(cds,
                         marker_set,
                         name,
                         fData_col = "gene_short_name",
                         method = c("totals", "corrected"),
                         reduction_method = "UMAP",
                         cell_size = 0.5) {
  method <- match.arg(method)

  if (!methods::is(cds, "cell_data_set")) {
    stop("cds must be a cell_data_set object from the monocle3 package.")
  }

  if (method == "totals") {
    colData(cds)[[name]] <- estimate_score(cds, marker_set, fData_col = fData_col)
  } else if (method == "corrected") {
    colData(cds)[[name]] <- estimate_corrected_score(cds, marker_set, fData_col = fData_col)
  }

  fontsize <- ifelse(nchar(name) > 50, 10, 14)
  loca <- ifelse(method == "totals", "log(sums)", "log(corr.)")

  plot_cells(
    cds,
    color_cells_by = name,
    label_cell_groups = FALSE,
    cell_size = cell_size,
    reduction_method = reduction_method
  ) +
    theme(legend.position = "top") +
    ggtitle(name) +
    theme(
      plot.title = element_text(size = fontsize, face = "bold"),
      legend.text = element_text(size = 9, angle = 90, vjust = 0.5, hjust = 0.3)
    ) +
    labs(color = loca) +
    scale_color_gradientn(colors = c("darkblue", "skyblue", "white", "red", "darkred"))
}
