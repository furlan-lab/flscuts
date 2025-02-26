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

#' Generate a Dot Plot with ADT Data from a Seurat Object
#'
#' This function creates a dot plot using Antibody-Derived Tag (ADT) data from a Seurat object, applying log transformation and optional noise to the data. The function also provides customization options for colors, downsampling, and buffer adjustments for the plot axes.
#'
#' @param seu A Seurat object containing ADT data.
#' @param x A string representing the name of the feature to plot on the x-axis.
#' @param y A string representing the name of the feature to plot on the y-axis.
#' @param color (Optional) A string representing the name of the feature to use for coloring the dots. Defaults to \code{NULL}.
#' @param addvalue A numeric value added to each feature value before log transformation. Defaults to \code{1}.
#' @param downsample A numeric value between 0 and 1 specifying the fraction of the data to randomly sample. Defaults to \code{0.2}.
#' @param cols (Optional) A vector of colors to use for coloring the dots. Can be used for both continuous and categorical color scales. Defaults to \code{NULL}.
#' @param add_noise A logical indicating whether to add noise (via jittering) to the data before plotting. Defaults to \code{TRUE}.
#' @param size A numeric value specifying the size of the dots in the plot. Defaults to \code{0.8}.
#' @param xpbuff A numeric value specifying the positive buffer for the x-axis limits. Defaults to \code{0.5}.
#' @param ypbuff A numeric value specifying the positive buffer for the y-axis limits. Defaults to \code{0.5}.
#' @param xnbuff A numeric value specifying the negative buffer for the x-axis limits. Defaults to \code{0}.
#' @param ynbuff A numeric value specifying the negative buffer for the y-axis limits. Defaults to \code{0}.
#' @param remove_negs A logical value indicating whether to remove points with negative values on the x or y axes after log transformation. Defaults to \code{TRUE}.
#'
#' @return A ggplot object representing the dot plot.
#'
#' @examples
#' # Example usage
#' adt_dotplot(seu = my_seurat_object,
#'             x = "ADT_feature1",
#'             y = "ADT_feature2",
#'             color = "ADT_feature3")
#'
#' @importFrom Seurat FetchData
#' @importFrom ggplot2 ggplot aes geom_point theme_bw xlim ylim scale_color_gradient scale_color_manual xlab ylab guides guide_legend
#' @export
adt_dotplot <- function(seu,
                        x,
                        y,
                        color = NULL,
                        addvalue = 1,
                        downsample = 0.2,
                        cols = NULL,
                        add_noise = TRUE,
                        size = 0.8,
                        xpbuff = 0.5,
                        ypbuff = 0.5,
                        xnbuff = 0,
                        ynbuff = 0,
                        remove_negs = TRUE
                        ){
  vars <- c(x, y, color)
  df<-FetchData(seu, vars=vars, layer = "counts")
  if(add_noise){
    amount<-0.8
    factor<-1
    nbsize = 3
    nbprob = 0.5
    input<-rbinom(length(df[[vars[2]]]), size = nbsize, prob = nbprob)+df[[vars[2]]]
    df$HTO_Y<-log10(jitter(input+addvalue, amount=amount, factor = factor))
    input<-rbinom(length(df[[vars[1]]]), size = nbsize, prob = nbprob)+df[[vars[1]]]
    df$HTO_X<-log10(jitter(input+addvalue, amount=amount, factor = factor))
  } else {
    df$HTO_Y<-log10(df[[vars[2]]]+addvalue)
    df$HTO_X<-log10(df[[vars[2]]]+addvalue)
  }
  if(remove_negs){
    df <- df[!(df$HTO_Y < 0),]
    df <- df[!(df$HTO_X < 0),]
  }
  df$color <- df[[color]]
  ymax = max(df$HTO_Y)+ ypbuff
  xmax = max(df$HTO_X)+ xpbuff
  ymin = min(df$HTO_Y)- ynbuff
  xmin = min(df$HTO_X)- xnbuff
  g <- ggplot(df, aes(x=HTO_X, y=HTO_Y, color=color))+
    geom_point(size = size)+
    xlim(xmin, xmax)+
    ylim(ymin, ymax)+
    theme_bw()+
    xlab(x)+
    ylab(y)+
    guides(colour = guide_legend(override.aes = list(size=3)))

  if(!is.null(cols)){
    if(is.numeric(df$color)){
      g+scale_color_gradient(cols)
    } else {
      g+scale_color_manual(values=cols)
    }

  } else {
    g
  }
}

#' Generate a Panel of Dot Plots for AML Markers
#'
#' This function creates a panel of dot plots for Acute Myeloid Leukemia (AML) marker expression, using ADT data from a Seurat object. Each plot compares pairs of AML-related surface markers, and the plots are arranged in a grid. The color of the points can be customized based on a specified feature.
#'
#' @param seu A Seurat object containing ADT data.
#' @param color A string representing the name of the feature to use for coloring the dots. Defaults to \code{"vmR_pred"}.
#' @param cols A vector of colors to use for the coloring the dots. Can be used for both continuous and categorical color scales.
#'
#' @return A patchwork object representing a grid of dot plots.
#'
#' @examples
#' # Example usage
#' aml_panel(seu = my_seurat_object, color = "vmR_pred", cols = my_colors)
#'
#' @importFrom Seurat FetchData
#' @importFrom ggplot2 NoLegend
#' @importFrom patchwork +
#' @export
aml_panel <- function(seu, color = "vmR_pred", cols){
  p1<-adt_dotplot(seu, "CD45", "CD34", color, cols = cols)+NoLegend()
  p2<-adt_dotplot(seu, "CD117", "CD34", color, cols = cols)+NoLegend()
  p3<-adt_dotplot(seu, "CD38", "HLA-DRA", color, cols = cols)+NoLegend()
  p4<-adt_dotplot(seu, "CD33", "CD13", color, cols = cols)+NoLegend()
  p5<-adt_dotplot(seu, "CD71", "CD15", color, cols = cols)+NoLegend()
  p6<-adt_dotplot(seu, "CD14", "CD64", color, cols = cols)+NoLegend()
  p7<-adt_dotplot(seu, "CD34", "CD123", color, cols = cols)+NoLegend()
  p8<-adt_dotplot(seu, "CD11b", "CD11c", color, cols = cols)+NoLegend()
  p9<-adt_dotplot(seu, "CD56", "CD34", color, cols = cols)+NoLegend()
  p1+p2+p3+p4+p5+p6+p7+p8+p9
}


#' Generate a Panel of Dot Plots for ALL Markers
#'
#' This function creates a panel of dot plots for Acute Myeloid Leukemia (ALL) marker expression, using ADT data from a Seurat object. Each plot compares pairs of AML-related surface markers, and the plots are arranged in a grid. The color of the points can be customized based on a specified feature.
#'
#' @param seu A Seurat object containing ADT data.
#' @param color A string representing the name of the feature to use for coloring the dots. Defaults to \code{"vmR_pred"}.
#' @param cols A vector of colors to use for the coloring the dots. Can be used for both continuous and categorical color scales.
#'
#' @return A patchwork object representing a grid of dot plots.
#'
#' @examples
#' # Example usage
#' aml_panel(seu = my_seurat_object, color = "vmR_pred", cols = my_colors)
#'
#' @importFrom Seurat FetchData
#' @importFrom ggplot2 NoLegend
#' @importFrom patchwork +
#' @export

all_panel <- function(seu, color = "vmR_pred", cols){
  p1<-adt_dotplot(seu, "CD45", "CD10", color, cols = cols)+NoLegend()
  p2<-adt_dotplot(seu, "CD10", "CD19", color, cols = cols)+NoLegend()
  p3<-adt_dotplot(seu, "CD34", "CD22", color, cols = cols)+NoLegend()
  p4<-adt_dotplot(seu, "CD20", "CD38", color, cols = cols)+NoLegend()
  #p5<-adt_dotplot(seu, "CD79a", "CD13", color, cols = cols)+NoLegend()
  #p6<-adt_dotplot(seu, "CD58", "CD9", color, cols = cols)+NoLegend()
  p5<-adt_dotplot(seu, "CD19", "CD33", color, cols = cols)+NoLegend()
  #p8<-adt_dotplot(seu, "Lambda", "Kappa", color, cols = cols)+NoLegend()
  p6<-adt_dotplot(seu, "CD3D", "CD7", color, cols = cols)+NoLegend()
  p7<-adt_dotplot(seu, "CD4", "CD8A", color, cols = cols)+NoLegend()
  p8<-adt_dotplot(seu, "CD56", "CD16", color, cols = cols)+NoLegend()
  #p12<-adt_dotplot(seu, "CD5", "CD1a", color, cols = cols)+NoLegend()
  p9<-adt_dotplot(seu, "CD34", "CD2", color, cols = cols)+NoLegend()
  p10<-adt_dotplot(seu, "TCRab", "TCRVd2", color, cols = cols)+NoLegend()
  #p22<-adt_dotplot(seu, "CD7", "CD1a", color, cols = cols)+NoLegend()
  p11<-adt_dotplot(seu, "CD34", "CD5", color, cols = cols)+NoLegend()
  p1+p2+p3+p4+p5+p6+p7+p8+p9+p10+p11
}
