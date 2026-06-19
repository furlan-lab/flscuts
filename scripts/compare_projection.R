#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Head-to-head comparison: project_data() vs project_data2()
# ---------------------------------------------------------------------------
# project_data()  -> serial `lapply` pseudo-cell simulation (current flscuts)
# project_data2() -> m3addon-faithful parallel `pbmcapply::pbmclapply` simulation
#
# Both now read the SAME flscuts object layout (int_metadata$LSI_model and the
# in-memory reduce_dim_aux[["UMAP"]]$model$umap_model), so they can be run on the
# same projector / projectee objects.
#
# HOW TO USE
#   1. Build / load your `projector` and `projectee` objects (see the EDIT block).
#      The projector must have been run through iterative_LSI(..., run_umap = TRUE).
#   2. Set the parameters in the EDIT block.
#   3. Run:  Rscript scripts/compare_projection.R     (or source() it in RStudio)
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(flscuts)
  library(ggplot2)
})

## ======================= EDIT THIS BLOCK ==================================
# Provide your two objects. Either readRDS/qs_read them here, or assign objects
# that already exist in your session before sourcing this script.
#
# Example (LSI workflow):
#   projector <- readRDS("/path/to/reference_with_LSI_UMAP.RDS")
#   projectee <- readRDS("/path/to/query.RDS")

projector <- if (exists("projector")) projector else stop("Define `projector` (reference with LSI + UMAP model).")
projectee <- if (exists("projectee")) projectee else stop("Define `projectee` (query / bulk to project).")

params <- list(
  reduced_dim               = "LSI",        # "LSI" for monocle3, "lsi" for Seurat
  embedding                 = "UMAP",       # "UMAP" for monocle3, "umap" for Seurat
  features                  = "annotation-based",  # or "range-based" (monocle3 + GRanges)
  make_pseudo_single_cells  = FALSE,        # TRUE for bulk projectee
  n                         = 250,
  ncells_coembedding        = 5000,
  threads                   = 6,
  seed                      = 2020,
  force                     = FALSE
)

# Also run project_data2 with the Annoy-index reload (save->load umap model).
# Set TRUE to test whether a stale in-memory nn_index pointer is hurting the
# projection (the m3addon-faithful behavior).
test_saved_umap_model <- TRUE

# Optional: a metadata column to color the projectee by, present in projectee.
projectee_label_col <- NULL                 # e.g. "celltype" or "ct2"
# Optional: a projector metadata column used for nearest-neighbor label transfer.
projector_label_col <- NULL                 # e.g. "celltype" / "BioClassification"

out_prefix <- "projection_comparison"       # output file prefix (png written to wd)
## ==========================================================================


run_one <- function(fun, label, extra = list()) {
  message("\n==== Running ", label, " ====")
  t0 <- Sys.time()
  res <- do.call(fun, c(
    list(projector = projector, projectee = projectee,
         projectee_label_col = projectee_label_col),
    params, extra
  ))
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  message(sprintf("%s finished in %.1f s", label, elapsed))
  list(res = res, secs = elapsed)
}

runs <- list()
runs[["project_data"]]  <- run_one(project_data,  "project_data (serial lapply)")
runs[["project_data2"]] <- run_one(project_data2, "project_data2 (parallel pbmclapply, in-memory model)")
if (isTRUE(test_saved_umap_model)) {
  runs[["project_data2_reload"]] <- run_one(
    project_data2, "project_data2 (Annoy reload via save/load_umap_model)",
    extra = list(use_saved_umap_model = TRUE)
  )
}

a <- runs[["project_data"]]
b <- runs[["project_data2"]]

## --- Quantitative agreement ------------------------------------------------
# How consistent are the two projected embeddings for the same projectee points?
proj_corr <- function(r1, r2) {
  u1 <- as.data.frame(r1$projectedUMAP)
  u2 <- as.data.frame(r2$projectedUMAP)
  # align by rowname (pseudo-cells: "<sample>#<k>"; single: projectee colname)
  common <- intersect(rownames(u1), rownames(u2))
  if (length(common) < 3) return(c(UMAP1 = NA_real_, UMAP2 = NA_real_))
  c(UMAP1 = cor(u1[common, "UMAP1"], u2[common, "UMAP1"]),
    UMAP2 = cor(u1[common, "UMAP2"], u2[common, "UMAP2"]))
}
cc <- proj_corr(a$res, b$res)
message(sprintf("\nProjected-coordinate correlation between methods: UMAP1 = %.3f, UMAP2 = %.3f",
                cc["UMAP1"], cc["UMAP2"]))

# Optional nearest-neighbor label-transfer concordance per method.
label_transfer_accuracy <- function(res) {
  if (is.null(projector_label_col) || is.null(projectee_label_col)) return(NA_real_)
  ref  <- as.data.frame(res$singleCellUMAP)
  qry  <- as.data.frame(res$projectedUMAP)
  nn <- RANN::nn2(data = ref[, c("UMAP1", "UMAP2")],
                  query = qry[, c("UMAP1", "UMAP2")], k = 1)
  if (methods::is(projector, "Seurat")) {
    ref_lab <- projector@meta.data[[projector_label_col]]
  } else {
    ref_lab <- as.character(SummarizedExperiment::colData(projector)[[projector_label_col]])
  }
  predicted <- ref_lab[nn$nn.idx[, 1]]
  truth <- as.character(qry$projectee_labels)
  mean(predicted == truth, na.rm = TRUE)
}
summary_tbl <- data.frame(
  method       = names(runs),
  seconds      = vapply(runs, function(r) r$secs, numeric(1)),
  nn_label_acc = vapply(runs, function(r) label_transfer_accuracy(r$res), numeric(1)),
  row.names    = NULL
)
message("\n--- Summary ---")
print(summary_tbl, row.names = FALSE)

## --- Side-by-side plots ----------------------------------------------------
plots <- lapply(names(runs), function(nm) {
  plot_projection(runs[[nm]]$res, projector, projectee,
                  projectee_col = projectee_label_col) +
    ggtitle(nm) + theme_void()
})

combined <- cowplot::plot_grid(plotlist = plots, nrow = 1)
png_file <- paste0(out_prefix, ".png")
ggsave(png_file, combined, width = 7 * length(plots), height = 6, dpi = 150)
message("\nWrote side-by-side plot to: ", normalizePath(png_file))

# Objects left in the environment for interactive inspection:
#   runs (named list of results + timings), summary_tbl, plots, combined
invisible(list(runs = runs, summary = summary_tbl))
