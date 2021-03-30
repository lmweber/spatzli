#' rankSVGsBRISC
#' 
#' Rank spatially variable genes (SVGs) using BRISC approach
#' 
#' Calculate ranking of spatially variable genes (SVGs) using the BRISC approach
#' (bootstrap for rapid inference on spatial covariances) developed by Saha and
#' Datta (2018).
#' 
#' 
#' @param spe Input object (SpatialExperiment). Assumed to contain "logcounts"
#'   assay and spatial coordinates retrievable as "spatialCoords".
#' 
#' @param x Matrix of covariates for "BRISC_estimation()". Number of rows must
#'   equal number of spots. See "?BRISC_estimation" for details. Default = NULL
#'   (intercept-only model).
#' 
#' @param n.neighbors Number of nearest neighbors for BRISC. Default = 15.
#' 
#' @param filter_counts Filter out zero-expressed and low-expressed genes by
#'   keeping genes with at least this number of UMI counts in at least one spot.
#'   Default = 5.
#' 
#' @param filter_mito Whether to filter out mitochondrial genes. (Mitochondrial
#'   genes are very highly expressed and may rank highly as SVGs if included,
#'   but are often not of primary biological interest.) Default = TRUE.
#' 
#' @param gene_name Column in "rowData" containing gene names if "filter_mito =
#'   TRUE". Default = "gene_name".
#' 
#' @param n_threads Number of threads for parallelization. Default = 4.
#' 
#' @param ... Additional arguments to pass to "BRISC_estimation()".
#' 
#' 
#' @return Returns test statistics and rankings of detected SVGs as new columns
#'   in "rowData". Genes that were not tested due to filtering (zero expression,
#'   low expression, and/or mitochondrial genes) are given NA values.
#' 
#' 
#' @importFrom SpatialExperiment spatialCoords
#' @importFrom SingleCellExperiment logcounts counts
#' @importFrom SummarizedExperiment rowData 'rowData<-'
#' @importFrom BiocParallel bplapply MulticoreParam
#' @importFrom BRISC BRISC_estimation
#' 
#' @export
#' 
#' @examples
#' # to do
#' 
rankSVGsBRISC <- function(spe, x = NULL, n.neighbors = 15, 
                          filter_counts = 5, filter_mito = TRUE, gene_name = "gene_name", 
                          n_threads = 4, ...) {
  
  if (!("logcounts" %in% assayNames(spe))) stop("input object must contain 'logcounts' assay")
  if (!("counts" %in% assayNames(spe))) stop("input object must contain 'counts' assay")
  
  if (!is.null(x)) stopifnot(nrow(x) == ncol(spe))
  
  # ---------
  # filtering
  # ---------
  
  # filtering: identify low-expressed genes
  is_low <- !apply(counts(spe), 1, function(row) any(row >= filter_counts))
  # filtering: identify mitochondrial genes
  if (filter_mito) {
    is_mito <- grepl("(^MT-)|(^mt-)", rowData(spe)[, gene_name])
  }
  # filtering: combined set of genes to discard
  if (filter_mito) {
    discard <- is_low | is_mito
  } else {
    discard <- is_low
  }
  # indices of genes to keep
  stopifnot(length(discard) == nrow(spe))
  ix_keep <- seq_len(nrow(spe))[!discard]
  
  # ---------------------------------------------------------
  # fit models and calculate statistics and rankings per gene
  # ---------------------------------------------------------
  
  y <- logcounts(spe)
  
  coords <- spatialCoords(spe)
  # scale coordinates
  # to do: scale proportionally
  coords <- apply(coords, 2, function(col) (col - min(col)) / (max(col) - min(col)))
  
  # parallelized
  out_brisc <- bplapply(ix_keep, function(i) {
    # fit intercept-only model for each gene, using mostly default parameters
    out_i <- BRISC_estimation(coords = coords, y = y[i, ], x = x, n.neighbors = n.neighbors, 
                              n_omp = 1, verbose = FALSE, ...)
    # return estimated parameters and runtime
    c(out_i$Theta, runtime = out_i$estimation.time[["elapsed"]])
  }, BPPARAM = MulticoreParam(workers = n_threads))
  
  # collapse list
  stopifnot(length(out_brisc) == length(ix_keep))
  mat_brisc <- do.call("rbind", out_brisc)
  stopifnot(nrow(mat_brisc) == length(ix_keep))
  
  # calculate fraction spatial variance (FSV)
  mat_brisc <- cbind(
    mat_brisc, 
    fsv = mat_brisc[, "sigma.sq"] / (mat_brisc[, "sigma.sq"] + mat_brisc[, "tau.sq"])
  )
  
  # calculate reversed ranks on sigma.sq
  rank_sigma.sq <- rank(-1 * mat_brisc[, "sigma.sq"])
  mat_brisc <- cbind(mat_brisc, rank_sigma.sq = rank_sigma.sq)
  
  # match to correct rows and store in rowData
  mat_brisc_all <- matrix(NA, nrow = nrow(spe), ncol = ncol(mat_brisc))
  colnames(mat_brisc_all) <- colnames(mat_brisc)
  mat_brisc_all[ix_keep, ] <- mat_brisc
  
  rowData(spe) <- cbind(rowData(spe), mat_brisc_all)
  
  spe
}

