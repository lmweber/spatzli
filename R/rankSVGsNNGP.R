#' rankSVGsNNGP
#' 
#' Rank spatially variable genes (SVGs) using spNNGP approach
#' 
#' Calculate ranking of spatially variable genes (SVGs) using the spNNGP
#' approach (spatial nearest neighbors Gaussian processes) developed by Finley
#' et al. (2020), Finley et al. (2019), and Datta et al. (2016).
#' 
#' 
#' @param spe Input object (SpatialExperiment). Assumed to contain "logcounts"
#'   assay and spatial coordinates retrievable as "spatialCoords".
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
#' 
#' @return Returns test statistics and rankings of detected SVGs as new columns
#'   in "rowData". Genes that were not tested due to filtering (zero expression,
#'   low expression, and/or mitochondrial genes) are given NA values.
#' 
#' 
#' @importFrom SpatialExperiment spatialCoords
#' @importFrom SingleCellExperiment logcounts counts
#' @importFrom SummarizedExperiment rowData
#' @importFrom BiocParallel bplapply MulticoreParam
#' @importFrom spNNGP spNNGP
#' @importFrom matrixStats rowMedians
#' 
#' @export
#' 
#' @examples
#' # to do
#' 
rankSVGsNNGP <- function(spe, filter_counts = 5, 
                         filter_mito = TRUE, gene_name = "gene_name", 
                         n_threads = 4) {
  
  if (!("logcounts" %in% assayNames(spe))) stop("input object must contain 'logcounts' assay")
  if (!("counts" %in% assayNames(spe))) stop("input object must contain 'counts' assay")
  
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
  
  y <- logcounts(spe_sub)
  
  coords <- spatialCoords(spe)
  # scale coordinates
  coords <- apply(coords, 2, function(col) (col - min(col)) / (max(col) - min(col)))
  
  # default parameters
  n.samples <- 2000
  starting <- list("phi" = 3/0.5, "sigma.sq" = 50, "tau.sq" = 1)
  tuning <- list("phi" = 0.05, "sigma.sq" = 0.05, "tau.sq" = 0.05)
  priors <- list("beta.Norm" = list(rep(0, 1), diag(1000, 1)), 
                 "phi.Unif" = c(3/1, 3/0.1), "sigma.sq.IG" = c(2, 2), 
                 "tau.sq.IG" = c(2, 0.1))
  method <- "latent"
  n.neighbors <- 5
  cov.model <- "exponential"
  return.neighbor.info <- FALSE
  n.omp.threads <- 1
  
  # parallelized
  out_spnngp <- bplapply(ix_keep, function(i) {
    # fit intercept-only model for each gene
    out_i <- spNNGP(y[i, ] ~ 1, coords = coords, starting = starting, method = method, 
                    n.neighbors = n.neighbors, tuning = tuning, priors = priors, 
                    cov.model = cov.model, n.samples = n.samples, 
                    return.neighbor.info = return.neighbor.info, n.omp.threads = n.omp.threads)
    # return sum of absolute values of medians of posterior samples of spatial random effects
    list(
      stat = sum(abs(rowMedians(out_i$p.w.samples))), 
      runtime = out_i$run.time
    )
  }, BPPARAM = MulticoreParam(workers = n_threads))
  
  # collapse list
  mat_spnngp <- do.call("rbind", out_spnngp)
  mat_spnngp <- apply(mat_spnngp, 2, as.numeric)
  stopifnot(nrow(mat_spnngp) == length(ix_keep))
  
  # calculate reversed ranks
  rank_spnngp <- rank(-1 * mat_spnngp[, "stat"])
  mat_spnngp <- cbind(mat_spnngp, rank = rank_spnngp)
  
  # match to correct rows and store in rowData
  mat_spnngp_all <- matrix(NA, nrow = nrow(spe), ncol = ncol(mat_spnngp))
  colnames(mat_spnngp_all) <- colnames(mat_spnngp)
  mat_spnngp_all[ix_keep, ] <- mat_spnngp
  
  rowData(spe) <- cbind(rowData(spe), mat_spnngp_all)
  
  spe
}

