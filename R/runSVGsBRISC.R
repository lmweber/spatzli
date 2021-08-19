#' runSVGsBRISC
#' 
#' Run method to identify spatially variable genes (SVGs) using BRISC
#' 
#' Identify top SVGs using BRISC ("bootstrap for rapid inference on spatial
#' covariances") methodology developed by Saha and Datta (2018).
#' 
#' This function runs BRISC separately for each gene, using parallelization for
#' faster runtime using one core per BRISC run. The main outputs of interest are
#' the covariance parameter estimates stored in 'Theta' in the BRISC output. We
#' use these estimates to calculate summary values 'ratio_sv' defined as
#' 'sigma.sq / tau.sq' (ratio of spatial to non-spatial variance), and 'prop_sv'
#' defined as 'sigma.sq / (sigma.sq + tau.sq)' (proportion of spatial variance
#' out of total variance), which can be used to rank SVGs.
#' 
#' The current version does not run the BRISC bootstrap inference step on the
#' parameter estimates, since this is much slower.
#' 
#' Assumes the input object is a \code{SpatialExperiment} containing an assay
#' named \code{logcounts} and filtered to exclude low-expressed genes, e.g. as
#' prepared with \code{\link{preprocessSVGs}}.
#' 
#' 
#' @param spe \code{SpatialExperiment} Input object, assumed to be a
#'   \code{SpatialExperiment} containing an assay named \code{logcounts} and
#'   spatial coordinates accessible with \code{spatialCoords()}.
#' 
#' @param x \code{numeric matrix} Matrix of covariates, with number of rows
#'   (spots) matching the number of columns (spots) in \code{spe}. Default =
#'   NULL, which is an intercept-only model. See \code{BRISC} documentation for
#'   more details.
#' 
#' @param n_threads \code{integer} Number of threads for parallelization.
#'   Default = 1.
#' 
#' 
#' @return Returns summary statistics and SVG ranks as new columns in
#'   \code{rowData} in \code{spe} object.
#' 
#' 
#' @importFrom SpatialExperiment spatialCoords
#' @importFrom SingleCellExperiment logcounts
#' @importFrom SummarizedExperiment assayNames rowData 'rowData<-'
#' @importFrom BRISC BRISC_estimation
#' @importFrom BiocParallel bplapply MulticoreParam
#' @importFrom Matrix rowMeans
#' 
#' @export
#' 
#' @examples
#' # to do
#' 
runSVGsBRISC <- function(spe, x = NULL, n_threads = 1, ...) {
  
  stopifnot("logcounts" %in% assayNames(spe))
  
  if (!is.null(x)) stopifnot(nrow(x) == ncol(spe))
  
  y <- logcounts(spe)
  
  # ---------
  # run BRISC
  # ---------
  
  # scale coordinates proportionally
  coords <- spatialCoords(spe)
  range_all <- max(apply(coords, 2, function(col) diff(range(col))))
  coords <- apply(coords, 2, function(col) (col - min(col)) / range_all)
  
  # run BRISC using parallelization
  ix <- seq_len(nrow(y))
  out_brisc <- bplapply(ix, function(i) {
    # fit model (note: default if x is NULL is intercept-only model)
    y_i <- y[i, ]
    out_i <- BRISC_estimation(coords = coords, y = y_i, x = x, 
                              n.neighbors = 15, order = "AMMD", 
                              cov.model = "exponential", search.type = "cb", 
                              verbose = FALSE)
    res_i <- c(out_i$Theta, out_i$Beta, runtime = out_i$estimation.time[["elapsed"]])
    res_i
  }, BPPARAM = MulticoreParam(workers = n_threads))
  
  # collapse output list into matrix
  mat_brisc <- do.call("rbind", out_brisc)
  
  # ------------------------------
  # calculate statistics and ranks
  # ------------------------------
  
  # mean logcounts
  mat_brisc <- cbind(
    mat_brisc, 
    mean = rowMeans(y)
  )
  
  # spatial coefficient of variation
  mat_brisc <- cbind(
    mat_brisc, 
    spcov = sqrt(mat_brisc[, "sigma.sq"]) / mat_brisc[, "mean"]
  )
  
  # ratio of spatial to non-spatial variance
  mat_brisc <- cbind(
    mat_brisc, 
    ratio_sv = mat_brisc[, "sigma.sq"] / mat_brisc[, "tau.sq"]
  )
  
  # proportion of spatial variance (out of total variance)
  mat_brisc <- cbind(
    mat_brisc, 
    prop_sv = mat_brisc[, "sigma.sq"] / (mat_brisc[, "sigma.sq"] + mat_brisc[, "tau.sq"])
  )
  
  # return in rowData of spe object
  
  stopifnot(nrow(spe) == nrow(mat_brisc))
  
  rowData(spe) <- cbind(rowData(spe), mat_brisc)
  
  spe
}

