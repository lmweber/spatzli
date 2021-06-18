#' rankSVGsBRISC
#' 
#' Rank spatially variable genes (SVGs) using BRISC approach
#' 
#' Calculate ranking of spatially variable genes (SVGs) using the BRISC approach
#' (bootstrap for rapid inference on spatial covariances) developed by Saha and
#' Datta (2018).
#' 
#' 
#' @param spe Input object (SpatialExperiment). Assumed to contain assays named
#'   logcounts", and spatial coordinates accessible with "spatialCoords()".
#' 
#' @param x Matrix of covariates for "BRISC_estimation()". Number of rows must
#'   equal number of spots. See "?BRISC_estimation" for details. Default = NULL
#'   (intercept-only model).
#' 
#' @param n_threads Number of threads for parallelization. Default = 4.
#' 
#' @param ... Additional arguments to pass to "BRISC_estimation()".
#' 
#' 
#' @return Returns summary statistics from BRISC as new columns in "rowData",
#'   which can be used to rank SVGs.
#' 
#' 
#' @importFrom SpatialExperiment spatialCoords
#' @importFrom SingleCellExperiment logcounts counts
#' @importFrom SummarizedExperiment rowData 'rowData<-'
#' @importFrom Matrix rowMeans
#' @importFrom BiocParallel bplapply MulticoreParam
#' @importFrom BRISC BRISC_estimation
#' 
#' @export
#' 
#' @examples
#' # to do
#' 
rankSVGsBRISC <- function(spe, x = NULL, n_threads = 4, ...) {
  
  if (!("logcounts" %in% assayNames(spe))) stop("input object must contain 'logcounts' assay")
  
  if (!is.null(x)) stopifnot(nrow(x) == ncol(spe))
  
  y <- logcounts(spe)
  
  # scale coordinates proportionally
  coords <- spatialCoords(spe)
  range_all <- max(apply(coords, 2, function(col) diff(range(col))))
  coords <- apply(coords, 2, function(col) (col - min(col)) / range_all)
  
  # parallelized
  ix <- seq_len(nrow(y))
  out_brisc <- bplapply(ix, function(i) {
    # fit model (default if x is NULL is intercept-only model)
    out_i <- BRISC_estimation(coords = coords, y = y[i, ], x = x, n_omp = 1, verbose = FALSE, ...)
    res_i <- c(out_i$Theta, out_i$Beta, runtime = out_i$estimation.time[["elapsed"]])
    res_i
  }, BPPARAM = MulticoreParam(workers = n_threads))
  
  # collapse list
  mat_brisc <- do.call("rbind", out_brisc)
  
  # include mean logcounts
  mat_brisc <- cbind(
    mat_brisc, 
    mean = rowMeans(y)
  )
  
  # calculate spatial coefficient of variation
  mat_brisc <- cbind(
    mat_brisc, 
    spcov = sqrt(mat_brisc[, "sigma.sq"]) / mat_brisc[, "mean"]
  )
  
  # calculate fraction spatial variance (FSV)
  mat_brisc <- cbind(
    mat_brisc, 
    fsv = mat_brisc[, "sigma.sq"] / (mat_brisc[, "sigma.sq"] + mat_brisc[, "tau.sq"])
  )
  
  # store in rowData
  stopifnot(nrow(rowData(spe)) == nrow(mat_brisc))
  rowData(spe) <- cbind(rowData(spe), mat_brisc)
  
  spe
}

