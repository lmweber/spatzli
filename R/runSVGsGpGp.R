#' runSVGsGpGp
#' 
#' Run method to identify spatially variable genes (SVGs) using GpGp
#' 
#' Identify top SVGs using GpGp (Guinness 2018).
#' 
#' This function runs GpGp separately for each gene, using parallelization for
#' faster runtime using one core per GpGp run. The main outputs of interest are
#' the covariance parameter estimates stored in 'covparms' in the GpGp output
#' (variance, range, nugget, if using 'exponential_isotropic' covariance
#' function). We fix the 'range' parameter (i.e. phi) so it is equal across
#' genes. In addition, the output 'loglik' can be used for likelihood ratio
#' tests.
#' 
#' Note parameterization used by GpGp: total variance = sigmasq + (tausq *
#' sigmasq); i.e. the nugget is 'tausq * sigmasq', not 'tausq' itself.
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
#'   NULL, which is an intercept-only model. See \code{GpGp} documentation for
#'   more details.
#' 
#' @param fix_param_range \code{numeric} Fixed parameter value to use for
#'   'range' covariance function parameter in 'exponential_isotropic' covariance
#'   function (corresponding to 'phi' in other parameterizations). Default =
#'   0.5. See \code{GpGp} documentation for details on parameterization.
#' 
#' @param n_neighbors \code{numeric} Number of nearest neighbors. See
#'   \code{GpGp} documentation for details.
#' 
#' @param n_threads \code{integer} Number of threads for parallelization.
#'   Default = 1.
#' 
#' @param verbose \code{logical} Whether to display verbose output from
#'   \code{GpGp}. Default = FALSE.
#' 
#' 
#' @return Returns summary statistics and SVG ranks as new columns in
#'   \code{rowData} in \code{spe} object.
#' 
#' 
#' @importFrom SpatialExperiment spatialCoords
#' @importFrom SingleCellExperiment logcounts
#' @importFrom SummarizedExperiment assayNames rowData 'rowData<-'
#' @importFrom GpGp order_middleout find_ordered_nn fit_model
#' @importFrom BiocParallel bplapply MulticoreParam
#' @importFrom Matrix rowMeans
#' @importFrom matrixStats rowVars
#' 
#' @export
#' 
#' @examples
#' library(SpatialExperiment)
#' library(STexampleData)
#' library(spatzli)
#' 
#' spe <- Visium_humanDLPFC()
#' 
#' spe <- preprocessSVGs(spe)
#' 
#' # subset 1 gene
#' spe_1 <- spe[1, ]
#' system.time({
#'   spe_1 <- runSVGsGpGp(spe_1, x = NULL, n_threads = 1, verbose = TRUE)
#' })
#' 
#' # subset 100 genes
#' # spe_100 <- spe[1:100, ]
#' # spe_100 <- runSVGsGpGp(spe_100, x = NULL, n_threads = 4)
#' 
runSVGsGpGp <- function(spe, x = NULL, fix_param_range = 0.5, n_neighbors = 15, 
                        n_threads = 1, verbose = FALSE) {
  
  stopifnot("logcounts" %in% assayNames(spe))
  
  if (!is.null(x)) stopifnot(nrow(x) == ncol(spe))
  
  y <- logcounts(spe)
  
  # --------
  # run GpGp
  # --------
  
  coords <- spatialCoords(spe)
  # scale coordinates proportionally
  # note: not needed when calculating ordering manually below
  #range_all <- max(apply(coords, 2, function(col) diff(range(col))))
  #coords <- apply(coords, 2, function(col) (col - min(col)) / range_all)
  
  # prepare for parallelized loop
  
  # calculate ordering only once
  ord <- order_middleout(coords)
  # re-parameterize ordering as ranking (if needed)
  #ord_rank <- match(seq_len(nrow(coords)), ord)
  # calculate nearest neighbors only once
  nn <- find_ordered_nn(coords[ord, ], m = n_neighbors)
  
  # run GpGp using parallelization
  ix <- seq_len(nrow(y))
  out_gpgp <- bplapply(ix, function(i) {
    # fit model (note: default if x is NULL is intercept-only model)
    y_i <- y[i, ]
    runtime <- system.time({
      # note: fixed parameter phi, manual reordering, using pre-calculated nearest neighbors
      out_i <- fit_model(y = y_i[ord], locs = coords[ord, ], X = x[ord, ], 
                         covfun_name = "exponential_isotropic", 
                         start_parms = c(0.1, fix_param_range, 0.1), fixed_parms = 2, 
                         NNarray = nn, reorder = FALSE, m_seq = n_neighbors, 
                         silent = !verbose)
    })
    res_i <- c(
      sigmasq = out_i$covparms[1], 
      tausq_sigmasq = out_i$covparms[3], 
      loglik = out_i$loglik, 
      runtime = runtime[["elapsed"]]
    )
    res_i
  }, BPPARAM = MulticoreParam(workers = n_threads))
  
  # collapse output list into matrix
  mat_gpgp <- do.call("rbind", out_gpgp)
  
  # add column containing standard parameterization of nugget, i.e. tausq
  mat_gpgp <- cbind(
    mat_gpgp, 
    tausq = mat_gpgp[, "tausq_sigmasq"] * mat_gpgp[, "sigmasq"]
  )
  
  # ------------------------------
  # calculate statistics and ranks
  # ------------------------------
  
  # mean logcounts
  mat_gpgp <- cbind(
    mat_gpgp, 
    mean = rowMeans(y)
  )
  
  # variance of logcounts
  mat_gpgp <- cbind(
    mat_gpgp, 
    var = rowVars(as.matrix(y))
  )
  
  # spatial coefficient of variation
  mat_gpgp <- cbind(
    mat_gpgp, 
    spcov = sqrt(mat_gpgp[, "sigmasq"]) / mat_gpgp[, "mean"]
  )
  
  # ratio of spatial to non-spatial variance
  mat_gpgp <- cbind(
    mat_gpgp, 
    ratio_sv = mat_gpgp[, "sigmasq"] / mat_gpgp[, "tausq"]
  )
  
  # proportion of spatial variance (out of total variance)
  mat_gpgp <- cbind(
    mat_gpgp, 
    prop_sv = mat_gpgp[, "sigmasq"] / (mat_gpgp[, "sigmasq"] + mat_gpgp[, "tausq"])
  )
  
  # return in rowData of spe object
  
  stopifnot(nrow(spe) == nrow(mat_gpgp))
  
  rowData(spe) <- cbind(rowData(spe), mat_gpgp)
  
  spe
}

