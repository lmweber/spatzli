#' calcSpatialAutoCov
#' 
#' Calculate spatial autocovariance
#' 
#' Calculate spatial autocovariance for each gene, for use in ranking spatially
#' variable genes (SVGs), either directly or as part of the formula for Moran's
#' I statistic.
#' 
#' This is a fast implementation that makes use of sparsity and vectorized
#' calculations. Sparsity is due to both (i) spots with zero expression and (ii)
#' weights below a minimum threshold. This greatly speeds up runtime compared to
#' simpler implementations of Moran's I statistic.
#' 
#' 
#' @param spe Input object (SpatialExperiment). Assumed to contain an assay
#'   named "logcounts" containing log-transformed normalized counts in sparse
#'   matrix format, and "spatialCoords" slot containing spatial coordinates.
#' 
#' @param l_prop Value to set characteristic length parameter in squared
#'   exponential kernel used to calculate weights matrix. The characteristic
#'   length parameter is set to "l_prop" times the maximum range of the x or y
#'   coordinates. Default = 0.2.
#' 
#' @param weights_min Minimum weights threshold. Weights (in the spatial
#'   covariance matrix) that are below "weights_min" times the maximum weights
#'   value are assumed to be zero. Default = 0.05.
#' 
#' @param x_coord Name of column in spatialCoords slot containing x-coordinates.
#'   Default = "pxl_row_in_fullres".
#' 
#' @param y_coord Name of column in spatialCoords slot containing y-coordinates.
#'   Default = "pxl_col_in_fullres".
#' 
#' @param max_cores Maximum number of cores to use for parallelized evaluation.
#'   Default = 4.
#' 
#' @param verbose Whether to print messages. Default = FALSE.
#' 
#' 
#' @return Returns a list containing output values (one value per gene).
#' 
#' 
#' @importFrom SingleCellExperiment logcounts
#' @importFrom SummarizedExperiment assayNames
#' @importFrom kernlab rbfdot kernelMatrix
#' @importFrom parallel detectCores
#' @importFrom BiocParallel bpparam bplapply
#' @importFrom Matrix rowSums
#' @importFrom matrixStats rowVars
#' @importFrom methods as
#' 
#' @export
#' 
#' @examples
#' # to do
#' 
calcSpatialAutoCov <- function(spe, l_prop = 0.2, weights_min = 0.05, 
                               x_coord = "pxl_row_in_fullres", y_coord = "pxl_col_in_fullres", 
                               max_cores = 4, verbose = FALSE) {
  
  stopifnot("logcounts" %in% assayNames(spe))
  
  # ------------------------
  # calculate weights matrix
  # ------------------------
  
  # calculate weights matrix using squared exponential kernel with characteristic 
  # length parameter equal to 'l_prop * max range of x or y coordinates'
  
  # get x-y coordinates
  xy_coords <- as.matrix(spatialCoords(spe)[, c(x_coord, y_coord)])
  # calculate characteristic length parameter
  l_default <- max(c(abs(diff(range(xy_coords[, 1]))), abs(diff(range(xy_coords[, 2]))))) * l_prop
  # change to alternative parameterization for RBF kernel used in kernlab package
  sigma_default <- 1 / (2 * l_default^2)
  # define kernel function
  kernel_rbf <- rbfdot(sigma = sigma_default)
  # calculate kernel weights matrix using kernlab package (fast)
  weights <- kernelMatrix(kernel_rbf, xy_coords)
  
  # sparsity-preserving assumption: assume weights below the threshold (i.e.
  # 'weights_min * max weights value') are equal to zero
  thresh <- weights_min * max(weights)
  if (verbose) {
    message(paste0("sparsity preserving assumption: ", 
                   round(mean(weights < thresh) * 100), 
                   "% of values in weights matrix assumed to be zero"))
  }
  weights[weights < thresh] <- 0
  
  # convert weights matrix to flattened vector
  weights_vec <- matrix(as.vector(weights), nrow = 1)
  
  
  # --------------------------------
  # calculate spatial autocovariance
  # --------------------------------
  
  # several tricks here to speed up runtime
  
  # parallelizable function to calculate statistic for gene i
  calc_i <- function(i, x, w) {
    if (verbose) print(i)
    # subset gene i and convert to non-sparse format (note: drop = FALSE can 
    # keep sparse format if required)
    x_i <- x[i, ]
    # mean expression for gene i (including zeros)
    mean_i <- mean(x_i)
    # subtract mean
    y_i <- x_i - mean_i
    # calculate Kronecker product (mathematically equivalent to flattened 
    # version of outer product but much faster to calculate; note this is the 
    # slowest part of the calculation)
    yy_i <- as.numeric(kronecker(y_i, y_i))
    # multiply by flattened weights vector
    yy_i_w <- yy_i * w
    # sum values (don't divide by sum of weights here, since values approach 
    # machine precision)
    stat_i <- sum(yy_i_w)
    stat_i
  }
  
  # expression values
  x <- logcounts(spe)
  # remove genes with all values equal to zero (for faster runtime)
  zeros <- rowSums(x) == 0
  x_nonzero <- x[!zeros, , drop = FALSE]
  n_genes_nonzero <- nrow(x_nonzero)
  
  # number of cores
  n_cores <- min(detectCores(), max_cores)
  BPPARAM <- bpparam()
  if (BPPARAM$workers < n_cores) {
    BPPARAM$workers <- n_cores
  }
  
  # calculate values using parallelized function
  runtime <- system.time({
    stats_nonzero <- bplapply(1:n_genes_nonzero, calc_i, 
                              x = x_nonzero, w = weights_vec, 
                              BPPARAM = BPPARAM)
  })
  
  # return values in full-length vectors (including zeros for genes with all 
  # values equal to zero)
  n_genes <- nrow(x)
  stats <- rep(0, n_genes)
  stats[!zeros] <- unlist(stats_nonzero)
  
  # other values to return
  means <- rowMeans(as.matrix(x))
  vars <- rowVars(as.matrix(x))
  weights_sum <- sum(weights_vec)
  
  # display runtime
  if (verbose) message(paste0("runtime: ", round(runtime[["elapsed"]]), " seconds"))
  
  # return output
  list(
    n_genes = n_genes, 
    means = means, 
    vars = vars, 
    stats = stats, 
    weights = weights, 
    weights_sum = weights_sum, 
    runtime = runtime[["elapsed"]]
  )
  
}

