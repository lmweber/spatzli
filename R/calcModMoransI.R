#' calcModMoransI
#' 
#' Calculate modified Moran's I statistic to rank spatially variable genes
#' (SVGs).
#' 
#' Fast implementation of modified Moran's I statistic for ranking spatially
#' variable genes (SVGs).
#' 
#' We modify the definition of Moran's I statistic to make two
#' sparsity-preserving assumptions, due to the sparse nature of spatial
#' transcriptomics data and to speed up runtime.
#' 
#' - We assume that most genes are not detected in most spots (spatial
#' coordinates), and perform calculations using only the values from the
#' non-zero spots. For example, mean expression of each gene (which is used
#' inside the Moran's I formula) is calculated as the mean of the non-zero
#' spots.
#' 
#' - We assume that weights (in the spatial covariance matrix) below some
#' threshold (e.g. below 1% of the maximum weights value) are equal to zero.
#' 
#' 
#' @param spe Input object (SpatialExperiment). Assumed to contain an assay
#'   named "logcounts" containing log-transformed normalized counts in sparse
#'   matrix format, and "spatialCoords" slot containing spatial coordinates.
#' 
#' @param l_prop Value to set characteristic length parameter in squared
#'   exponential kernel used to calculate weights matrix. The characteristic
#'   length parameter is set to "l_prop" times the maximum range of the x or y
#'   coordinates. Default = 0.1.
#' 
#' @param weights_min Minimum weights threshold. Weights (in the spatial
#'   covariance matrix) that are below "weights_min" times the maximum weights
#'   value are assumed to be zero.
#' 
#' @param x_coord Name of column in spatialCoords slot containing x-coordinates.
#'   Default = "pxl_row_in_fullres".
#' 
#' @param y_coord Name of column in spatialCoords slot containing x-coordinates.
#'   Default = "pxl_col_in_fullres".
#' 
#' @param verbose Whether to print messages. Default = FALSE.
#' 
#' 
#' @return Returns a list containing output values (one value per gene).
#' 
#' 
#' @importFrom SingleCellExperiment logcounts
#' @importFrom SummarizedExperiment assayNames
#' @importFrom kernlab laplacedot kernelMatrix
#' @importFrom methods as
#' 
#' @export
#' 
#' @examples
#' paste0("to do")
#' 
calcModMoransI <- function(spe, l_prop = 0.1, weights_min = 0.01, 
                           x_coord = "pxl_row_in_fullres", y_coord = "pxl_col_in_fullres", 
                           verbose = FALSE) {
  
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
  # change to alternative parameterization for Laplacian kernel in kernlab package
  sigma_default <- 1 / l_default
  # define kernel function
  kernel_rbf <- laplacedot(sigma = sigma_default)
  # calculate kernel weights matrix using kernlab package (fast)
  weights <- kernelMatrix(kernel_rbf, xy_coords)
  
  # sparsity-preserving assumption: assume weights below below the threshold 
  # (i.e. 'weights_min * max weights value') are equal to zero
  thresh <- weights_min * max(weights)
  if (verbose) {
    message(paste0("sparsity preserving assumption: ", 
                   round(mean(weights < thresh) * 100), 
                   "% of values in weights matrix assumed to be zero"))
  }
  weights[weights < thresh] <- 0
  
  # convert weights matrix to sparse flattened vector
  weights <- matrix(as.vector(weights), nrow = 1)
  weights <- as(weights, "dgCMatrix")
  
  
  # --------------------------------------
  # calculate modified Moran's I statistic
  # --------------------------------------
  
  # lots of tricks here to speed up runtime
  
  x <- logcounts(spe)
  
  # initialize with zeros since many genes have all values equal to zero; skip
  # calculations for these genes
  n_nonzero <- rep(0, nrow(x))
  means <- rep(0, nrow(x))
  vars <- rep(0, nrow(x))
  stats_novar <- rep(0, nrow(x))
  stats <- rep(0, nrow(x))
  
  # convert weights to non-sparse format (but still in flattened format) for
  # faster multiplication inside loop
  weights_mx <- as.matrix(weights)
  stopifnot(nrow(weights_mx) == 1)
  
  # calculate for each gene
  runtime <- system.time({
    for (i in seq_len(nrow(x))) {
      if (verbose) print(i)
      # subset gene i; drop = FALSE is required to keep as sparse matrix
      x_i <- x[i, , drop = FALSE]
      # if all values equal to zero then skip to next gene (for faster runtime)
      if (sum(x_i) == 0) next
      # number of nonzero values
      n_nonzero_i <- length(x_i@x)
      n_nonzero[i] <- n_nonzero_i
      # mean expression that takes sparsity into account (i.e. mean for spots with 
      # non-zero values only)
      mean_i <- mean(x_i@x)
      means[i] <- mean_i
      y_i <- x_i
      y_i@x <- y_i@x - mean_i
      # calculate Kronecker product (mathematically equivalent to flattened version 
      # of outer product; much faster to calculate)
      yy_i <- kronecker(y_i, y_i)
      # multiply by flattened weights vector (non-sparse format for faster multiplication, 
      # but extracting only non-zero elements)
      # note zero-based indexing so need to add 1
      idx <- yy_i@j + 1
      weights_keep <- weights_mx[idx]
      wyy_i <- yy_i@x * weights_keep
      # sum of weights taking sparsity into account
      tot_weights_keep <- sum(weights_keep)
      # sum and divide by weights
      wyy_i_scaled <- sum(wyy_i) / tot_weights_keep
      stats_novar[i] <- wyy_i_scaled
      # variance (for non-zero points only)
      var_i <- sum((yy_i@x)^2) / n_nonzero_i
      vars[i] <- var_i
      # divide by variance
      stats[i] <- wyy_i_scaled / var_i
    }
  })
  
  # display runtime
  if (verbose) message(paste0("runtime: ", round(runtime[["elapsed"]]), " seconds"))
  
  # return output
  list(
    n_nonzero = n_nonzero, 
    means = means, 
    stats_novar = stats_novar, 
    vars = vars, 
    stats = stats, 
    runtime = runtime[["elapsed"]]
  )
  
}

