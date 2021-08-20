% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/runSVGsGpGp.R
\name{runSVGsGpGp}
\alias{runSVGsGpGp}
\title{runSVGsGpGp}
\usage{
runSVGsGpGp(
  spe,
  x = NULL,
  fix_param_range = 0.5,
  n_neighbors = 15,
  lr_test = FALSE,
  n_threads = 1,
  verbose = FALSE
)
}
\arguments{
\item{spe}{\code{SpatialExperiment} Input object, assumed to be a
\code{SpatialExperiment} containing an assay named \code{logcounts} and
spatial coordinates accessible with \code{spatialCoords()}.}

\item{x}{\code{numeric matrix} Matrix of covariates, with number of rows
(spots) matching the number of columns (spots) in \code{spe}. Default =
NULL, which is an intercept-only model. See \code{GpGp} documentation for
more details.}

\item{fix_param_range}{\code{numeric} Fixed parameter value to use for
'range' covariance function parameter in 'exponential_isotropic' covariance
function (corresponding to 'phi' in other parameterizations). Default =
0.5. Set to NULL to estimate parameter instead. See \code{GpGp}
documentation for details on parameterization.}

\item{n_neighbors}{\code{numeric} Number of nearest neighbors. See
\code{GpGp} documentation for details.}

\item{lr_test}{\code{logical} Whether to calculate log likelihoods for model
without spatial terms for likelihood ratio test. Default = FALSE.}

\item{n_threads}{\code{integer} Number of threads for parallelization.
Default = 1.}

\item{verbose}{\code{logical} Whether to display verbose output from
\code{GpGp}. Default = FALSE.}
}
\value{
Returns summary statistics and SVG ranks as new columns in
  \code{rowData} in \code{spe} object.
}
\description{
Run method to identify spatially variable genes (SVGs) using GpGp
}
\details{
Identify top SVGs using GpGp (Guinness 2018).

This function runs GpGp separately for each gene, using parallelization for
faster runtime using one core per GpGp run. The main outputs of interest are
the covariance parameter estimates stored in 'covparms' in the GpGp output
(variance, range, nugget, if using 'exponential_isotropic' covariance
function). We fix the 'range' parameter (i.e. phi) so it is equal across
genes. In addition, the output 'loglik' can be used for likelihood ratio
tests.

Note parameterization used by GpGp: total variance = sigmasq + (tausq *
sigmasq); i.e. the nugget is 'tausq * sigmasq', not 'tausq' itself.

Assumes the input object is a \code{SpatialExperiment} containing an assay
named \code{logcounts} and filtered to exclude low-expressed genes, e.g. as
prepared with \code{\link{preprocessSVGs}}.
}
\examples{
library(SpatialExperiment)
library(STexampleData)
library(spatzli)

spe <- Visium_humanDLPFC()

spe <- preprocessSVGs(spe)

# subset 1 gene
spe_1 <- spe[1, ]
system.time({
  spe_1 <- runSVGsGpGp(spe_1, x = NULL, n_threads = 1, verbose = TRUE)
})

# subset 100 genes
# spe_100 <- spe[1:100, ]
# spe_100 <- runSVGsGpGp(spe_100, x = NULL, n_threads = 4)

}