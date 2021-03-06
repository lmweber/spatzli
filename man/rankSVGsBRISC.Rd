% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/preprocessSVGs.R, R/rankSVGsBRISC.R
\name{rankSVGsBRISC}
\alias{rankSVGsBRISC}
\title{preprocessSVGs}
\usage{
rankSVGsBRISC(spe, x = NULL, n_threads = 4, ...)

rankSVGsBRISC(spe, x = NULL, n_threads = 4, ...)
}
\arguments{
\item{spe}{Input object (SpatialExperiment). Assumed to contain assays named
logcounts", and spatial coordinates accessible with "spatialCoords()".}

\item{x}{Matrix of covariates for "BRISC_estimation()". Number of rows must
equal number of spots. See "?BRISC_estimation" for details. Default = NULL
(intercept-only model).}

\item{n_threads}{Number of threads for parallelization. Default = 4.}

\item{...}{Additional arguments to pass to "BRISC_estimation()".}

\item{filter_counts}{\code{integer} Filter out zero-expressed and
low-expressed genes by keeping genes with at least this number of UMI
counts in at least one spot. Default = 5.}

\item{filter_mito}{\code{logical} Whether to filter out mitochondrial genes.
(Mitochondrial genes are highly expressed and may rank highly as SVGs if
included, but are often not of primary biological interest.) Default =
TRUE.}

\item{gene_name}{\code{character} Column in \code{rowData} gene names if
\code{filter_mito = TRUE}. Default = \code{gene_name}.}
}
\value{
Returns a \code{SpatialExperiment} ready for \code{rankSVGsBRIC()}.

Returns summary statistics from BRISC as new columns in "rowData",
  which can be used to rank SVGs.
}
\description{
Preprocessing steps to run SVGs methods

Rank spatially variable genes (SVGs) using BRISC approach
}
\details{
Calculate ranking of spatially variable genes (SVGs) using the BRISC approach
(bootstrap for rapid inference on spatial covariances) developed by Saha and
Datta (2018).
}
\examples{
# to do

# to do

}
