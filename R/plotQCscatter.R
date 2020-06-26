#' plotQCscatter
#'
#' Plots for quality control (QC)
#'
#' Functions to generate plots for quality control (QC) purposes.
#'
#' This function generates a scatterplot showing a QC metric (e.g. library size
#' or number of expressed features) vs. number of cells per spot.
#'
#'
#' @param spe Input object (SpatialExperiment or SingleCellExperiment).
#'
#' @param cell_count Name of column in colData containing number of cells per
#'   spot to be plotted on x axis. Default = "cell_count".
#'
#' @param metric Name of column in colData containing QC metric to be plotted on
#'   y axis. Common options are "sum" (usually refers to total library size) and
#'   "detected" (usually refers to number of expressed features). Default =
#'   "detected".
#'
#' @param threshold If provided, a horizontal line will be drawn at this
#'   threshold for the metric.
#'
#' @param trend Whether to include smoothed trend (loess). Default = FALSE.
#'
#' @param marginal Whether to include marginal histograms. Default = FALSE.
#'
#'
#' @return Returns ggplot object containing plot. Returning the plot as an
#'   object allows users to add additional ggplot elements (e.g. new title or
#'   different formatting).
#'
#'
#' @importFrom rlang sym "!!"
#' @importFrom SingleCellExperiment colData
#' @importFrom magrittr "%>%"
#' @importFrom ggplot2 ggplot aes geom_point geom_hline geom_smooth ggtitle
#'   theme_bw
#' @importFrom ggExtra ggMarginal
#'
#' @export
#'
#' @examples
#' TO DO
#' 
plotQCscatter <- function(spe, 
                          cell_count = "cell_count", metric = "detected", 
                          threshold = NULL, trend = FALSE, marginal = FALSE) {
  
  # note: using quasiquotation to allow custom variable names in ggplot ("sym" and "!!")
  
  cell_count <- sym(cell_count)
  metric <- sym(metric)
  
  p <- as.data.frame(colData(spe)) %>% 
    ggplot(aes(x = !!cell_count, y = !!metric)) + 
    geom_point(size = 0.8) + 
    ggtitle("QC metric vs. number of cells per spot") + 
    theme_bw()
  
  if (!is.null(threshold)) {
    p <- p + geom_hline(yintercept = threshold, color = "red")
  }
  
  if (trend) {
    p <- p + geom_smooth(method = "loess")
  }
  
  if (marginal) {
    p <- ggMarginal(p, type = "histogram")
  }
  
  p
  
}

