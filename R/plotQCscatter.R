#' plotQCscatter
#' 
#' Plots for quality control (QC)
#' 
#' Functions to generate plots for quality control (QC) purposes.
#' 
#' This function generates a scatterplot comparing two QC metrics, e.g. number
#' of detected features vs. number of cells per spot.
#' 
#' 
#' @param spe Input object (SingleCellExperiment).
#' 
#' @param metric_x Name of column in colData containing QC metric to plot on
#'   x-axis (e.g. "cell_count" for number of cells per spot). Default =
#'   "cell_count".
#' 
#' @param metric_y Name of column in colData containing QC metric to plot on
#'   y-axis (e.g. "sum" for total library size, "detected" for number of
#'   expressed features). Default = "detected".
#' 
#' @param threshold_x If provided, a vertical line will be drawn at this
#'   x-value. Default = NULL.
#' 
#' @param threshold_y If provided, a horizontal line will be drawn at this
#'   y-value. Default = NULL.
#' 
#' @param trend Whether to include a smoothed trend (loess). Default = TRUE.
#' 
#' @param marginal Whether to include marginal histograms. Default = TRUE.
#' 
#' 
#' @return Returns a ggplot object. Additional plot elements can be added as
#'   ggplot elements (e.g. title, formatting).
#' 
#' 
#' @importFrom rlang sym "!!"
#' @importFrom SingleCellExperiment colData
#' @importFrom ggplot2 ggplot aes geom_point geom_hline geom_smooth ggtitle
#'   theme_bw
#' 
#' @export
#' 
#' @examples
#' # to do
#' 
plotQCscatter <- function(spe, 
                          metric_x = "cell_count", metric_y = "detected", 
                          threshold_x = NULL, threshold_y = NULL, 
                          trend = TRUE, marginal = TRUE) {
  
  # note: using quasiquotation to allow custom variable names in ggplot ("sym" and "!!")
  
  metric_x <- sym(metric_x)
  metric_y <- sym(metric_y)
  
  df <- as.data.frame(colData(spe))
  
  p <- ggplot(df, aes(x = !!metric_x, y = !!metric_y)) + 
    geom_point(size = 0.8) + 
    ggtitle("QC metrics") + 
    theme_bw()
  
  if (!is.null(threshold_x)) p <- p + geom_vline(xintercept = threshold_x, color = "red")
  if (!is.null(threshold_y)) p <- p + geom_hline(yintercept = threshold_y, color = "red")
  
  if (trend) p <- p + geom_smooth(method = "loess")
  
  # removed for now since does not return ggplot object
  # if (marginal) p <- ggMarginal(p, type = "histogram")
  
  p
  
}

