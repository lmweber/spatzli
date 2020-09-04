#' plotQCbar
#'
#' Quality control (QC) plots for spatial transcriptomics datasets.
#'
#' Functions to generate quality control (QC) plots for spatial transcriptomics
#' datasets.
#'
#' This function generates a barplot for one quality control (QC) metric, e.g.
#' number of cells per spot. For number of cells per spot, the barplot
#' highlights spots with zero cells, which are special values.
#'
#'
#' @param spe Input object (SpatialExperiment).
#'
#' @param metric_x Name of column in colData containing QC metric to plot on
#'   x-axis (e.g. "cell_count" for number of cells per spot). Default =
#'   "cell_count".
#'
#' @param highlight_zeros Whether to highlight bar for x = 0 (e.g. zero cells
#'   per spot, which is a special value).
#'
#'
#' @return Returns a ggplot object. Additional plot elements can be added as
#'   ggplot elements (e.g. title, formatting).
#'
#'
#' @importFrom rlang sym "!!"
#' @importFrom SingleCellExperiment colData
#' @importFrom ggplot2 ggplot aes geom_bar xlab ggtitle theme_bw
#'   scale_fill_manual
#'
#' @export
#'
#' @examples
#' # to do
#' 
plotQCbar <- function(spe, metric_x = "cell_count", 
                      highlight_zeros = TRUE) {
  
  # note: using quasiquotation to allow custom variable names in ggplot ("sym" and "!!")
  
  metric_x <- sym(metric_x)
  
  df <- as.data.frame(colData(spe))
  
  # identify spots with metric_x == 0 (e.g. zero cells per spot)
  if (highlight_zeros) {
    df$is_zero <- factor(df[, as.character(metric_x)] == 0)
  }
  
  p <- ggplot(df, aes(x = !!metric_x)) + 
    geom_bar(fill = "gray70") + 
    xlab(as.character(metric_x)) + 
    ggtitle("QC metrics") + 
    theme_bw()
  
  if (highlight_zeros) {
    p <- p + geom_bar(aes(fill = is_zero)) + 
      scale_fill_manual(values = c("gray70", "firebrick2"))
  }
  
  p
  
}

