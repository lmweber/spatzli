#' plotQCspots
#' 
#' Plots for quality control (QC)
#' 
#' Functions to generate plots for quality control (QC) purposes.
#' 
#' This function generates a plot identifying spots that do not meet QC metrics
#' (i.e. spots to be discarded before downstream analysis) on the physical x-y
#' coordinates of the tissue slide.
#' 
#' 
#' @param spe Input object (SingleCellExperiment).
#' 
#' @param x_coord Name of column in colData containing x-coordinates. Default =
#'   "x_coord".
#' 
#' @param y_coord Name of column in colData containing y-coordinates. Default =
#'   "y_coord".
#' 
#' @param discard Name of column in colData identifying spots to be discarded
#'   (TRUE/FALSE values). Default = "discard".
#' 
#' 
#' @return Returns a ggplot object. Additional plot elements can be added as
#'   ggplot elements (e.g. title, formatting).
#' 
#' 
#' @importFrom rlang sym "!!"
#' @importFrom SingleCellExperiment colData
#' @importFrom ggplot2 ggplot aes geom_point coord_fixed scale_color_manual
#'   ggtitle theme_bw theme element_blank
#' 
#' @export
#' 
#' @examples
#' # to do
#' 
plotQCspots <- function(spe, 
                        x_coord = "x_coord", y_coord = "y_coord", 
                        discard = "discard") {
  
  # note: using quasiquotation to allow custom variable names in ggplot ("sym" and "!!")
  
  x_coord <- sym(x_coord)
  y_coord <- sym(y_coord)
  discard <- sym(discard)
  
  df <- as.data.frame(colData(spe))
  
  p <- ggplot(df, aes(x = !!x_coord, y = !!y_coord, color = !!discard)) + 
    geom_point(size = 0.5) + 
    coord_fixed() + 
    scale_color_manual(values = c("gray85", "red")) + 
    ggtitle("QC: discarded spots") + 
    theme_bw() + 
    theme(panel.grid = element_blank(), 
          axis.title = element_blank(), 
          axis.text = element_blank(), 
          axis.ticks = element_blank())
  
  p
  
}

