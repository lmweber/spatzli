#' plotQCspots
#'
#' Plots for quality control (QC)
#'
#' Functions to generate plots for quality control (QC) purposes.
#'
#' This function generates a plot showing spots to be discarded on the physical
#' x-y coordinates of the tissue slide.
#'
#'
#' @param spe Input object (SpatialExperiment or SingleCellExperiment).
#'
#' @param x_coord Name of column in colData containing x coordinates. Default =
#'   x_coord. Note that the column name needs to be provided as a variable name,
#'   not a character string (due to the way ggplot2 works).
#'
#' @param y_coord Name of column in colData containing y coordinates. Default =
#'   y_coord. Note that the column name needs to be provided as a variable name,
#'   not a character string (due to the way ggplot2 works).
#'
#' @param discard Name of column containing logical entries identifying spots to
#'   be discarded. Default = discard. Note that the column name needs to be
#'   provided as a variable name, not a character string (due to the way ggplot2
#'   works).
#'
#'
#'
#' @return Returns ggplot object containing plot. Returning the plot as an
#'   object allows users to add additional ggplot elements (e.g. new title or
#'   different formatting).
#'
#'
#' @importFrom SingleCellExperiment colData
#' @importFrom magrittr "%>%"
#' @importFrom ggplot2 ggplot aes geom_point coord_fixed scale_color_manual
#'   ggtitle theme_bw theme element_blank
#'
#' @export
#'
#' @examples
#' TO DO
#' 
plotQCspots <- function(spe, 
                        x_coord = x_coord, y_coord = y_coord, 
                        discard = discard) {
  
  # note: using quasiquotation to allow custom variable names in ggplot ("enquo"
  # and "!!" below, and arguments provided as variable names instead of
  # character strings)
  
  x_coord <- enquo(x_coord)
  y_coord <- enquo(y_coord)
  discard <- enquo(discard)
  
  as.data.frame(colData(spe)) %>% 
    ggplot(aes(x = !!x_coord, y = !!y_coord, color = !!discard)) + 
    geom_point(size = 0.9) + 
    coord_fixed() + 
    scale_color_manual(values = c("gray90", "red")) + 
    ggtitle("Quality control: discarded spots") + 
    theme_bw() + 
    theme(panel.grid = element_blank(), 
          axis.title = element_blank(), 
          axis.text = element_blank(), 
          axis.ticks = element_blank())
  
}

