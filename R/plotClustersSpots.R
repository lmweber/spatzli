#' plotClustersSpots
#'
#' Plots for clustering
#'
#' Functions to generate plots for clustering results.
#'
#' This function generates a plot showing cluster labels on the spots in the
#' physical x-y coordinates of the tissue slide.
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
#' @param cluster_id Name of column containing cluster IDs. Default =
#'   cluster_id. Note that the column name needs to be provided as a variable
#'   name, not a character string (due to the way ggplot2 works).
#'
#' @param palette Color palette. Available options are "libd_layer_colors" and
#'   "Okabe-Ito". To use a custom palette, provide a vector of hex codes for the
#'   colors.
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
plotClustersSpots <- function(spe, 
                              x_coord = x_coord, y_coord = y_coord, 
                              cluster_id = cluster_id, 
                              palette = "libd_layer_colors") {
  
  # note: using quasiquotation to allow custom variable names in ggplot ("enquo"
  # and "!!" below, and arguments provided as variable names instead of
  # character strings)
  
  x_coord <- enquo(x_coord)
  y_coord <- enquo(y_coord)
  cluster_id <- enquo(cluster_id)
  
  if (palette == "libd_layer_colors") {
    palette <- c("#F0027F", "#377EB8", "#4DAF4A", "#984EA3", "#FFD700", 
                 "#FF7F00", "#1A1A1A", "#666666" )
  }
  else if (palette == "Okabe-Ito") {
    palette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", 
                 "#0072B2", "#D55E00", "#CC79A7")
  }
  
  as.data.frame(colData(spe)) %>% 
    ggplot(aes(x = !!x_coord, y = !!y_coord, color = !!cluster_id)) + 
    geom_point(size = 0.9) + 
    coord_fixed() + 
    scale_color_manual(values = palette) + 
    ggtitle("Clustering: x-y space") + 
    theme_bw() + 
    theme(panel.grid = element_blank(), 
          axis.title = element_blank(), 
          axis.text = element_blank(), 
          axis.ticks = element_blank())
  
}

