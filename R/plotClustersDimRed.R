#' plotClustersDimRed
#'
#' Plots for clustering
#'
#' Functions to generate plots for clustering results.
#'
#' This function generates a plot showing cluster labels on the spots in the
#' dimension-reduced space (PCA or UMAP).
#'
#'
#' @param spe Input object (SpatialExperiment or SingleCellExperiment).
#'
#' @param type Type of dimension reduction to use. Options are "PCA" or "UMAP".
#'   Default = "UMAP".
#'
#' @param x_axis Name of column in reducedDim to use for x-axis. Default =
#'   "UMAP1".
#'
#' @param y_axis Name of column in reducedDim to use for y-axis. Default =
#'   "UMAP2".
#'
#' @param cluster_id Name of column containing cluster IDs. Default =
#'   "cluster_id".
#'
#' @param palette Color palette. Available options are "libd_layer_colors" and
#'   "Okabe-Ito". To use a custom palette, provide a vector of hex codes for the
#'   colors.
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
#' @importFrom ggplot2 ggplot aes geom_point coord_fixed scale_color_manual
#'   ggtitle theme_bw theme element_blank
#'
#' @export
#'
#' @examples
#' TO DO
#' 
plotClustersDimRed <- function(spe, 
                               type = "UMAP", 
                               x_axis = "UMAP1", y_axis = "UMAP2", 
                               cluster_id = "cluster_id", 
                               palette = "libd_layer_colors") {
  
  # note: using quasiquotation to allow custom variable names in ggplot ("sym" and "!!")
  
  x_axis <- sym(x_axis)
  y_axis <- sym(y_axis)
  cluster_id <- sym(cluster_id)
  
  if (palette == "libd_layer_colors") {
    palette <- c("#F0027F", "#377EB8", "#4DAF4A", "#984EA3", "#FFD700", 
                 "#FF7F00", "#1A1A1A", "#666666" )
  }
  else if (palette == "Okabe-Ito") {
    palette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", 
                 "#0072B2", "#D55E00", "#CC79A7")
  }
  
  as.data.frame(cbind(colData(spe), reducedDim(spe, type))) %>% 
    ggplot(aes(x = !!x_axis, y = !!y_axis, color = !!cluster_id)) + 
    geom_point(size = 0.8) + 
    scale_color_manual(values = palette) + 
    ggtitle(paste0("Clustering: ", type, " space")) + 
    theme_bw() + 
    theme(panel.grid = element_blank(), 
          axis.title = element_blank(), 
          axis.text = element_blank(), 
          axis.ticks = element_blank())
  
}

