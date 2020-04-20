#' plotSpotQC
#' 
#' Plots for spot-level quality control (QC)
#' 
#' Plotting functions for spot-level quality control (QC)
#' 
#' 
#' @param sce Input object (SingleCellExperiment). Currently assumed to have the
#'   structure from our DLPFC collaboration, i.e. specific rowData and colData
#'   slots. Will generalize it later.
#' 
#' @param sample_name Which sample to subset from the SingleCellExperiment.
#'   Assumes structure from our DLPFC collaboration, e.g. "151673" for sample
#'   151673.
#' 
#' @param plot_type Plot type: either "UMI" or "gene", to plot either total UMI
#'   counts per spot or total gene counts per spot vs. number of cells per spot.
#' 
#' @param save_plot Whether to save plot. If FALSE, plot will be displayed to
#'   screen instead.
#' 
#' @param filename Filename and path to save plot.
#' 
#' @param width Plot width (inches).
#' 
#' @param height Plot height (inches).
#' 
#' 
#' @return Saves a plot.
#' 
#' 
#' @importFrom SingleCellExperiment SingleCellExperiment colData
#' @importFrom magrittr '%>%'
#' @importFrom ggplot2 ggplot aes geom_point geom_smooth scale_x_continuous xlab ylab ggtitle theme_bw theme element_blank ggsave
#' @importFrom ggExtra ggMarginal
#' 
#' 
#' @export
#' 
#' @examples
#' # load("../../../K99_grant_application/data/Human_DLPFC_Visium_processedData_sce_scran.Rdata")
#' # plotSpotQC(sce, sample_name = "151673")
#' 
plotSpotQC <- function(sce, sample_name = "151673", plot_type = c("UMI", "gene"), 
                       save_plot = FALSE, filename = NULL, width = 5, height = 5) {
  
  plot_type <- match.arg(plot_type)
  
  stopifnot(sample_name %in% names(table(colData(sce)$sample_name)))
  
  # subset sample of interest
  sce_sub <- sce[, colData(sce)$sample_name == sample_name]
  
  # create plot (note: assumes SingleCellExperiment structure from our DLPFC collaboration)
  max_count <- max(colData(sce_sub)$cell_count)
  
  p1 <- colData(sce_sub) %>% 
    as.data.frame %>% 
    ggplot(aes(x = cell_count, y = sum_umi)) + 
    geom_point() + 
    geom_smooth(method = "loess") + 
    scale_x_continuous(breaks = seq(0, max_count, by = 2)) + 
    xlab("Cells per spot") + 
    ylab("UMI counts per spot") + 
    ggtitle("UMI counts vs. number of cells per spot", 
            subtitle = paste0("Sample ", sample_name)) + 
    theme_bw() + 
    theme(panel.grid = element_blank())
  
  p1 <- ggMarginal(p1, type = "histogram")
  
  if (is.null(filename)) {
    filename <- paste0("plot_umis_vs_ncells_", sample_name, ".pdf")
  }
  
  if (save_plot) {
    ggsave(filename, p1, width = width, height = height)
  } else {
    p1
  }
}


