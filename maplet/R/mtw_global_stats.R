#' Wrapper function to run global statistics visualization
#'
#' Generating PCA, UMAP and heatmap plots of dataset
#'
#' @param D \code{SummarizedExperiment} input
#' @param cols_for_color List of names in colData to be colored in PCA and UMAP, column names should be quoted by quote()
#' @param do_pca Logical. Whether to generate PCA plot
#' @param pca_options A list of parameters for \code{mt_plots_pca}, expressions should be quoted by quote()
#' @param do_umap Logical. Whether to generate UMAP plot 
#' @param umap_options A list of parameters for \code{mt_plots_umap}, expressions should be quoted by quote()
#' @param do_heatmap Logical. Whether to generate heatmap
#' @param heatmap_options A list of parameters for \code{mt_plots_heatmap}, expressions should be quoted by quote()
#'
#' @return D with PCA, UMAP and heatmap plots
#' @examples
#' \dontrun{... %>% mtw_global_stats(cols_for_color = list(quote(Cohort)),
#' pca_options = list(scale_data = T, size = 2, ggadd=quote(scale_size_identity())),
#' umap_options = list(scale_data = T, size = 2, ggadd=quote(scale_size_identity())),
#' heatmap_options = list(scale_data = T, fontsize = 5)) %>% ...}
#'
#' @author Zeyu Wang
#'
#' @importFrom magrittr %>% %<>%
#' @import SummarizedExperiment
#' @export

mtw_global_stats <-
  function(D,
           cols_for_color=NULL,
           do_pca = T,
           pca_options = list(),
           do_umap = T,
           umap_options = list(),
           do_heatmap = T,
           heatmap_options = list()) {
    stopifnot("SummarizedExperiment" %in% class(D))
    
    # function to merge default and user input
    map_lists <- function(def, user) {
      # check if there are any illegal entries
      bad_args <- setdiff(names(user), names(def))
      if (length(bad_args) > 0) {
        stop(sprintf(
          "Invalid argument(s): %s",
          paste0(bad_args, collapse = ", ")
        ))
      }
      # now fill up entries in def
      def[names(user)] <- user
      # return
      def
    }
    
    pca_options_def = list(
      title = "PCA",
      scale_data = F,
      pc1 = 1,
      pc2 = 2,
      data_type = "scores",
      label_col = "",
      text_repel = T,
      ellipse = NA,
      exp_var_plot = F,
      store_matrices = F,
      ggadd = NULL,
      size = 2
    )
    
    umap_options_def = list(
      title = "UMAP",
      scale_data = F,
      label_col = "",
      text_repel = T,
      store_matrices = F,
      ggadd = NULL,
      n_neighbors = 15,
      size = 2
    )
    
    heatmap_options_def = list(
      scale_data = F,
      sym_zero = F,
      annotation_row = NA,
      annotation_col = NA,
      fontsize = 10,
      silent = TRUE,
      ggadd = NULL,
      return_gg = T,
      gg_scale = 1,
      # gg_ymin = quote(1 - gg_scale),
      # gg_xmin = quote(1 - gg_scale),
      # gg_xmax = quote(gg_scale),
      # gg_ymax = quote(gg_scale),
      cutree_rows = NA, 
      cutree_cols = NA,
      clustering_distance_rows = "euclidean",
      clustering_distance_cols = "euclidean", 
      clustering_method = "complete",
      color = (grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name =
                                                                          "RdYlBu"))))(100)
    )
    
    
    # PCA
    if(do_pca == T){
      pca_options <- map_lists(pca_options_def, pca_options)
      if(!missing(cols_for_color)){
        for(s in cols_for_color){
          pca_options$D <- D
          pca_options$color <- s
          D <- do.call("mt_plots_pca", pca_options)
        }
      }else{
        pca_options$D <- D
        D <- do.call("mt_plots_pca", pca_options)
      }
    }
    
    # UMAP
    if(do_umap == T){
      umap_options <- map_lists(umap_options_def, umap_options)
      if(!missing(cols_for_color)){
        for(s in cols_for_color){
          umap_options$D <- D
          umap_options$color <- s
          D <- do.call("mt_plots_umap",umap_options)
        }
      }else{
        umap_options$D <- D
        D <- do.call("mt_plots_umap",umap_options)
      }
    }
    # Heatmap
    if(do_heatmap == T){
      heatmap_options <- map_lists(heatmap_options_def,heatmap_options)
      heatmap_options$D <- D
      D <- do.call("mt_plots_heatmap",heatmap_options)
    }
    D
  }
