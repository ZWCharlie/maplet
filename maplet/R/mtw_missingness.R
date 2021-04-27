#' Wrapper function to filter missingness
#'
#' Filtering missingness of both features and samples, and generating 3 missingness plots: before filtering, after filter features, and after filter samples
#'
#' @param D \code{SummarizedExperiment} input
#' @param plot_options A list of parameters for \code{mt_plots_qc_missingness}
#' @param filter_options A list of parameters for \code{mt_pre_filtermiss}, default thresholds for features and samples are 0.2 and 0.1
#'
#' @return D with filtered data
#'  column "Feat_Missing" will be added to rowData with missing information on features
#'  column "Sample_Missing" will be added to colData with missing information on samples
#' 
#' @examples
#' \dontrun{... %>% mtw_missingness() %>% ...}
#' \dontrun{... %>% mtw_missingness(plot_options = list(feat_max = 0.1, samp_max = 0.05), filter_options = list(feat_max = 0.1, samp_max = 0.05)) %>% ...}
#'
#' @author Zeyu Wang
#'
#' @importFrom magrittr %>% %<>%
#' @import SummarizedExperiment
#'
#' @export

## Add all options of each function used inside
mtw_missingness <- function (D,
                             plot_options = list(),
                             filter_options = list()) {
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
  
  plot_options_def = list(
    feat_max = NA,
    samp_max = NA,
    plot_features = T,
    plot_samples = F,
    sec_axis_features = F,
    sec_axis_samples = F,
    samp_label_col = NA,
    plot_data = T
  )
  
  filter_options_def = list(
    feat_max = 0.2,
    samp_max = 0.1,
    group_col = NA,
    report_filtered = F,
    report_sample_col = ""
  )
  
  # missingness plot before filtering
  plot_options <- map_lists(plot_options_def, plot_options)
  samp_max <- plot_options$samp_max
  plot_options$samp_max <- NULL
  plot_options$D <- D
  D <- do.call("mt_plots_missingness", plot_options)
  
  # filtering features
  filter_options <- map_lists(filter_options_def, filter_options)
  samp_max <- filter_options$samp_max
  filter_options$samp_max <- NULL
  filter_options$D <- D
  D <- do.call("mt_pre_filter_missingness", filter_options)
  
  # plot missingness after filtering features
  plot_options$D <- D
  D <- do.call("mt_plots_missingness", plot_options)
  
  # filtering samples
  filter_options$samp_max <- samp_max
  filter_options$feat_max <- NULL
  filter_options$D <- D
  D <- do.call("mt_pre_filter_missingness", filter_options)
  
  # plot missingness after filtering samples
  plot_options$plot_features <- F
  plot_options$plot_samples <- T
  plot_options$samp_max <- samp_max
  plot_options$feat_max <- NULL
  plot_options$D <- D
  D <- do.call("mt_plots_missingness", plot_options)
  
  #annotations
  D %<>%
  mt_anno_missingness(anno_type = "features", out_col = "Met_Missing") %>%
  mt_anno_missingness(anno_type = "samples", out_col = "Sample_Missing")
  
  D
}
