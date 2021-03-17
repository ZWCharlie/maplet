#' mtw_preprocess
#'
#' Preprocess on data, including batch correction, quotient normalization, log transformation and kNN imputations
#'
#' @param D \code{SummarizedExperiment} input
#' @param batch_col Batch column
#' @param batch_fun Batch function to use, should be rather median or combat method. Default is \code{mt_pre_batch_median}
#' @param batch_loc Where should batch normalization take place, should be 'begin'/'middle'/'end', 'begin' means before quotient normalization, 'middle' means between quotient and kNN imputation, 'end' means after imputation. Default is 'begin'
#' @param batch_ref_samples Expression to filter out reference samples to use from rowData. Default is NULL
#' @param boxplot_options A list of parameters for \code{mt_plots_sample_boxplot}
#' @param quot_options A list of parameters for \code{mt_pre_norm_quot}
#' @param dilut_fac_cols List of names in colData column to correlate dilution factors with
#' @param dilution_options A list of parameters for \code{mt_plots_qc_dilutionplot}
#' @param log_base Base for log function, default is 2
#' @param do_impute Logical. Whether to run imputation, default is TRUE
#' @param knn_options A list of parameters for \code{mt_pre_impute_knn}
#'
#' @return D that has been preprocessed
#' @examples
#' \dontrun{... %>% mtw_preprocess() %>% ...}
#' \dontrun{... %>% mtw_preprocess(batch_col = "BATCH_MOCK", 
#' batch_loc = 'begin', dilut_fac_cols = c("num1","Group"), 
#' quot_options = list(ref_samples = quote(Group == 'Vehicle'))) %>% ...}
#'
#' @author Zeyu Wang
#'
#' @importFrom magrittr %>% %<>%
#' @import SummarizedExperiment
#'
#' @export


### Add choice to batch function -> when to run
mtw_preprocess <-
  function(D,
           batch_col,
           batch_fun = mt_pre_batch_median,
           batch_loc = 'begin',
           batch_ref_samples = NULL,
           boxplot_options = list(),
           quot_options = list(),
           dilut_fac_cols = c(),
           dilution_options = list(),
           log_base = 2,
           do_impute = T,
           knn_options = list()) {
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
    
    batch_norm <- function(D){
      # Batch normalization based on ref_sample
      if(!is.null(batch_ref_samples)){
        D %<>%
          batch_fun(batches = batch_col, ref_samples = dplyr::enquo(batch_ref_samples))
      } else{
        D %<>%
          batch_fun(batches = batch_col)
      }
    }
    
    run_batch = F
    
    boxplot_options_def = list(
      title = "Sample boxplot",
      show_legend = T,
      ylabel = "Metabolite concentrations",
      plot_logged = F,
      ggadd = NULL
    )
    
    quot_options_def = list(
      vars = 1:dim(D)[1],
      na_err = F,
      ref_samples = NULL,
      feat_max = 1
    )
    
    dilution_options_def = list(boxplot = T, ggadd = NULL)
    
    knn_options_def = list(method = "knn.obs.euc.sel", k = 10, 
                           verbose = F, use_multicore = F, n_cores = 5)
    
    # Boxplot at first
    boxplot_options <- map_lists(boxplot_options_def, boxplot_options)
    boxplot_options$D <- D
    D <- do.call("mt_plots_sample_boxplot", boxplot_options)
    
    # Batch function check
    if (!missing(batch_col)) {
      if(!(batch_loc == 'begin' || batch_loc == 'middle' || batch_loc == 'end')){
        stop('Batch location must be "begin"/"middle"/"end" ')
      }
      run_batch = T
    }
    
    # Batch at the beginning
    if (run_batch == T && batch_loc == 'begin') {
      # Batch normalization
      D %<>% batch_norm()
      
      # Boxplot after batch normalization
      boxplot_options$D <- D
      D <- do.call("mt_plots_sample_boxplot", boxplot_options)
      run_batch = F
    }
    
    # Quotient normalization
    quot_options <- map_lists(quot_options_def, quot_options)
    quot_options$D <- D
    if (!is.null(quot_options$ref_samples)) {
      D <- do.call("mt_pre_norm_quot", quot_options)
    } else{ # if NULL, remove the ref_sample
      quot_options$ref_samples <- NULL
      D <- do.call("mt_pre_norm_quot", quot_options)
    }
    
    # Dilution plot
    # check if there is any correlation between dilution factors and outcomes (bad sign if so)
    if (!missing(dilut_fac_cols)) {
      dilution_options <- map_lists(dilution_options_def, dilution_options)
      
      for (s in dilut_fac_cols) {
        dilution_options$D <- D
        dilution_options$in_col <- s
        D <- do.call("mt_plots_dilution_factor", dilution_options)
      }
    }
    
    # Boxplot after quotient
    boxplot_options$D <- D
    D <- do.call("mt_plots_sample_boxplot", boxplot_options)
    
    # Batch after quot norm
    if (run_batch == T && batch_loc == 'middle') {
      D %<>% batch_norm()
      
      # Boxplot after batch normalization
      boxplot_options$D <- D
      D <- do.call("mt_plots_sample_boxplot", boxplot_options)
      run_batch = F
    } 
    
    
    # logging
    D %<>%  mt_pre_trans_log(base = log_base)
    
    # KNN imputation
    if (do_impute == T) {
      knn_options <- map_lists(knn_options_def, knn_options)
      knn_options$D <- D
      D <- do.call("mt_pre_impute_knn", knn_options)
      boxplot_options$D <- D
      D <- do.call("mt_plots_sample_boxplot", boxplot_options)
    }
    
    # Batch at end
    if (run_batch == T && batch_loc == 'end') {
      # Batch normalization
      D %<>% batch_norm()
      
      # Boxplot after batch normalization
      boxplot_options$D <- D
      D <- do.call("mt_plots_sample_boxplot", boxplot_options)
    }
    D
  }