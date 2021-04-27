#' mtw_apply_analysis
#'
#' Running standard analyses for multiple phenotypes, analysis steps including lm/cor analysis, 
#' multi-test with histgram, QQ-plot and volcano plot showing p-values, and box/scatter plot for the 
#' visualization of significant metabolites.
#' 
#'
#' @param D \code{SummarizedExperiment} input
#' @param outcomes A data frame, should contain two columns: outcomes$outcome: name of phenotype columns in colData; outcomes$type: data type of pheotype columns, must be one of numeric / binary / ordinal
#' @param stat_feat_filter Filter formula for metabolites, used as "stat_filter" in \code{mt_logging_statsinfo()} and "feat_filter" in \code{mt_plots_volcano()}, should be quoted by \code{quote()}
#' @param lm_options A list of parameters for \code{mt_stats_univ_lm()}, expressions should be quoted by \code{quote()}
#' @param cor_options A list of parameters for \code{mt_stats_univ_cor()}, expressions should be quoted by \code{quote()}
#' @param multtest_options A list of parameters for \code{mt_post_multtest()}, expressions should be quoted by \code{quote()}
#' @param volcano_options A list of parameters for \code{mt_plots_volcano()}, expressions and variable names should be quoted by \code{quote()}
#' @param box_scatter_options A list of parameters for \code{mt_plots_box_scatter()}, expressions and variable names should be quoted by \code{quote()}
#' 
#'
#' @return D with analysis results and plots of multiple phenotypes
#' @examples
#' \dontrun{... %>% mtw_apply_analysis(outcomes1, multtest_options = list(method = 'BH')) %>% ...}
#'
#' @author Zeyu Wang
#'
#' @importFrom magrittr %>% %<>%
#' @import SummarizedExperiment
#'
#' @export

mtw_apply_analysis <-   function(D,
                                 outcomes,
                                 stat_feat_filter = quote(p.adj < 0.05),
                                 lm_options = list(),
                                 cor_options = list(),
                                 multtest_options = list(),
                                 volcano_options = list(),
                                 box_scatter_options = list()) {
  
  stopifnot("SummarizedExperiment" %in% class(D))
  
  # help function to merge default and user input
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

  #nrow(outcomes)
  for (i in 1:nrow(outcomes)) {
    outcome <- outcomes$outcome[i]
    outcome_type <- outcomes$type[i]
    outcome_name <- paste0(outcome, "_", outcome_type)
    
    input_list <- list(
      stat_name = quote(outcome_name),
      is_bin = F,
      multtest_options = list(method = 'bonferroni'),
      stat_feat_filter = quote(stat_feat_filter),
      volcano_options = list(x = quote(statistic), 
                             feat_filter = stat_feat_filter,
                             #xlabel = quote(gsub("~", "", as.character(x))), 
                             vline = NA,
                             #hline,
                             ggadd = NULL,
                             colour = stat_feat_filter)
    )
    
    input_user <- list(
      lm_options = lm_options,
      cor_options = cor_options,
      multtest_options = multtest_options,
      volcano_options = volcano_options,
      box_scatter_options = box_scatter_options
    )
    
    #browser()
    # ordinal variables are actually modeled by their numeric counterparts and Kendall's tau
    if (outcome_type == 'ordinal') {
      D %<>% mt_reporting_heading(heading = paste0(outcome,'\n'), lvl = 2) 
      
      input_list_add <- list(
        D = D,
        test_type = 'cor',
        lm_options = list(),
        cor_options = list(method = 'kendall', in_col = outcome, stat_name = outcome_name, samp_filter = NULL,
                           exact = NULL),
        box_scatter_options = list(plot_type = 'box',
                                   x                  = sym(outcome),
                                   fill               = sym(outcome),
                                   #correct_confounder,
                                   feat_filter       = stat_feat_filter,
                                   feat_sort         = quote(p.value),
                                   annotation         = "{sprintf('P-value: %.2e', p.value)}\nPadj: {sprintf('%.2e', p.adj)}",
                                   full_info = F,
                                   text_size = 3.88,
                                   jitter = "beeswarm",
                                   restrict_to_used_samples = T,
                                   ylabel = NULL,
                                   fit_line = T,
                                   fit_line_se = T,
                                   ggadd = NULL)
      )
      
      # analysis wrapper
      input_list[names(input_list_add)] <- input_list_add
      
      for (s in names(input_user)) {
        input_list[[s]] <- map_lists(input_list[[s]], input_user[[s]])
      }
      # [TODO] Merge the following lists from user input
      # lm_options = list(),
      # multtest_options = list(),
      # volcano_options = list(),
      # box_scatter_options = list())
      # [TODO] Merge user input with input_list
      D <- do.call('mtw_analysis', input_list) # for numeric, twofactor, or multifactor: use lm
      
    } else if (outcome_type == 'numeric') {
      
      D %<>% mt_reporting_heading(heading = paste0(outcome,'\n'), lvl = 2)
      
      input_list_add <- list(
        D = D,
        test_type = 'lm',
        lm_options = list(formula = formula(sprintf("~ %s",outcome)),samp_filter = NULL, n_cores = 1),
        cor_options = list(),
        box_scatter_options = list(plot_type = 'scatter',
                                   x                  = sym(outcome),
                                   fill               = NULL,
                                   #correct_confounder,
                                   feat_filter       = stat_feat_filter,
                                   feat_sort         = quote(p.value),
                                   annotation         = "{sprintf('P-value: %.2e', p.value)}\nPadj: {sprintf('%.2e', p.adj)}",
                                   full_info = F,
                                   text_size = 3.88,
                                   jitter = "beeswarm",
                                   restrict_to_used_samples = T,
                                   ylabel = NULL,
                                   fit_line = T,
                                   fit_line_se = T,
                                   ggadd = NULL)
      )
      
      input_list[names(input_list_add)] <- input_list_add
      
      for (s in names(input_user)) {
        input_list[[s]] <- map_lists(input_list[[s]], input_user[[s]])
      }

      # [TODO] Merge user input with input_list
      # analysis wrapper
      D <- do.call('mtw_analysis', input_list)
      
    } else if (outcome_type == 'binary'){
      
      D %<>% mt_reporting_heading(heading = paste0(outcome,'\n'), lvl = 2)
      
      input_list_add <- list(
        D = D,
        test_type = 'lm',
        is_bin = T,
        lm_options = list(formula = formula(sprintf("~ %s",outcome)),samp_filter = NULL, n_cores = 1),
        cor_options = list(),
        box_scatter_options = list(plot_type = 'box',
                                   x                  = sym(outcome),
                                   fill               = sym(outcome),
                                   #correct_confounder,
                                   feat_filter       = stat_feat_filter,
                                   feat_sort         = quote(p.value),
                                   annotation         = "{sprintf('P-value: %.2e', p.value)}\nPadj: {sprintf('%.2e', p.adj)}",
                                   full_info = F,
                                   text_size = 3.88,
                                   jitter = "beeswarm",
                                   restrict_to_used_samples = T,
                                   ylabel = NULL,
                                   fit_line = T,
                                   fit_line_se = T,
                                   ggadd = NULL,
                                   fill = NULL)
      )
      
      # now fill up entries in def
      input_list[names(input_list_add)] <- input_list_add
      
      for (s in names(input_user)) {
        input_list[[s]] <- map_lists(input_list[[s]], input_user[[s]])
      }
      # [TODO] Merge the following lists from user input
      # lm_options = list(),
      # multtest_options = list(),
      # volcano_options = list(),
      # box_scatter_options = list())
      # [TODO] Merge user input with input_list
      # analysis wrapper
      D <- do.call('mtw_analysis', input_list)
    }
  }
  return(D)
}