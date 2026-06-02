#' exports_plots
#' 
#' function to generate plots (no hierarchy) from list of labelled dataframes.
#' @param list_gated_data list of dataframes. Each dataframe has 3 columns: marker 1 values, marker 2 values and label column.
#' @param n_cores Number of cores to use. Default to 1.
#' @param path_output Path to the directory where to export the plots.
#' @param type_plot the user can choose between density (="dens) or label assignment visualization (="ML")
#' @param show_legend If True it shows the legend for the label assignment visualization. Default to True.
#' @param x_lab x-axis label.
#' @param y_lab y-axis label.
#' @param size_title_x Size x axis label.
#' @param size_title_y Size y axis label.
#' @param aspect_ratio Set aspect ratio. Default to NULL> If = 1, y and x axis ticks have same distance.
#' @param w_val width value. Default to 16.
#' @param h_val height value. Default to 10.
#' @param size_axis_text Size of ticks labels.
#' @param export_csv Export plot data as csv files. Default to False.
#' @param side_by_side Arrange dens plot and ML plot side-by-side. Default to False.
#' @return NULL
#' @keywords flowMagic
#' @export
#' @examples 
#' \donttest{exports_plots()}



exports_plots<-function(list_gated_data,path_output,n_cores=1,type_plot="dens",show_legend=T,x_lab="x",
                        y_lab="y",size_title_x=23,size_title_y=23,aspect_ratio=NULL,w_val=16,h_val=10,
                        size_axis_text=25,export_csv=F,side_by_side=F,...){
  start<-Sys.time()
  all_names<-names(list_gated_data)
 if(export_csv==T){
    list_n_gates_all_data <- parallel::mclapply(1:length(list_gated_data), 
                                      function(i) {
                                        name_current_file <- all_names[i]
                                        name_current_file <- stringr::str_remove(name_current_file, 
                                                                        ".csv")
                                        print(name_current_file)
                                        if ("df_test_original" %in% names(list_gated_data[[i]])) {
                                          df_p <- list_gated_data[[i]]$df_test_original
                                        }
                                        else {
                                          df_p <- list_gated_data[[i]]
                                        }
                                        
                                        print("---- export csv file")
                                        path_output_file <- paste0(path_output, sprintf("/%s.csv", 
                                                                                        name_current_file))
                                        write.csv(df_p,file = path_output_file,row.names = F)
                                                
                                        return(NULL)
                                      }, mc.cores = n_cores)
    
  }else{
    list_n_gates_all_data <- parallel::mclapply(1:length(list_gated_data), 
                                      function(i) {
                                        name_current_file <- all_names[i]
                                        name_current_file <- stringr::str_remove(name_current_file, 
                                                                        ".csv")
                                        print(name_current_file)
                                        if ("df_test_original" %in% names(list_gated_data[[i]])) {
                                          df_p <- list_gated_data[[i]]$df_test_original
                                        }
                                        else {
                                          df_p <- list_gated_data[[i]]
                                        }
                                        all_classes <- unique(df_p[, 3])
                                        all_classes <- all_classes[all_classes != 0]

                                        if(side_by_side==F){
                                          if (length(all_classes) == 0) {
                                            type_plot <- "ML"
                                          }
                                          plot_name <- tryCatch(magicPlot(df_p, type = type_plot, 
                                          show_legend = show_legend, x_lab = x_lab, y_lab = y_lab, 
                                          size_title_x = size_title_x, size_title_y = size_title_y, 
                                          aspect_ratio = aspect_ratio, size_axis_text = size_axis_text, 
                                          ...), error = function(e) {
                                            return(NULL)
                                          })
                                        }else if(side_by_side==T){
                                            if (!requireNamespace("patchwork", quietly = TRUE)) {
                                              stop("The 'patchwork' package is required for side by side export. Please install it.")
                                            }
                                            if (!"package:patchwork" %in% search()){
                                              library(patchwork)
                                            }
                                            plot_dens<-magicPlot(df_p, type = "dens",...)
                                            plot_ml<-magicPlot(df_p, type = "ML",...)
                                            plot_name<-plot_dens + plot_ml
                                            w_val<-30
                                            h_val<-20
                                        }
                                        print("---- export plot")
                                        path_output_file <- paste0(path_output, sprintf("/%s.png", 
                                                                                        name_current_file))
                                        ggsave(filename = path_output_file, plot = plot_name, 
                                               width = w_val, height = h_val)
                                        return(NULL)
                                      }, mc.cores = n_cores)
  }
  end<-Sys.time()
  time_taken<-end-start
  print("Execution time:")
  print(time_taken)
  print("Done")
}

#' export_raw_gs_plots
#' 
#' function to generate ungated plots from selected gs node
#' @param gs GatingSet
#' @param n_cores Number of cores to use. Default to 1.
#' @param path_output Path to the directory where to export the plots.
#' @param x_lab x-axis label. Default to NULL (inherited from GatingSet metadata).
#' @param y_lab y-axis label. Default to NULL (inherited from GatingSet metadata).
#' @param w_val width value. Default to 7 inches.
#' @param h_val height value. Default to 7 inches.
#' @param size_points Size points scatter plot.
#' @param return_data  If TRUE, return the list of dataframes used to generate the plots. Default to FALSE.
#' @return NULL
#' @keywords flowMagic
#' @export
#' @examples 
#' \donttest{export_raw_gs_plots()}

export_raw_gs_plots<-function(gs,node_name,channel_x,channel_y,path_output,n_cores=1,x_lab = NULL, y_lab = NULL, 
                               w_val = 7, h_val = 7,size_points=1,return_data=F,...){
  start <- Sys.time()
  samples_names <- sampleNames(gs)
 if(return_data==T){
    list_n_gates_all_data <- parallel::mclapply(1:length(samples_names), 
                                      function(i) {
                                        s <- samples_names[i]
                                        print(s)
                                        ff_temp <- gh_pop_get_data(gs[[s]], node_name)
                                        expr_matrix_ff <- exprs(ff_temp)
                                        df_exprs <- as.data.frame(expr_matrix_ff)
                                        df_exprs_selected_channels <- df_exprs[, c(channel_x, 
                                                                                   channel_y)]
                                        return(df_exprs_selected_channels)
                                      },mc.cores=n_cores)
    
    names(list_n_gates_all_data)<-samples_names
    return(list_n_gates_all_data)
  }else{
    list_n_gates_all_data <- parallel::mclapply(1:length(samples_names), 
                                      function(i) {
                                        s <- samples_names[i]
                                        print(s)
                                        ff_temp <- gh_pop_get_data(gs[[s]], node_name)
                                        expr_matrix_ff <- exprs(ff_temp)
                                        df_exprs <- as.data.frame(expr_matrix_ff)
                                        df_exprs_selected_channels <- df_exprs[, c(channel_x, 
                                                                                   channel_y)]
                                        
                                        plot_name <- tryCatch(magicPlot(df_exprs_selected_channels, 
                                                                        type = "no_gate", x_lab = x_lab, y_lab = y_lab, 
                                                                        size_points = size_points, ...), error = function(e) {
                                                                          return(NULL)
                                                                        })
                                        path_output_file <- paste0(path_output, sprintf("/%s.png", 
                                                                                        s))
                                        ggsave(filename = path_output_file, plot = plot_name, 
                                               width = w_val, height = h_val)
                                        return(NULL)
                                      }, mc.cores = n_cores)
  }

  end <- Sys.time()
  time_taken <- end - start
  print("Execution time:")
  print(time_taken)
  print("Done")
}





#' magicPlot_gs_hierarchy_all
#'
#' Plot and export GatingSet hierarchy figures for multiple samples.
#'
#' This function is a wrapper around `magicPlot_gs_hierarchy()`. It applies
#' `magicPlot_gs_hierarchy()` to multiple samples in a `GatingSet` and exports
#' one hierarchy-style plot per sample.
#'
#' For each sample, the output file name is generated automatically from the
#' sample name by removing the file extension and appending the selected output
#' extension.
#'
#' For example:
#'
#' \preformatted{
#' sample name: Flow 1_BIS002.fcs
#' output file: Flow 1_BIS002.tiff
#' }
#'
#' @param gs A `GatingSet` object.
#'
#' @param samples_id Character vector or `NULL`. Names of the samples to plot.
#'   If `NULL`, all samples in `sampleNames(gs)` are plotted. If a character
#'   vector is provided, all values must be present in `sampleNames(gs)`.
#'   Default is `NULL`.
#'
#' @param path_output Character. Directory where the exported hierarchy plots
#'   should be saved. This argument is different from `path_output` in
#'   `magicPlot_gs_hierarchy()`: here it is a directory, while in
#'   `magicPlot_gs_hierarchy()` it is the full output file path.
#'
#' @param extension Character. File extension used for exported plots. The value
#'   can be supplied with or without a leading dot, for example `".tiff"` or
#'   `"tiff"`. The extension is passed indirectly to `ggplot2::ggsave()` through
#'   the generated output file path. Default is `".tiff"`.
#'
#' @param return_plot_list Logical. If `FALSE`, return only a named character
#'   vector containing the exported file paths. If `TRUE`, return a list
#'   containing both the exported file paths and the plot objects returned by
#'   `magicPlot_gs_hierarchy()`. Default is `FALSE`.
#'
#' @param ... Additional arguments passed to `magicPlot_gs_hierarchy()`. Common
#'   examples include `n_col_wrap`, `size_points`, `auto_size`, `add_labels`,
#'   `size_pol_name`, `plot_width`, `plot_height`, `dpi`, `type`,
#'   `size_axis_text`, `size_title_x`, and `size_title_y`.
#'
#' @return If `return_plot_list = FALSE`, a named character vector containing
#'   the exported file paths, with names corresponding to sample names. If
#'   `return_plot_list = TRUE`, a list with the following elements:
#'   \describe{
#'     \item{paths}{A named character vector containing the exported file paths.}
#'     \item{plots}{A named list of plot objects returned by `magicPlot_gs_hierarchy()`.}
#'   }
#'
#' @details
#' The function first determines which samples should be plotted. If
#' `samples_id = NULL`, all samples in the `GatingSet` are used. Otherwise, the
#' requested sample names are checked against `sampleNames(gs)`.
#'
#' The output directory is created automatically if it does not already exist.
#' Each sample is then plotted by calling `magicPlot_gs_hierarchy()` with a
#' sample-specific full output file path.
#'
#' This wrapper is useful when the same hierarchy plot should be generated for
#' many or all samples in a `GatingSet`.
#'
#' @keywords flowMagic plotting GatingSet hierarchy export
#' @export
#'
#' @examples
#' \donttest{
#' # Export hierarchy plots for all samples
#' paths <- magicPlot_gs_hierarchy_all(
#'   gs = gs,
#'   path_output = "~/main/Results/Plots_gated_data_hierarchy"
#' )
#'
#' # Export hierarchy plots for selected samples only
#' paths <- magicPlot_gs_hierarchy_all(
#'   gs = gs,
#'   samples_id = sampleNames(gs)[1:3],
#'   path_output = "~/main/Results/Plots_gated_data_hierarchy"
#' )
#'
#' # Pass plotting options to magicPlot_gs_hierarchy()
#' paths <- magicPlot_gs_hierarchy_all(
#'   gs = gs,
#'   samples_id = sampleNames(gs)[1:3],
#'   path_output = "~/main/Results/Plots_gated_data_hierarchy",
#'   extension = ".tiff",
#'   type = "dens",
#'   add_labels = TRUE,
#'   dpi = 300
#' )
#'
#' # Return both exported paths and plot objects
#' out <- magicPlot_gs_hierarchy_all(
#'   gs = gs,
#'   samples_id = sampleNames(gs)[1:3],
#'   path_output = "~/main/Results/Plots_gated_data_hierarchy",
#'   return_plot_list = TRUE
#' )
#' }


magicPlot_gs_hierarchy_all <- function(gs,
                                       samples_id = NULL,
                                       path_output,
                                       extension = ".tiff",
                                       return_plot_list = FALSE,
                                       ...) {
  
  # =========================================================================
  # 1. Start workflow and validate input
  # =========================================================================
  # This function is a wrapper around magicPlot_gs_hierarchy().
  #
  # It runs magicPlot_gs_hierarchy() for multiple samples in a GatingSet and
  # exports one hierarchy plot per sample.
  #
  # The output file name is generated automatically from the sample name.
  #
  # Example:
  #   sample name: Flow 1_BIS002.fcs
  #   output file: Flow 1_BIS002.tiff
  
  message("$$$ Plot GatingSet hierarchy for multiple samples $$$")
  
  if (missing(path_output) || is.null(path_output)) {
    stop("path_output must be provided.")
  }
  
  if (!dir.exists(path_output)) {
    dir.create(path_output, recursive = TRUE)
  }
  
  
  # =========================================================================
  # 2. Select samples to plot
  # =========================================================================
  # If samples_id is NULL, plot all samples in the GatingSet.
  # Otherwise, plot only the selected samples.
  
  all_names <- sampleNames(gs)
  
  if (is.null(samples_id)) {
    
    samples_id <- all_names
    
  } else {
    
    missing_samples <- setdiff(samples_id, all_names)
    
    if (length(missing_samples) > 0) {
      stop(
        "These samples are not present in sampleNames(gs): ",
        paste(missing_samples, collapse = ", ")
      )
    }
  }
  
  message(sprintf("Number of samples to plot: %s", length(samples_id)))
  
  
  # =========================================================================
  # 3. Standardize file extension
  # =========================================================================
  # The user can provide extension either as ".tiff" or "tiff".
  # Internally, make sure the extension starts with ".".
  
  if (!startsWith(extension, ".")) {
    extension <- paste0(".", extension)
  }
  
  
  # =========================================================================
  # 4. Loop over samples and export one hierarchy plot per sample
  # =========================================================================
  # For each sample:
  #   1. remove the .fcs extension from the sample name,
  #   2. build the output file path,
  #   3. call magicPlot_gs_hierarchy(),
  #   4. optionally store the returned plot object.
  
  list_plots <- list()
  vec_paths <- character(0)
  
  for (n in samples_id) {
    
    message(sprintf("---- Plot hierarchy for sample: %s ----", n))
    
    n_final <- tools::file_path_sans_ext(n)
    path_file <- file.path(path_output, paste0(n_final, extension))
    
    message(sprintf("Output file: %s", path_file))
    
    p <- flowMagic::magicPlot_gs_hierarchy(
      gs = gs,
      sample_id = n,
      path_output = path_file,
      ...
    )
    
    list_plots[[n]] <- p
    vec_paths[n] <- path_file
  }
  
  
  # =========================================================================
  # 5. Return output
  # =========================================================================
  # By default, return only the vector of exported file paths.
  #
  # If return_plot_list = TRUE, return both:
  #   - exported file paths
  #   - plot objects
  
  message("Done")
  
  if (return_plot_list == TRUE) {
    return(
      list(
        paths = vec_paths,
        plots = list_plots
      )
    )
  }
  
  return(vec_paths)
}