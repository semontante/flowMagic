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
#' @param w_val width value. Default to 7.
#' @param h_val height value. Default to 7.
#' @param size_axis_text Size of ticks labels.
#' @param export_csv Export plot data as csv files. Default to False.
#' @return NULL
#' @export
#' @examples 
#' \donttest{exports_plots()}



exports_plots<-function(list_gated_data,path_output,n_cores=1,type_plot="dens",show_legend=T,x_lab="x",
                        y_lab="y",size_title_x=23,size_title_y=23,aspect_ratio=NULL,w_val=7,h_val=7,
                        size_axis_text=25,export_csv=F,...){
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
#' @param x_lab x-axis label.
#' @param y_lab y-axis label.
#' @param w_val width value. Default to 7 inches.
#' @param h_val height value. Default to 7 inches.
#' @param size_points Size points scatter plot.
#' @param return_data  If TRUE, return the list of dataframes used to generate the plots. Default to FALSE.
#' @return NULL
#' @export
#' @examples 
#' \donttest{export_raw_gs_plots()}

export_raw_gs_plots<-function(gs,node_name,channel_x,channel_y,path_output,n_cores=1,x_lab = "x", y_lab = "y", 
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