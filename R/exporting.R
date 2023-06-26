#' exports_plots
#' 
#' function to generate plots (no hierarchy).
#' @param list_gated_data list of dataframes. Each dataframe has 3 columns: marker 1 values, marker 2 values and label column.
#' @param n_cores Path to directory containing the expression data of the files analyzed.
#' @param path_output Path to the directory where to export the plots.
#' @param type_plot the user can choose between density (="dens) or label assignment visualization (="ML")
#' @param show_legend If True it shows the legend for the label assignment visualization. Default to True.
#' @param x_lab x-axis label.
#' @param y_lab y-axis label.
#' @param size_title_x Size x axis label.
#' @param size_title_y Size y axis label.
#' @return NULL
#' @export
#' @examples 
#' \donttest{exports_plots()}



exports_plots<-function(list_gated_data,path_output,n_cores=1,type_plot="dens",show_legend=T,x_lab="x",
                             y_lab="y",size_title_x=23,size_title_y=23){
  start<-Sys.time()
  all_names<-names(list_gated_data)
  list_n_gates_all_data<-mclapply(1:length(list_gated_data),function(i){
    name_current_file<-all_names[i]
    name_current_file<-str_remove(name_current_file,".csv")
    print(name_current_file)
    if("df_test_original" %in% names(list_gated_data[[i]])){
      df_p<-list_gated_data[[i]]$df_test_original
    }else{
      df_p<-list_gated_data[[i]]
    }
    plot_name<-tryCatch(magicPlot(df_p,type = type_plot,show_legend = show_legend,x_lab = x_lab,y_lab = y_lab,
                                  size_title_x = size_title_x,size_title_y=size_title_y),error=function(e){return(NULL)})
    # export plot in correct folder
    print("---- export plot")
    path_output_file<-paste0(path_output,sprintf("/%s.png",name_current_file))
    ggsave(filename = path_output_file,plot=plot_name)
    return(NULL)
  },mc.cores=n_cores)
  end<-Sys.time()
  time_taken<-end-start
  print("Execution time:")
  print(time_taken)
  print("Done")
}
