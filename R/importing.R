#' import_sample_gated
#' 
#' convert a flowWorkspace or a GatingMl file of the train data folder into a gated gh (gh_train)
#' @param path path to gs or GatingML object.
#' @param type type of object.
#' @param group_wsp Group of wsp to import.
#' @return GatingSet object
#' @export
#' @examples 
#' \donttest{import_sample_gated()}

import_gating_info<-function(path,type="gs",group_wsp=NULL){
  if(type=="gs"){
    gs<-load_gs(path)
  }else if(type=="ws"){
    ws<-CytoML::open_flowjo_xml(path)
    gs<-flowjo_to_gatingset(ws,name=group_wsp)
  }
  return(gs)
}

#' import_reference_csv
#' 
#' function to import plain gold standards data (no hierarchy)
#' @param path_results path to directory containing the csv files  to read (with third column of labels).
#' @param n_cores Number of cores to use. Default to 1.
#' @return list of dataframes
#' @export
#' @examples 
#' \donttest{import_reference_csv()}


import_reference_csv<-function(path_results,n_cores=1){
  start<-Sys.time()
  path_data_expr<-list.files(path = path_results,full.names = T,recursive = F)
  
  names_plot<-list.files(path = path_results,full.names = F,recursive = F)
  list_data_plot<-parallel::mclapply(names_plot,function(n){
      ind_expr<-grep(n,path_data_expr,fixed=T)
      path_n_expr<-path_data_expr[ind_expr]
      data_plot<-read.csv(path_n_expr,check.names = F)
      return(data_plot)
  },mc.cores = n_cores)
  names(list_data_plot)<-names_plot
  end<-Sys.time()
  time_taken<-end-start
  print("Execution time:")
  print(time_taken)
  print("Done")
  return(list_data_plot)

}

#' import_test_set
#' 
#' read the ungated fcs files into a flowSet.
#' The ungated fcs are assumed to be already cleaned,compensated,and transformed.
#' @param path path of directory containig the fcs files.
#' @param n_samples Number of samples. Default to All.
#' @param ref_f_n Set reference flowFrame to match channel names. Default to 1(first flowFrame).
#' @return flowSet.
#' @export
#' @examples 
#' \donttest{import_test_set()}


import_test_set_fcs<-function(path,n_samples="All",ref_f_n=1){
  start<-Sys.time()
  paths_files<-list.files(path,full.names = T)
  if(is.character(n_samples)==F){
    paths_files<-paths_files[n_samples]
  }
  # get references parameters from the ref flowFrame (the first one by default)
  f_ref<-read.FCS(paths_files[ref_f_n])
  m_expr<-exprs(f_ref)
  channel_names_ref<-colnames(m_expr)
  # check if the paramaters of all input flowframes are identical to the ref parameters
  list_all_f<-list()
  c<-0
  for(i in 1:length(paths_files)){
    f<-read.FCS(paths_files[i])
    sample_name<-identifier(f)
    print(sprintf("importing sample %s",sample_name))
    m_expr_current_f<-exprs(f)
    channel_names_f<-colnames(m_expr_current_f)
    check_out<-all(channel_names_f==channel_names_ref)
    if(check_out==F){
      c<-c+1
      f<-NULL
      warning(sprintf("%s has different colnames from ref colnames. Return NULL",sample_name))
    }
    list_all_f[[sample_name]]<-f
  }
  fs <- as(list_all_f,"flowSet")
  print(sprintf("Excluded samples: %d",c))
  print(sprintf("Final samples: %d",length(paths_files)-c))
  end<-Sys.time()
  time_taken<-end-start
  print("Time of execution:")
  print(time_taken)
  print("Done")
  return(fs)
}


#' import_test_set_csv
#' 
#' function to import test set in csv format.
#' @param path_data path to directory containing csv files to read (third column is ignored).
#' @param n_cores Number cores. Default to 1.
#' @param xy_col Colnames equal to x and y. Default to True.
#' @return List of dataframes.
#' @export
#' @examples 
#' \donttest{import_test_set_csv()}

import_test_set_csv<-function(path_data,n_cores=1,xy_col=T){
  path_data_all<-list.files(path = path_data,full.names = T,recursive = F)
  names_plot<-list.files(path = path_data,full.names = F,recursive = F)
  list_test_data<-parallel::mclapply(path_data_all,function(p){
    df_expr<-read.csv(p,check.names = F)
    if(nrow(df_expr)==0){
      return(NULL)
    }
    df_expr<-df_expr[,c(1,2)]
    if(xy_col==T){
      colnames(df_expr)<-c("x","y")
    }
    return(df_expr)
  },mc.cores = n_cores)
  names(list_test_data)<-names_plot
  vec_check<-sapply(list_test_data,function(x){
    check_x<-is.null(x)
  })
  ind<-which(vec_check==T)
  if(length(ind)!=0){
    list_test_data<-list_test_data[-ind]
    
  }
  return(list_test_data)
}


#' get_train_data
#' 
#' function to import training data based on paths to files.
#' @param paths_file Vector of paths. Each path points toward a single csv file containin training info (labels and bivariate expression). paths_file can be also directly the list of dataframes containing labels and bivariate expression.
#' @param df_paths Dataframe containing the paths of file to read. The paths to data must be in the first column. The associated paths to classes are in second column.
#' @param n_cores Number of cores. Default to 1.
#' @param prop_down Proportion of events (downsampling). Default to NULL (downsampling using number of points).
#' @param n_points_per_plot Number of points for downsampling.
#' @param remove_class Vector of classes to ignore. Default to NULL.
#' @param normalize_data If True, data is normalized to 0-1 range. Default to True.
#' @param vec_col vector of columns names if the input dataframes have more than 3 columns. The third column name must always refer to the column with the gate label of each event. Default to NULL.
#' @return Dataframe.
#' @export
#' @examples 
#' \donttest{get_train_data()}


get_train_data<-function(paths_file=NULL,df_paths=NULL,n_cores=1,prop_down=NULL,remove_class=NULL,
                         n_points_per_plot=NULL,normalize_data=T,vec_col=NULL){
  start<-Sys.time()
  if(is.null(df_paths)==F){
    paths_file<-df_paths[,1]
  }
  
  names_dens_features<-c("n_peaks_m1","h_peak_m1_1","pos_peak_m1_1","start_peak_m1_1","end_peak_m1_1",
                         "h_peak_m1_2","pos_peak_m1_2","start_peak_m1_2","end_peak_m1_2",
                         "h_peak_m1_3","pos_peak_m1_3","start_peak_m1_3","end_peak_m1_3",
                         "h_peak_m1_4","pos_peak_m1_4","start_peak_m1_4","end_peak_m1_4",
                         "n_peaks_m2","h_peak_m2_1","pos_peak_m2_1","start_peak_m2_1","end_peak_m2_1",
                         "h_peak_m2_2","pos_peak_m2_2","start_peak_m2_2","end_peak_m2_2",
                         "h_peak_m2_3","pos_peak_m2_3","start_peak_m2_3","end_peak_m2_3",
                         "h_peak_m2_4","pos_peak_m2_4","start_peak_m2_4","end_peak_m2_4")
  
  list_dfs<-parallel::mclapply(1:length(paths_file),function(i){
    print(sprintf("plot_num:%s",i))
    #print("----- get or import dataframe with classes")
    if(is.list(paths_file)==F && is.null(df_paths)==T){
      # import df
      current_path<-paths_file[i]
      df<-read.csv(current_path)
    }else if(is.list(paths_file)==T && is.null(df_paths)==T){
      df<-paths_file[[i]]
    }else if(is.null(df_paths)==F){
      current_path_data<-df_paths[i,1]
      current_path_classes<-df_paths[i,2]
      df_data<-read.csv(current_path_data)
      df_classes<-read.csv(current_path_classes)
      if("cor_labels" %in% colnames(df_classes)){
        df<-cbind(df_data,df_classes$cor_labels)
      }else{
        df<-cbind(df_data,df_classes)
      }
    }else{
      stop("input not valid")
    }
    if(ncol(df)>3){
      warning("dataframe has more than three columns, checking vec_col argument")
      if(is.null(vec_col)==T || length(vec_col)!=3){
        stop("the input dataframes has length > 3 and vec_col format is not valid. 
             Please either make dataframes of 3 columns or indicate 3 valid columns names in vec_col argument. 
             Third column must contain the classes.")
      }
      df<-df[,vec_col]
      }
    colnames(df)<-c("x1_expr","x2_expr","classes")
    #show(magicPlot(df = df,type = "dens",size_points = 1))
    if(is.null(prop_down)==T & is.null(n_points_per_plot)==T){
      prop_down<-1
    }else if(is.null(prop_down)==T & is.null(n_points_per_plot)==F){
      prop_down<-(n_points_per_plot/nrow(df))
      if(prop_down>1){
        prop_down<-1
      }
    }
    # downsample df
    out_part<-caret::createDataPartition(y=factor(df[,"classes"]),times = 1,p = prop_down)
    df<-df[out_part$Resample1,]
    # remove some classes if needed
    if(is.null(remove_class)==F){
      inds_to_remove<-which((df$classes %in% remove_class)==T)
      if(length(inds_to_remove)!=0){
        new_df$classes[inds_to_remove]<-0
      }
    }
    # get density features
    if(normalize_data==T){
      df$x1_expr<-range01(df$x1_expr)
      df$x2_expr<-range01(df$x2_expr)
    }
    df$x1_expr<-round(df$x1_expr,2)
    df$x2_expr<-round(df$x2_expr,2)
    df_dens<-csv_to_dens(df = df,with_classes = F,n_coord = 50)
    if(normalize_data==F){
      vec_info_dens<-get_density_features(df_dens = df_dens,min_height = 0.00)
    }else{
      vec_info_dens<-get_density_features(df_dens = df_dens)
    }
    #show(magicPlot(df = new_df,type = "ML",size_points = 1))
    # add density features
    m_info_dens<-matrix(vec_info_dens,length(vec_info_dens),nrow(df))
    m_info_dens<-t(m_info_dens)
    df_info_dens<-as.data.frame(m_info_dens)
    colnames(df_info_dens)<-names_dens_features
    df<-cbind(df,df_info_dens)
    # add other info
    df$plot_num<-as.character(rep(i,nrow(df)))
    all_classes<-unique(df$classes)
    all_classes<-all_classes[all_classes!=0]
    n_gates<-length(all_classes)
    df$n_gates_info<-rep(n_gates,nrow(df))
    return(df)
    },mc.cores = n_cores)
  gc()
  df_train<-do.call(rbind,list_dfs)
  row.names(df_train)<-NULL
  end<-Sys.time()
  time_taken<-end-start
  print("Time of execution:")
  print(time_taken)
  print("Done")
  return(df_train)
  }


