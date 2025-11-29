
#' get_pops_hierarchy_list
#' 
#' function to get info from a list of hierarchical dataset (it reports all the pops for each level in a vector). 
#' @param hierarchical_list list of populations names for each level.
#' @return Vector of characters.
#' @keywords flowMagic
#' @export
#' @examples 
#' \donttest{get_pops_hierarchy_list()}

get_pops_hierarchy_list<-function(hierarchical_list){
  all_levels<-names(hierarchical_list)
  final_vec<-sapply(1:length(all_levels),function(i){
    level_name<-all_levels[i]
    current_level<-hierarchical_list[[level_name]]
    all_pops_current_level<-names(current_level)
    return(sprintf("%s;%s",level_name,all_pops_current_level))
  })
  final_vec<-unlist(final_vec)
  return(final_vec)
}

#' get_slot_hierarchy_list
#' 
#' function to access a slot of results from the hierarchy list based on the pop selected. 
#' @param hierarchical_list list of populations names for each level.
#' @param pop_selected Population name.
#' @return Element of a list.
#' @keywords flowMagic
#' @export
#' @examples 
#' \donttest{get_slot_hierarchy_list()}

get_slot_hierarchy_list<-function(hierarchical_list,pop_selected){
  vec_levelxpops<-get_pops_hierarchy_list(hierarchical_list)
  ind<-grep(pop_selected,vec_levelxpops,fixed = T)
  if(length(ind)>1){
    vec_pops<-c()
    for(element in vec_levelxpops){
      strsplitted<-strsplit(element,";")[[1]]
      pop_name<-strsplitted[2]
      vec_pops<-c(vec_pops,pop_name)
    }
    ind<-which((vec_pops %in% pop_selected) ==T)
  }
  if(length(ind)>1){
    stop("multiple matches for the pop_selected")
  }
  if(length(ind)==0){
    stop("pop not found in the hierarchy")
  }
  info_position<-vec_levelxpops[ind]
  strsplitted<-strsplit(info_position,";")[[1]]
  level_pop<-strsplitted[1]
  pop_name<-strsplitted[2]
  slot_selected_pop<-hierarchical_list[[level_pop]][[pop_name]]
  return(slot_selected_pop)
}

#' get_distance_loc_vs_test
#' 
#' function to compare the similiraty between current test data and local training set.
#' @param test_df Dataframe with bivariate marker expression.
#' @param loc_df Dataframe with bivariate marker expression.
#' @param show_plot show density comparison. Default to none.
#' @param nboot Number of permutations. Default to 50.
#' @return List of p-values.
#' @keywords flowMagic
#' @export
#' @examples 
#' \donttest{get_distance_loc_vs_test()}


get_distance_loc_vs_test<-function(test_df,loc_df,show_plot="none",nboot=50){
  set.seed(123)
  # get markers expression for test and train
  marker_1_expr_test<-test_df[,1]
  marker_2_expr_test<-test_df[,2]
  marker_1_expr_train<-loc_df[,1]
  marker_2_expr_train<-loc_df[,2]
  # get density of markers for test and train
  d1_test<-density(marker_1_expr_test,n=50)
  d2_test<-density(marker_2_expr_test,n=50)
  d1_train<-density(marker_1_expr_train,n=50)
  d2_train<-density(marker_2_expr_train,n=50)
  # get values correctly formatted
  vec_1<-rep(1,length(d1_train$y))
  vec_2<-rep(2,length(d1_test$y))
  vec_all_mark_1<-c(vec_1,vec_2)
  d_all_mark_1<-c(d1_train$y,d1_test$y)
  vec_1<-rep(1,length(d2_train$y))
  vec_2<-rep(2,length(d2_test$y))
  vec_all_mark_2<-c(vec_1,vec_2)
  d_all_mark_2<-c(d2_train$y,d2_test$y)
  # Execute significance test between the density for each marker
  output_mark_1<-sm.density.compare(d_all_mark_1,group = vec_all_mark_1,model = "equal",display=show_plot,nboot=nboot)
  output_mark_2<-sm.density.compare(d_all_mark_2,group = vec_all_mark_2,model = "equal",display=show_plot,nboot=nboot)
  pvalue_mark_1<-output_mark_1$p
  pvalue_mark_2<-output_mark_2$p
  # lower pvalue more significant differences between the two densities
  return(list(pvalue_mark_1=pvalue_mark_1,pvalue_mark_2=pvalue_mark_2))
  
}

#' get_weights_density_features
#' 
#' function to calculate final score based on density features.
#' @param df_scores Dataframe of distance scores.
#' @return Vector of numbers.
#' @keywords flowMagic
#' @export
#' @examples 
#' \donttest{get_weights_density_features()}


get_weights_density_features<-function(df_scores){
  score_test_data<-df_scores[nrow(df_scores),]
  models_scores<-df_scores[-nrow(df_scores),]
  vec_scores<-sapply(1:nrow(models_scores),function(m){
    models_scores_vec<-models_scores[m,]
    matrix_scores<-rbind(models_scores_vec,score_test_data)
    dist_score<-dist(matrix_scores,method = "euclidean")
    dist_score<-as.numeric(dist_score)
    return(dist_score)
  })
  return(vec_scores)
}

#' get_density_scores
#' 
#' function to get scores for distance template calculation.
#' @param df_template Dataframe of template markers expression values.
#' @param df_test Dataframe of test markers expression values.
#' @param select_density_features Select features to use. Default to NULL.
#' @return Matrix of numbers.
#' @keywords flowMagic
#' @export
#' @examples 
#' \donttest{get_density_scores()}


get_density_scores<-function(df_template,df_test,select_density_features=NULL){
  names_dens_features<-c("n_peaks_m1","h_peak_m1_1","pos_peak_m1_1","start_peak_m1_1","end_peak_m1_1",
                         "h_peak_m1_2","pos_peak_m1_2","start_peak_m1_2","end_peak_m1_2",
                         "h_peak_m1_3","pos_peak_m1_3","start_peak_m1_3","end_peak_m1_3",
                         "h_peak_m1_4","pos_peak_m1_4","start_peak_m1_4","end_peak_m1_4",
                         "n_peaks_m2","h_peak_m2_1","pos_peak_m2_1","start_peak_m2_1","end_peak_m2_1",
                         "h_peak_m2_2","pos_peak_m2_2","start_peak_m2_2","end_peak_m2_2",
                         "h_peak_m2_3","pos_peak_m2_3","start_peak_m2_3","end_peak_m2_3",
                         "h_peak_m2_4","pos_peak_m2_4","start_peak_m2_4","end_peak_m2_4")
  
  colnames(df_template)<-c("x1_expr","x2_expr")
  colnames(df_test)<-c("x1_expr","x2_expr")
  df_dens_template<-csv_to_dens(df = df_template,with_classes = F,n_coord = 50)
  df_dens_test<-csv_to_dens(df = df_test,with_classes = F,n_coord = 50)
  
  vec_info_dens_template<-get_density_features(df_dens = df_dens_template)
  vec_info_dens_test<-get_density_features(df_dens = df_dens_test)
  matrix_scores<-rbind(vec_info_dens_template,vec_info_dens_test)
  colnames(matrix_scores)<-names_dens_features
  if(is.null(select_density_features)==F){
    matrix_scores<-matrix_scores[,select_density_features]
  }
  return(matrix_scores)
  
}

#' get_dist_template
#' 
#' function to get distance between template and test data.
#' @param matrix_scores Matrix of density features generated by get_density_scores function.
#' @param dist_method Type of distance method calculation.
#' @return Number.
#' @keywords flowMagic
#' @export
#' @examples 
#' \donttest{get_dist_template()}

get_dist_template<-function(matrix_scores,dist_method="euclidean"){
  dist_score<-dist(matrix_scores,method = dist_method)
  dist_score<-as.numeric(dist_score)
  dist_score<-round(dist_score,2)
  return(dist_score)
}


#' get_density_features
#'
#' function to get density features only given a bivarite density csv.
#' @param df_dens Dataframe of density estimates for both markers.
#' @param min_height Minimum height of the peaks to consider.
#' @return Vector of numbers.
#' @keywords flowMagic
#' @export
#' @examples 
#' \donttest{get_density_features()}

get_density_features<-function(df_dens,min_height=0.06){
  density_m1<-df_dens[,1]
  density_m2<-df_dens[,2]
  # show(plot(density_m1,type="l"))
  # show(plot(density_m2,type="l"))

  # features first density
  matrix_peaks_m1<-pracma::findpeaks(density_m1,minpeakheight = min_height,minpeakdistance = 5)
  if(nrow(matrix_peaks_m1)>1){
    matrix_peaks_m1<-matrix_peaks_m1[order(matrix_peaks_m1[,2],decreasing = F),]
  }
  #print(matrix_peaks_m1)
  n_peaks_m1<-nrow(matrix_peaks_m1)
  if(n_peaks_m1<4){
    n_missing<-4-n_peaks_m1
    for(n in 1:n_missing){
      vec_0<-rep(0,4)
      matrix_peaks_m1<-rbind(matrix_peaks_m1,vec_0)
    }
  }else if(n_peaks_m1>4){
    matrix_peaks_m1<-matrix_peaks_m1[1:4,]
    n_peaks_m1<-nrow(matrix_peaks_m1)
  }
  info_all_peaks_m1<-c()
  vec_names_info_m1<-c()
  for(i in 1:nrow(matrix_peaks_m1)){
    info_peak_i_m1<-matrix_peaks_m1[i,]
    info_all_peaks_m1<-c(info_all_peaks_m1,info_peak_i_m1)
  }
  info_all_peaks_m1<-round(info_all_peaks_m1,2)
  info_all_peaks_m1<-c(n_peaks_m1,info_all_peaks_m1)
  # features second density
  matrix_peaks_m2<-pracma::findpeaks(density_m2,minpeakheight = min_height,minpeakdistance = 5)
  if(nrow(matrix_peaks_m2)>1){
    matrix_peaks_m2<-matrix_peaks_m2[order(matrix_peaks_m2[,2],decreasing = F),]
  }
  #print(matrix_peaks_m2)
  
  n_peaks_m2<-nrow(matrix_peaks_m2)
  if(n_peaks_m2<4){
    n_missing<-4-n_peaks_m2
    for(n in 1:n_missing){
      vec_0<-rep(0,4)
      matrix_peaks_m2<-rbind(matrix_peaks_m2,vec_0)
    }
  }else if(n_peaks_m2>4){
    matrix_peaks_m2<-matrix_peaks_m2[1:4,]
    n_peaks_m2<-nrow(matrix_peaks_m2)
  }
  info_all_peaks_m2<-c()
  vec_names_info_m2<-c()
  for(i in 1:nrow(matrix_peaks_m2)){
    info_peak_i_m2<-matrix_peaks_m2[i,]
    info_all_peaks_m2<-c(info_all_peaks_m2,info_peak_i_m2)
  }
  info_all_peaks_m2<-round(info_all_peaks_m2,2)
  info_all_peaks_m2<-c(n_peaks_m2,info_all_peaks_m2)
  
  # return final results
  final_vec<-c(info_all_peaks_m1,info_all_peaks_m2)
  return(final_vec)
}


#' csv_to_dens
#'
#' function to get density of events (with classes associated if present).
#' @param df Dataframe of marker expression values.
#' @param with_classes Consider classes. Default to True.
#' @param n_coord Grid size. Default to df.
#' @param normalize_data If True, data is normalized to 0-1 range. Default to True.
#' @return Dataframe.
#' @keywords flowMagic
#' @export
#' @examples 
#' \donttest{csv_to_dens()}

csv_to_dens<-function(df,with_classes=T,n_coord="df",normalize_data=T){

  if(with_classes==T){
      x1_expr<-df[,1]
      x1density<-(density(df[,1],n=nrow(df)))$y # density of x1 marker expression
      x2_expr<-df[,2]
      x2density<-(density(df[,2],n=nrow(df)))$y # density of x2 marker expression
      classes<- df[,3] 
      if(normalize_data==T){
        x1_expr<-range01(x1_expr)
        x2_expr<-range01(x2_expr)
      }
      new_df<-as.data.frame(cbind(x1_expr,x1density,x2_expr,x2density,classes))
  }else if(with_classes==F){
    if(n_coord=="df"){
      x1_expr<-df[,1]
      x1density<-(density(df[,1],n=nrow(df)))$y # density of x1 marker expression
      x2_expr<-df[,2]
      x2density<-(density(df[,2],n=nrow(df)))$y # density of x2 marker expression
      if(normalize_data==T){
        x1_expr<-range01(x1_expr)
        x2_expr<-range01(x2_expr)
      }
      new_df<-as.data.frame(cbind(x1_expr,x1density,x2_expr,x2density))
    }else{
      x1density<-(density(df[,1],n=n_coord))$y # density of x1 marker expression
      x2density<-(density(df[,2],n=n_coord))$y # density of x2 marker expression
      new_df<-as.data.frame(cbind(x1density,x2density))
    }
  }
  return(new_df)
}


#' get_classes_expr_df
#'
#' function to get classes of original expression df based on density df predictions.
#' @param dens_df Dataframe of density values.
#' @param original_df Original dataframe.
#' @return Dataframe.
#' @keywords flowMagic
#' @export
#' @examples 
#' \donttest{get_classes_expr_df()}


get_classes_expr_df<-function(dens_df,original_df){
  list_original_dfs_all_classes<-list()
  all_classes<-unique(dens_df$classes)
  all_classes<-all_classes[all_classes!=0]
  for(dens_class in all_classes){
    print(dens_class)
    dens_df_class<-dens_df[dens_df$classes==dens_class,]
    # get classes for x1 coordinates
    min_thr_x1<-min(dens_df_class$x1expr_coord)
    max_thr_x1<-max(dens_df_class$x1expr_coord)

    inds<-which((original_df$x >= min_thr_x1) & (original_df$x <= max_thr_x1))
    original_df_temp<-original_df[inds,]

    # get classes for x2 coordinates
    
    min_thr_x2<-min(dens_df_class$x2expr_coord)
    max_thr_x2<-max(dens_df_class$x2expr_coord)

    inds<-which((original_df_temp$y >= min_thr_x2) & (original_df_temp$y <= max_thr_x2))
    original_df_class<-original_df_temp[inds,]
    original_df_class$classes<-rep(dens_class,nrow(original_df_class))
    
    # return original df expression of current class
    list_original_dfs_all_classes[[dens_class]]<-original_df_class
    
    
  }
  df_original_all_gates<-do.call(rbind,list_original_dfs_all_classes)
  return(df_original_all_gates)
}

#' add_labels_column
#'
#' function to add label association column.
#' @param df Dataframe.
#' @param labels_assocation Vector of label association.
#' @return Dataframe.
#' @keywords flowMagic
#' @export
#' @examples 
#' \donttest{add_labels_column()}


add_labels_column<-function(df,labels_assocation){
  df$labels<-rep("none",nrow(df))
  for(l in labels_assocation){
    current_label<-strsplit(l,":")[[1]][1]
    current_class<-as.numeric(strsplit(l,":")[[1]][2])
    inds<-which(df[,3]==current_class)
    df$labels[inds]<-current_label
  }
  return(df)
}

#' update_label_association
#'
#' function to update label association.
#' @param df Dataframe.
#' @return Vector.
#' @keywords flowMagic
#' @export
#' @examples 
#' \donttest{update_label_association()}

update_label_association<-function(df){
  all_labels<-unique(df$labels)
  all_labels<-all_labels[all_labels!="none"]
  labels_association<-c()
  for(l in all_labels){
    inds<-which(df$labels==l)
    classes<-df$classes[inds]
    string<-sprintf("%s:%s",l,classes[1])
    labels_association<-c(labels_association,string)
  }
  return(labels_association)
}

#' range01
#'
#' function to put data in range 0-1.
#' @param x Vector of numbers to scale.
#' @return Vector of numbers.
#' @keywords flowMagic
#' @export
#' @examples 
#' \donttest{range01()}

range01 <- function(x){(x-min(x))/(max(x)-min(x))}

#' get_indices_cross_val
#'
#' function to get indices for cross val.
#' @param df_train Dataframe of training features generated by the get_train_data function.
#' @param n_cores Number of cores to use. Default to 1.
#' @param train_inds Type of method to extract training indices:  plot_num,rand_set_num,rand_set_n_gates_info.
#' @param val_inds Type of method to extract validation indices:  plot_num,rand_set_num,rand_set_n_gates_info.
#' @param n_train_plots Number of training plots to include in training data for each iteration.
#' @param n_folds Number of training iterations.
#' @param seed Seed to randomly extract data.
#' @param n_val_plots Number of training plots to include in validation data for each iteration.
#' @return List of integers.
#' @keywords flowMagic
#' @export
#' @examples 
#' \donttest{get_indices_cross_val()}

get_indices_cross_val<-function(df_train,n_cores=1,train_inds="plot_num",val_inds="none",n_train_plots=5,
                                n_folds=5,seed=40,n_val_plots=5){
  
  # get training indices
  print("----- get training indices")
  if(train_inds=="plot_num"){
    all_plot_num<-unique(df_train$plot_num)
    n_folds<-length(all_plot_num)
    list_inds_train<-parallel::mclapply(1:n_folds,function(i){
      current_plot_num<-all_plot_num[i]
      inds<-which(df_train$plot_num==current_plot_num)
      return(inds)
    },mc.cores = n_cores)
    all_names_list<-sapply(1:length(list_inds_train),function(i){
      name_i<-sprintf("fold_%s",i)
      return(name_i)
    })
    list_inds_train <- setNames(list_inds_train,all_names_list)
  }else if(train_inds=="rand_set_num"){
    list_inds_train<-list()
    all_plot_num<-unique(df_train$plot_num)
    for(f in 1:n_folds){
      print(sprintf("fold_%s",f))
      set.seed(seed+f)
      rand_ints<-sample.int(n = length(all_plot_num),size = n_train_plots)
      all_plot_num_selected<-all_plot_num[rand_ints]
      list_inds_train_temp<-parallel::mclapply(1:length(all_plot_num_selected),function(i){
        current_plot_num<-all_plot_num_selected[i]
        inds<-which(df_train$plot_num==current_plot_num)
        return(inds)
      },mc.cores = n_cores)
      inds_all_plot_num_selected<-unlist(list_inds_train_temp)
      list_inds_train[[sprintf("fold_%s",f)]]<-inds_all_plot_num_selected
    }
  }else if(train_inds=="rand_set_n_gates_info"){
    list_inds_train<-list()
    all_n_gates_info<-unique(df_train$n_gates_info)
    for(f in 1:n_folds){
      print(sprintf("fold_%s",f))
      set.seed(seed+f)
      list_inds_train_ngates<-list()
      for(n_gates in all_n_gates_info){
        print(sprintf("n_gates_info_%s",n_gates))
        df_train_current_ngates<-df_train[df_train$n_gates_info==n_gates,]
        all_plot_num<-unique(df_train_current_ngates$plot_num)
        if(length(all_plot_num)<n_train_plots){
          all_plot_num_selected<-all_plot_num
        }else{
          rand_ints<-sample.int(n = length(all_plot_num),size = n_train_plots)
          all_plot_num_selected<-all_plot_num[rand_ints]
        }
        list_inds_train_temp<-parallel::mclapply(1:length(all_plot_num_selected),function(i){
          current_plot_num<-all_plot_num_selected[i]
          inds<-which(df_train$plot_num==current_plot_num)
          return(inds)
        },mc.cores = n_cores)
        inds_all_plot_num_selected<-unlist(list_inds_train_temp)
        list_inds_train_ngates[[sprintf("n_gates_%s",n_gates)]]<-inds_all_plot_num_selected

      }
      inds_all_n_gates_selected<-unlist(list_inds_train_ngates)
      list_inds_train[[sprintf("fold_%s",f)]]<-inds_all_n_gates_selected
    }
  }else if(train_inds=="leave_one_out"){
    list_inds_train<-list()
    list_inds_val_one_out<-list()
    all_plot_num<-unique(df_train$plot_num)
    for (f in seq_along(all_plot_num)){
      # Identify the one group to leave out
      plot_num_out <- all_plot_num[f]
      
      # Validation indices: those with the left-out group
      inds_val_out <- which(df_train$plot_num == plot_num_out)
      
      # Training indices: all others
      inds_train <- which(df_train$plot_num != plot_num_out)
      
      # Store in list format that caret expects
      fold_name <- sprintf("fold_%s", f)
      list_inds_train[[fold_name]] <- inds_train
      list_inds_val_one_out[[fold_name]] <- inds_val_out
    }
  }
  #------- get validation indices
  print("----- get validation indices")
  if(val_inds=="plot_num"){
    all_plot_num<-unique(df_train$plot_num)
    n_folds<-length(all_plot_num)
    list_inds_val<-parallel::mclapply(1:n_folds,function(i){
      current_plot_num<-all_plot_num[i]
      inds<-which(df_train$plot_num==current_plot_num)
      return(inds)
    },mc.cores = n_cores)
    all_names_list<-sapply(1:length(list_inds_val),function(i){
      name_i<-sprintf("fold_%s",i)
      return(name_i)
    })
    list_inds_val <- setNames(list_inds_val,all_names_list)
  }else if(val_inds=="rand_set_num"){
    list_inds_val<-list()
    all_plot_num<-unique(df_train$plot_num)
    for(f in 1:n_folds){
      print(sprintf("fold_%s",f))
      set.seed(seed-f)
      rand_ints<-sample.int(n = length(all_plot_num),size = n_val_plots)
      all_plot_num_selected<-all_plot_num[rand_ints]
      list_inds_val_temp<-parallel::mclapply(1:length(all_plot_num_selected),function(i){
        current_plot_num<-all_plot_num_selected[i]
        inds<-which(df_train$plot_num==current_plot_num)
        return(inds)
      },mc.cores = n_cores)
      inds_all_plot_num_selected<-unlist(list_inds_val_temp)
      list_inds_val[[sprintf("fold_%s",f)]]<-inds_all_plot_num_selected
    }
  }else if(val_inds=="none"){
    list_inds_val<-NULL
  }else if(train_inds=="rand_set_n_gates_info"){
    list_inds_val<-list()
    all_n_gates_info<-unique(df_train$n_gates_info)
    for(f in 1:n_folds){
      print(sprintf("fold_%s",f))
      set.seed(seed-f)
      list_inds_val_ngates<-list()
      for(n_gates in all_n_gates_info){
        print(sprintf("n_gates_info_%s",n_gates))
        df_train_current_ngates<-df_train[df_train$n_gates_info==n_gates,]
        all_plot_num<-unique(df_train_current_ngates$plot_num)
        if(length(all_plot_num)<n_train_plots){
          all_plot_num_selected<-all_plot_num
        }else{
          rand_ints<-sample.int(n = length(all_plot_num),size = n_train_plots)
          all_plot_num_selected<-all_plot_num[rand_ints]
        }
        list_inds_val_temp<-parallel::mclapply(1:length(all_plot_num_selected),function(i){
          current_plot_num<-all_plot_num_selected[i]
          inds<-which(df_train$plot_num==current_plot_num)
          return(inds)
        },mc.cores = n_cores)
        inds_all_plot_num_selected<-unlist(list_inds_val_temp)
        list_inds_val_ngates[[sprintf("n_gates_%s",n_gates)]]<-inds_all_plot_num_selected
        
      }
      inds_all_n_gates_selected<-unlist(list_inds_val_ngates)
      list_inds_val[[sprintf("fold_%s",f)]]<-inds_all_n_gates_selected
    }
  }else if(val_inds=="leave_one_out"){
    list_inds_val<-list_inds_val_one_out
  }

  return(list(inds_train=list_inds_train,inds_val=list_inds_val))
}


#' process_test_data
#'
#' function to get test data correctly formatted.
#' @param test_data Dataframe of bivariate markers expression.
#' @param prop_down Proportion of events (downsampling). Default to NULL (downsampling using number of points).
#' @param n_points_per_plot Number of points for downsampling.
#' @param normalize_data If True, data is normalized to 0-1 range. Default to True.
#' @return Dataframe.
#' @keywords flowMagic
#' @export
#' @examples 
#' \donttest{process_test_data()}

process_test_data<-function(test_data,prop_down=NULL,n_points_per_plot=500,normalize_data=T){
  names_dens_features<-c("n_peaks_m1","h_peak_m1_1","pos_peak_m1_1","start_peak_m1_1","end_peak_m1_1",
                         "h_peak_m1_2","pos_peak_m1_2","start_peak_m1_2","end_peak_m1_2",
                         "h_peak_m1_3","pos_peak_m1_3","start_peak_m1_3","end_peak_m1_3",
                         "h_peak_m1_4","pos_peak_m1_4","start_peak_m1_4","end_peak_m1_4",
                         "n_peaks_m2","h_peak_m2_1","pos_peak_m2_1","start_peak_m2_1","end_peak_m2_1",
                         "h_peak_m2_2","pos_peak_m2_2","start_peak_m2_2","end_peak_m2_2",
                         "h_peak_m2_3","pos_peak_m2_3","start_peak_m2_3","end_peak_m2_3",
                         "h_peak_m2_4","pos_peak_m2_4","start_peak_m2_4","end_peak_m2_4")
  colnames(test_data)<-c("x1_expr","x2_expr")
  test_data$classes<-rep("0",nrow(test_data))

  # downsample
  if(is.null(prop_down)==T){
      if(n_points_per_plot>nrow(test_data)){
       prop_down<-1
      }else{
       prop_down<-(n_points_per_plot/nrow(test_data))
      }
  }
  out_part<-caret::createDataPartition(y=factor(test_data[,"classes"]),times = 1,p = prop_down)
  inds_new_df<-out_part$Resample1
  Xtest<-test_data[inds_new_df,c(1,2)]

  # get density features
  if(normalize_data==T){
    Xtest[,1]<-range01(Xtest[,1])
    Xtest[,2]<-range01(Xtest[,2])
  }
  Xtest[,1]<-round(Xtest[,1],2)
  Xtest[,2]<-round(Xtest[,2],2)
  df_dens<-csv_to_dens(df = Xtest,with_classes = F,n_coord = 50,normalize_data = normalize_data)
  if(normalize_data==F){
    vec_info_dens<-get_density_features(df_dens = df_dens,min_height = 0.00)
  }else{
    vec_info_dens<-get_density_features(df_dens = df_dens)
  }
  m_info_dens<-matrix(vec_info_dens,length(vec_info_dens),nrow(Xtest))
  m_info_dens<-t(m_info_dens)
  df_info_dens<-as.data.frame(m_info_dens)
  colnames(df_info_dens)<-names_dens_features
  Xtest<-cbind(Xtest,df_info_dens)
  return(Xtest)
}

#' get_list_df_gated_plots
#'
#' function to get labelled dataframe based on selected node for each sample in gs.
#' @param gs GatingSet
#' @param gate_name Name of the Gating tree node whose gating data needs to be extracted.
#' @return List.
#' @keywords flowMagic
#' @export
#' @examples 
#' \donttest{get_list_df_gated_plots()}

get_list_df_gated_plots<-function(gs,gate_name){
  list_gated_data <- list()
  # Loop over all samples in the GatingSet
  for(s in sampleNames(gs)){
    # Get data from the specified gated population
    gh<-gs[[s]]
    
    # df_pop_root = mapping of selected pop events to root
    df_pop_root<-map_to_root(gh = gh,pop = gate_name)
    
    # df_pop_mother = mapping of selected pop events to mother of selected pop

    df_pop_mother<-map_to_parent(gh = gh,binary_df = df_pop_root)
    
    # Store in list
    list_gated_data[[s]] <- df_pop_mother
  }
  return(list_gated_data)
}

#' get_flowframe_from_gs
#'
#' function to get flowframe from gs (converting cytoframe to flowframe)
#' @param gs GatingSet
#' @param node_name Name of the Gating tree node whose gating data needs to be extracted.
#' @param sample_id Name or index of the sample to extract.
#' @return List.
#' @keywords flowMagic
#' @export
#' @examples 
#' \donttest{get_flowframe_from_gs()}

get_flowframe_from_gs<-function(gs,node_name,sample_id){
  # Get flowFrame
  ff <- flowWorkspace::gh_pop_get_data(gs[[sample_id]], node_name)
  # ff is a cytoframe object not compatible with flowDensity
  # we need to build a flowFrame from the cytoframe
  # Note: as(ff, "flowFrame") does not work because the cytoframe 
  # is only a visual representation of data not real data
  
  # Extract expression data
  # exprs_mat <- flowCore::exprs(ff)  # or exprs(ff_temp)
  
  # Extract parameter data
  # param_data <- flowCore::pData(parameters(ff)) 
  
  # Extract metadata/keywords
  # desc <- flowCore::keyword(ff)
  
  # ff <- flowFrame(exprs = exprs_mat, parameters = Biobase::AnnotatedDataFrame(param_data), description = desc)
  
  ff<-cytoframe_to_flowFrame(cf = ff)

  return(ff)
}

#' magic_label_rectangle
#'
#' function to label points of dataframe based on rectangle coordinates
#' @param df Dataframe composed of two columns for marker expression of first (x axis = first column) and second marker (y axis = second column)
#' @param x_min x coordinate minimum
#' @param x_max  x coordinate maximum
#' @param y_min y coordinate minimum
#' @param y_max  y coordinate maximum
#' @param label_pol  Label for current polygon.
#' @return Dataframe
#' @keywords flowMagic
#' @export
#' @examples 
#' \donttest{magic_label_rectangle()}

magic_label_rectangle<-function(df,x_min,x_max,y_min,y_max,label_pol="1"){
  # Initialize class column if it doesn't exist
  if (!"class" %in% names(df)) {
    df$class <- "0"  
  }
  df$class[df[,1] > x_min & df[,1] < x_max & df[,2] > y_min & df[,2]] <- label_pol
  return(df)
}

#' magic_label_poly
#'
#' function to label points of dataframe based on polygon coordinates
#' @param df Dataframe composed of two columns for marker expression of first (x axis = first column) and second marker (y axis = second column)
#' @param polygon_df Dataframe containing the x coordinates (first column) and y coordinates (second column) of the current polygon to label.
#' @param label_pol  Label for current polygon.
#' @return Dataframe
#' @keywords flowMagic
#' @export
#' @examples 
#' \donttest{magic_label_poly()}

magic_label_poly<-function(df,polygon_df,label_pol="1"){
  # Initialize class column if it doesn't exist
  if (!"class" %in% names(df)) {
    df$class <- "0"  
  }
  inside <- sp::point.in.polygon(df[,1], df[,2], polygon_df[,1], polygon_df[,2])
  df$class[inside == 1] <- label_pol
  return(df)
}

#' flowmagic_pred_to_poly_gates
#'
#' function to label points of dataframe based on polygon coordinates
#' @param list_df List of dataframes with three columns: marker 1, marker 2, class. Can get input generated by the magicPred function.
#' @param pred_label Select prediction label to extract polygon.
#' @param gate_label  Label for the selected polygon within the GatingSet framework.
#' @param n_cores  Number of cores to use. Default to 1.
#' @param concavity_val  Concavity value for drawing polygon. Default to 10.
#' @return list
#' @keywords flowMagic
#' @export
#' @examples 
#' \donttest{flowmagic_pred_to_poly_gates()}

flowmagic_pred_to_poly_gates<-function(list_df,pred_label,gate_label,n_cores=1,concavity_val=10){
  start <- Sys.time()
  all_names<-names(list_df)
  message("Get selected polygon gate...")
  list_poly_gates<-parallel::mclapply(seq_along(list_df),function(i){
    current_sample<-all_names[i]
    message(current_sample)
    if ("df_test_original" %in% names(list_df[[i]])) {
      df <- list_df[[i]]$df_test_original
    }else{
      df <- list_df[[i]]
    }
    names_columns<-colnames(df)
    list_coords<-flowMagic::extract_polygon_gates(gated_df = df,concavity_val = concavity_val)
    df_coord_label_selected<-list_coords[[pred_label]]
    coords<-as.matrix(df_coord_label_selected[,c(1,2)])
    colnames(coords)<-colnames(df[,c(1,2)])
    pg <- flowCore::polygonGate(filterId = gate_label,
                      .gate = coords,
                      parameters = c(names_columns[1], names_columns[2]))
    return(pg)
  },mc.cores=n_cores)
  names(list_poly_gates)<-all_names
  end <- Sys.time()
  time_taken <- end - start
  print("Execution time:")
  print(time_taken)
  print("Done")
  return(list_poly_gates)
}

#' flowmagic_pred_to_gs
#'
#' function to label points of dataframe based on polygon coordinates
#' @param list_poly_gates List of polygonGate objects generated using the flowmagic_pred_to_poly_gates
#' @param gs GatingSet to update with the flowMagic polygons.
#' @param parent_node  Parent population name of the new gate.
#' @return GatingSet
#' @keywords flowMagic
#' @export
#' @examples 
#' \donttest{flowmagic_pred_to_gs()}

flowmagic_pred_to_gs<-function(list_poly_gates,gs,parent_node){
  all_names<-names(list_poly_gates)
  for(i in seq_along(list_poly_gates)){
    current_sample<-all_names[i]
    message(current_sample)
    current_poly<-list_poly_gates[[i]]
    print(current_poly)
    suppressWarnings({
      flowWorkspace::gs_pop_add(gs[[current_sample]], current_poly, parent = parent_node, name = current_poly@filterId)
    })
    
  }
  recompute(gs)
  return(gs)
}

#' convert_to_integers_chr
#'
#' function to text classes to integers-like classes if necessary.
#' @param df Dataframe with classes (third column) to check and eventually modify
#' @return Dataframe
#' @keywords flowMagic
#' @export
#' @examples 
#' \donttest{convert_to_integers_chr()}

convert_to_integers_chr<-function(df){
  classes<-df[,3]
  has_real_text_labels <- any(!grepl("^[0-9]$", classes))
  if(has_real_text_labels){
    message("converting to integers-like characters")
    classes_vec<-df$classes
    # Create a logical index for non-zero (i.e., character) labels
    is_label <- classes_vec != "0"

    # Get unique labels (excluding "0")
    unique_labels <- unique(classes_vec[is_label])

    # Create a named mapping: text label -> integer (starting from 1)
    label_map <- stats::setNames(seq_along(unique_labels), unique_labels)

    classes_vec[is_label] <- label_map[classes_vec[is_label]]

    df$classes<-classes_vec
    }
  return(df)
}

#' magicGating
#'
#' function to manually gate samples in a flowSet.
#' @param fs An object of class flowSet. Can be also an object of class cytoset (it will be converted to flowSet) or a GatingSet object (in this case, gs_node is mandatory)
#' @param sample_id Names of the samples to gate. It can also be the numerical index of the sample. Default to sample 1.
#' @param channel_x Name of the channel (x-axis).
#' @param channel_y Name of the channel (y-axis).
#' @param gs_node Name of the node to extract data, if fs is a GatingSet object this a mandatory argument.

#' @param label_pol Label of the gate polygon. Default to "1".
#' @return List of two objects of class List. 
#'         list_poly_gates: List of polygon coordinates. Each element refers to the coordinates of one sample. 
#'         list_gated_data: List of gated data. Each element refers to the gated data of one sample. 
#' @keywords flowMagic
#' @export
#' @examples 
#' \donttest{magicGating()}

magicGating<-function(fs,sample_id=1,channel_x,channel_y, gs_node=NULL, label_pol="1",...){

  # capture all extra arguments
  args_list <- list(...)
  
  # Arguments for magicPlot_template
  args_plot <- args_list[names(args_list) %in% names(formals(magicPlot_template))]
  
  # Arguments for flowmagic_pred_to_poly_gates
  args_gates <- args_list[names(args_list) %in% names(formals(flowmagic_pred_to_poly_gates))]

  if(class(fs)=="cytoset"){
    fs<-cytoset_to_flowSet(cs = fs)
  }
  if(class(fs)=="GatingSet"){
    if(is.null(gs_node)==T){
      stop("fs is a GatingSet object, please provide a valid node name in the gs_node argument")
    }
    fs<-gh_pop_get_data(obj = gs,y = gs_node)
    fs<-cytoset_to_flowSet(cs = fs)
  }
  if(class(fs)=="cytoframe" || class(fs)=="flowFrame"){
    stop("fs need to be a flowSet,cytoset or GatingSet object")
  }
  if(length(sample_id)==1){
    sample_id<-c(sample_id)
  }
  list_gated_data<-rep(list(0), length(sample_id))
  names(list_gated_data)<-as.character(sample_id)
  for(id in sample_id){
    print("--- Gating sample with id:")
    print(id)
    ff <- fs[[id]]
    expr_matrix_ff <- exprs(ff)
    df_exprs <- as.data.frame(expr_matrix_ff)
    df_exprs_selected_channels <- df_exprs[, c(channel_x, channel_y)]
    
    # get gate coordinates
    print("Get gate coordinates")
    
    polygon_df <- do.call(magicPlot_template, c(list(df = df_exprs_selected_channels), args_plot))
    
    df_gated<-flowMagic::magic_label_poly(df = df_exprs_selected_channels,polygon_df = polygon_df,label_pol = label_pol)
    # add gated data to list
    print("Add gated data to list")
    list_gated_data[[as.character(id)]]<-df_gated
    

  }
  
  # add gates coordinates to list of polygon gates
  print("Add gates coordinates to list")
  list_poly_gates <- do.call(flowmagic_pred_to_poly_gates,
                             c(list(list_df = list_gated_data,
                                    pred_label = label_pol,
                                    gate_label = label_pol),
                               args_gates))
  
  return(list(list_poly_gates=list_poly_gates,list_gated_data=list_gated_data))
}

