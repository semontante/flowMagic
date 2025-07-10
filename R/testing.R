

#' magicPred_hierarchy
#' 
#' function to predict the gates on the ungated .fcs samples.
#' @param list_test_sets contains the list of root dataframe for each ungated fcs file imported.
#' @param list_models_local contains the optimized local models pre-generated using the magicTrain_local function.
#' @param df_tree contains the info related to the populations hierarchy.
#' @param n_cores Number of cores to use. Default to 1.
#' @return List of Dataframes.
#' @export
#' @examples 
#' \donttest{magicPred_hierarchy()}


magicPred_hierarchy<-function(list_test_sets,list_models_local,df_tree,n_cores=1){
  n_samples<-length(list_test_sets)
  name_samples<-names(list_test_sets)
  start<-Sys.time()
  list_gated_data<-parallel::mclapply(1:n_samples,function(i){
    sample_name<-name_samples[i]
    sample_name<-strsplit(sample_name,".fcs")[[1]]
    print(sprintf("########## gating sample %s #########",sample_name))
    root_sample_i<-list_test_sets[[i]]
    all_levels<-names(list_models_local)
    list_gated_data_sample_i<-list() # gated data with the correct dims
    list_temp_gated_data_all_dims<-list() # gated data with all dims
    for(level in all_levels){
      print(sprintf("########## gating pops in level %s #########",level))
      current_level<-list_models_local[[level]]
      all_pops_current_level<-names(current_level)
      list_gated_data_level<-list()
      list_temp_gated_data_all_dims_level<-list()
      for(pop in all_pops_current_level){
        print(sprintf("########## gating pop %s #########",pop))
        model_current_pops_to_gate<-current_level[[pop]]
        # check if pop to gate is binary or multiclass
        stringsplitted<-strsplit(pop,"_")[[1]]
        if(length(stringsplitted)>1){
          print("multiclass pop")
          # check hierarchy position of the pops
          logic<-df_tree$Children %in% stringsplitted
          inds<-which(logic==T)
          mother_current_pop<-df_tree$Mother[inds]
          mother_current_pop<-unique(mother_current_pop)
        }else if(length(stringsplitted)==1){
          print("binary pop")
          # check hierarchy position of this pop
          logic<-df_tree$Children %in% stringsplitted
          ind<-which(logic==T)
          mother_current_pop<-df_tree$Mother[ind]
        }
        print(sprintf("mother current pop to gate:%s",mother_current_pop))
        # get dimensions current pop (or pops) to gate
        dims<-model_current_pops_to_gate[["Dimensions_set"]]
        # get labels association current pop to gate
        label_set<-model_current_pops_to_gate[["Labels_set"]]
        print(sprintf("Dimensions_set:%s",paste0(dims,collapse = " ")))
        print(sprintf("Labels_set: %s",paste0(label_set,collapse = " ")))
        ########################## get test df #######################
        # get appropriate test df based on the current hierarchy level
        if(mother_current_pop=="root"){
          print(sprintf("the mother of %s is root",pop))
          print("---- get root df with correct dims ----")
          Xtest<-root_sample_i[,dims]
        }else if(mother_current_pop!="root"){
          print(sprintf("the mother of %s is %s",pop,mother_current_pop))
          ########## look for the correct dataset in the already gated pops list with all dims#######
          print("----- get associated test mother df with correct dims ------")
          out<-get_slot_hierarchy_list(list_temp_gated_data_all_dims,mother_current_pop)
          df_containing_mother<-out[[1]]
          if(is.character(df_containing_mother)==T){
            print("Mother pop is None. Hierarchy is broken. No Xtest")
            Xtest<-"none_mother"
          }else{
            label_assoc<-out[[2]]
            print("all labels association selected df")
            label_assoc<-label_assoc[[1]]
            print(label_assoc)
            ind<-grep(mother_current_pop,label_assoc,fixed=T)
            label_assoc<-label_assoc[ind]
            print("selected label association:")
            print(label_assoc)
            s<-strsplit(label_assoc,":")[[1]]
            label_mother_pop<-s[2]
            print(sprintf("the associated label of the mother,%s,is:%s",mother_current_pop,label_mother_pop))
            all_labels_df<-df_containing_mother[,ncol(df_containing_mother)]
            all_dims_df<-df_containing_mother[,-ncol(df_containing_mother)]
            inds<-which(all_labels_df==label_mother_pop)
            df_only_mother_all_dims<-all_dims_df[inds,]
            df_only_mother<-df_only_mother_all_dims[,dims]
            Xtest<-df_only_mother
          }
        }
        ######################### get local df #######################
        # print("---- get local train data---- ")
        # local_train_df<-list_local_train[[level]][[pop]]
        ########################## check number of points in test mother df ##############################
        print("######### check number of points in test mother df (Xtest) ######### ")
        n_points<-nrow(Xtest)
        if(is.null(n_points)==T){
          warning("mother is None. Hierarchy is broken for this sample.")
          check_1<-F
          check_2<-F
          ####### report results in lists ###########
          print("------ report results ----------")
          gated_data_current_pop_all_dims<-"None"
          gated_data_current_pop<-"None"
          polygons_coords_list<-"None"
          df_gates_indices_on_root<-"None"
          type_gate<-"Hierarchy_broken"
          list_gated_data_level[[pop]]<-list(gated_data_current_pop,label_set,polygons_coords_list,df_gates_indices_on_root,type_gate)
          names(list_gated_data_level[[pop]])<-c("gated_data","labels","poly_coordinates","inds_on_root","info")
          list_temp_gated_data_all_dims_level[[pop]]<-list(gated_data_current_pop_all_dims,label_set)
        }else if(n_points<3){ # not enough points to make a gate (min 2)
          warning("mother df contains less than 3 points. No gate calculable")
          check_1<-F
          check_2<-F
          ####### report results in lists ###########
          print("------ report results ----------")
          gated_data_current_pop_all_dims<-"None"
          gated_data_current_pop<-"None"
          polygons_coords_list<-"None"
          df_gates_indices_on_root<-"None"
          type_gate<-"None"
          list_gated_data_level[[pop]]<-list(gated_data_current_pop,label_set,polygons_coords_list,df_gates_indices_on_root,type_gate)
          names(list_gated_data_level[[pop]])<-c("gated_data","labels","poly_coordinates","inds_on_root","info")
          list_temp_gated_data_all_dims_level[[pop]]<-list(gated_data_current_pop_all_dims,label_set)
        }else{ # enough points to make a gate
          ####### predict the gates ######
          reference_model_local<-model_current_pops_to_gate$ref_model_info
          print("----- predicting gates -----")
          out_pred<-magicPred(test_data = Xtest,ref_model_info=reference_model_local,n_cores=1,
                              prop_down=0.9,thr_dist=0.05,normalize_data=F)
          final_df<-out_pred$test_data_original
          ######## generate gated df current pop ######################
          print("------ generating gated data current pop -----")
          yhat_final_current_pop<-final_df[,3]
          gated_data_current_pop<-cbind(Xtest,yhat_final_current_pop)
          # get indices of the gates based on the root df
          indices_pop_on_root<-row.names(gated_data_current_pop)
          
          df_gates_indices_on_root<-as.data.frame(cbind(indices_pop_on_root,yhat_final_current_pop))
          df_gates_indices_on_root$indices_pop_on_root<-as.character(df_gates_indices_on_root$indices_pop_on_root)
          ####### report results in lists ###########
          print("------ report results ----------")
          if(mother_current_pop=="root"){
            gated_data_current_pop_all_dims<-cbind(root_sample_i,yhat_final_current_pop)
          }else{
            gated_data_current_pop_all_dims<-cbind(df_only_mother_all_dims,yhat_final_current_pop)
          }
          list_gated_data_level[[pop]]<-list(gated_data_current_pop,label_set,df_gates_indices_on_root,final_df)
          names(list_gated_data_level[[pop]])<-c("gated_data","labels","inds_on_root","final_df")
          list_temp_gated_data_all_dims_level[[pop]]<-list(gated_data_current_pop_all_dims,label_set)
          
        }
      }
      list_temp_gated_data_all_dims[[level]]<-list_temp_gated_data_all_dims_level
      list_gated_data_sample_i[[level]]<-list_gated_data_level
    }
    return(list_gated_data_sample_i)
  },mc.cores = n_cores)
  gc()
  names(list_gated_data)<-name_samples
  end<-Sys.time()
  time_taken<-end-start
  print("prediction time:")
  print(time_taken)
  return(list_gated_data)
}

#' magicPred
#' 
#' function to predict on plain test data (no hierarchy)
#' @param test_data Dataframe of test data to gate. It has only the two columns of marker expression.
#' @param magic_model Global trained model to predict gates. It can be a single model or list of named models (each model trained on selected number of gates).
#' @param magic_model_n_gates  Global trained model to predict number of gates. If different from NULL, magic_model is expected to be a list of models to predict certain gates (e.g., 5 models for 2,3,4,5 or 6 gates).
#' @param ref_model_info Template model to predic gates.
#' @param n_cores Number of cores to use. Default to 1.
#' @param ref_data_train Template data used to generate ref_model_info. Needed to calculate target-template distance.
#' @param prop_down Proportion for downsampling. Default to NULL (automatic downsampling using n_points_per_plot).
#' @param n_points_per_plot Number of points to consider for downsampling. Default to 500.
#' @param normalize_data If True, data is normalized to 0-1 range. Default to True.
#' @param include_zero_val considering events labeled as 0 as an additional gate when there is only one gate. Default to  True.
#' @return List of Dataframes.
#' @export
#' @examples 
#' \donttest{magicPred()}


magicPred<-function(test_data,magic_model=NULL,magic_model_n_gates=NULL,ref_model_info=NULL,n_cores=1,ref_data_train=NULL,
                    prop_down=NULL,thr_dist=0.05,n_points_per_plot=NULL,normalize_data=T,include_zero_val=T){
  set.seed(40)
  start<-Sys.time()
  if(ncol(test_data)>2){
    stop("only 2 columns must be present in data to gate")
  }
  # --------- prepare test data -------------
  message("----- prepare test data -------")
  if(is.null(ref_model_info)==F && is.null(prop_down)==T && is.null(n_points_per_plot)==T){
    prop_down<-1
  }else if(is.null(ref_model_info)==T  && is.null(prop_down)==T && is.null(n_points_per_plot)==T){
    n_points_per_plot<-500
  }
  Xtest<-process_test_data(test_data = test_data,prop_down = prop_down,n_points_per_plot = n_points_per_plot,
                           normalize_data = normalize_data)
  #show(magicPlot(Xtest[,c(1,2)],type = "no_gate",size_points = 2))
  #---------- get predictions based on provided model ------------
  message("---------- get predictions based on provided model ------------")
  if(is.null(magic_model)==F && is.null(ref_model_info)==T){
    if(is.null(magic_model_n_gates)==F){
      if(is.list(magic_model)==F){
        stop("List of magic models required.")
      }
      message("Using list of general models and n_gates model")
      #---- list of models and n_gates model
      if(class(magic_model_n_gates)=="train"){
        n_gates<-caret::predict.train(magic_model_n_gates,Xtest)
        inds_max_gates<-which.max(table(n_gates))
        max_gate<-names(table(n_gates))[inds_max_gates]
        message("Number of gates predicted")
        message(max_gate)
      }else{
        max_gate<-as.character(magic_model_n_gates)
        message("Number of gates selected")
        message(max_gate)
      }
      magic_model_selected<-magic_model[[max_gate]]
      classes<-caret::predict.train(magic_model_selected,Xtest)
      vec_dist<-0
    }else if(is.null(magic_model_n_gates)==T){
      message("Using general purpose model")
      #---- only  PD models predictions
      classes<-caret::predict.train(magic_model,Xtest)
      # get distance templates - test data
      vec_dist<-0
    }
    
  }else if(is.null(ref_model_info)==F){
    message("Using template model")
    #---- only reference predictions
    classes<-caret::predict.train(ref_model_info,Xtest)
    # get distance templates - test data
    if(is.null(ref_data_train)==F){
      all_plot_num<-unique(ref_data_train$plot_num)
      list_scores_ref<-lapply(1:length(all_plot_num),function(i){
        ref_data_train_plot_num_i<-ref_data_train[which(ref_data_train$plot_num==all_plot_num[i]),]
        check_col<-colnames(ref_data_train_plot_num_i) %in% c("x1_expr","x2_expr","classes","plot_num","n_gates_info")
        inds<-which(check_col==T)
        vec_scores<-ref_data_train_plot_num_i[1,-inds]
        return(vec_scores)
      })
      
      df_scores_ref<-do.call(rbind,list_scores_ref)
      df_scores<-rbind(df_scores_ref,Xtest[1,-c(1,2)])
      vec_dist<-get_weights_density_features(df_scores = df_scores)
    }else{
      vec_dist<-0
    }
  }else{
    stop("InputError: Either general purpose model or template model must be provided")
  }
  test_data_final<-cbind(Xtest,classes)
  final_df<-test_data_final[,c("x1_expr","x2_expr","classes")]
  final_df$classes<-as.character(final_df$classes)
  #------------------- post-processing ----------------
  message("---------- post-processing ------------")
  all_classes<-unique(final_df$classes)
  all_classes<-all_classes[all_classes!="0"]
  if(length(all_classes)!=0){
    if(is.null(ref_model_info)==T){
      # ------- no template model ---------
      if(length(all_classes)==1){
        final_df<-post_process_gates(gated_df=final_df,n_cores=n_cores,include_zero = include_zero_val,thr_dist = thr_dist,
                                       type="dist")
        
      }else{
        final_df<-post_process_gates(gated_df=final_df,n_cores=n_cores,include_zero = F,thr_dist = thr_dist,type="dist")
      }
    }else{
      # ------- Yes template model ---------
      if(length(all_classes)<=4){
      message("Yes template model with less than 4 polygons")
        final_df<-post_process_gates(gated_df=final_df,n_cores=n_cores,type = "polygon",normalize_data = normalize_data)
      }else{
      message("Yes template model with more than 4 polygons")
        final_df<-post_process_gates(gated_df=final_df,n_cores=n_cores,include_zero = F,thr_dist = thr_dist,type="dist")
      }
    }
  }
  
  # get polygons after post-processing
  message("------ get polygons after post-processing ------")
  list_df_hull<-extract_polygon_gates(gated_df = final_df,concavity_val=5)
  
  # compute gates on original data
  message("------ compute gates ------")
  test_data_temp<-test_data
  test_data_temp$classes<-rep("0",nrow(test_data_temp))
  if(is.null(list_df_hull)==F){
    if(normalize_data==T){
      test_data_temp[,1]<-range01(test_data_temp[,1])
      test_data_temp[,2]<-range01(test_data_temp[,2])
    }
    test_data_temp[,1]<-round(test_data_temp[,1],2)
    test_data_temp[,2]<-round(test_data_temp[,2],2)
    test_data_temp_original<-compute_gates(gated_df=test_data_temp,list_final_polygons_coords =  list_df_hull)
    test_data_temp_original$classes<-as.character(test_data_temp_original$classes)
    test_data_original<-cbind(test_data,test_data_temp_original$classes)
  }else{
    test_data_temp_original<-test_data_temp
    test_data_original<-cbind(test_data,test_data_temp_original$classes)
  }
  
  # return gated results
  vec_dist<-round(vec_dist,2)
  end<-Sys.time()
  time_taken<-end-start
  message("Execution time:")
  message(time_taken)
  message("Done")
  return(list(test_data_original=test_data_original,test_data_temp_original=test_data_temp_original,
              final_df=final_df,vec_dist=vec_dist,test_data_final=test_data_final))
}


#' magicPred_all
#' 
#' function to predict on plain test data (no hierarchy)
#' @param list_test_data List of unlabeled  dataframes. It has only the two columns of marker expression.
#' @param magic_model Global trained model to predict gates. It can be a single model or list of named models (each model trained on selected number of gates).
#' @param magic_model_n_gates  Global trained model to predict number of gates. If different from NULL, magic_model is expected to be a list of models to predict certain gates (e.g., 5 models for 2,3,4,5 or 6 gates).
#' @param ref_model_info Template model to predic gates.
#' @param n_cores Number of cores to use to process one sample. Default to 1.
#' @param ref_data_train Template data used to generate ref_model_info. Needed to calculate target-template distance.
#' @param n_points_per_plot Number of points to consider for downsampling. Default to 500.
#' @param normalize_data If True, data is normalized to 0-1 range. Default to True.
#' @param include_zero_val considering events labeled as 0 as an additional gate when there is only one gate. Default to  True.
#' @param n_cores_all Number of cores to use across all samples. Default to 1.
#' @param verbose If True, print all message and disable tryCatch (any error will stop the execution). Default to False.
#' @return List of Dataframes.
#' @export
#' @examples 
#' \donttest{magicPred_all()}

magicPred_all<-function(list_test_data,magic_model=NULL,ref_model_info=NULL,magic_model_n_gates=NULL,
                        ref_data_train=NULL,prop_down=NULL,n_points_per_plot=NULL,
                        thr_dist=0.05,n_cores=1,normalize_data=T,include_zero_val=T,n_cores_all=1,verbose=F){
  if (("package:dplyr" %in% search())==F) {
  library(dplyr)
  }                      
  set.seed(40)
  start<-Sys.time()
  all_names_test_data<-names(list_test_data)
  # prediction for each test data
  message("------------- Prediction for each test data")
  list_all_dfs_pred<-parallel::mclapply(1:length(list_test_data),function(i){
    message(sprintf("########### %s ##########",all_names_test_data[i]))
    df_test<-list_test_data[[i]]
    if(verbose==T){
    out_pred<-magicPred(test_data = df_test,magic_model=magic_model,
                                            ref_model_info=ref_model_info,n_cores=n_cores,
                                            ref_data_train=ref_data_train,prop_down=prop_down,thr_dist=thr_dist,
                                            magic_model_n_gates = magic_model_n_gates,n_points_per_plot=n_points_per_plot,
                                            normalize_data=normalize_data,include_zero_val=include_zero_val)
    }else{
    out_pred<-tryCatch(suppressMessages(magicPred(test_data = df_test,magic_model=magic_model,
                                                  ref_model_info=ref_model_info,n_cores=n_cores,
                                                  ref_data_train=ref_data_train,prop_down=prop_down,thr_dist=thr_dist,
                                                  magic_model_n_gates = magic_model_n_gates,n_points_per_plot=n_points_per_plot,
                                                  normalize_data=normalize_data,include_zero_val=include_zero_val)),error=function(e){return(NULL)})
    }

    
    
    if(is.null(out_pred)==F){
      df_test_original<-out_pred$test_data_original
      final_df<-out_pred$final_df
      vec_dist<-out_pred$vec_dist
    }else{
      df_test_original<-NULL
      final_df<-NULL
      vec_dist<-NULL
    }
    
    return(list(df_test_original=df_test_original,final_df=final_df,vec_dist=vec_dist))
  },mc.cores = n_cores_all)
  names(list_all_dfs_pred)<-all_names_test_data
  end<-Sys.time()
  time_taken<-end-start
  print("Execution time:")
  print(time_taken)
  print("Done")
  return(list_all_dfs_pred)
}

