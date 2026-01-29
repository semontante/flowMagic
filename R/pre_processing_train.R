#' name_pop_gating
#' 
#' function to get the name of all the pops of the gating hierarchy.
#' @param gh GatingHierarchy.
#' @return Vector of characters.
#' @keywords flowMagic
#' @export
#' @examples 
#' \donttest{name_pop_gating()}


name_pop_gating<-function(gh){
  sample_name<-sampleNames(gh)
  nodes_list<-gs_get_pop_paths(gh)
  vector_pops<-sapply(1:length(nodes_list), function(i){
    splitted_nodes<-strsplit(gs_get_pop_paths(gh)[i],"/")[[1]]
    last_node<-splitted_nodes[length(splitted_nodes)]
    return(last_node)
  })
  return(vector_pops)
}

#' map_to_root
#' 
#' it generates a dataset indicating what cells belong to the selected pop (1) 0 otherwise. 
#' Old name: pre_process_manual_binary()
#' The dataset is always the dataset of the Root pop.
#' gh and pop are mandatory arguments, dim and mode have a default value.
#' @param gh GatingHierarchy.
#' @param pop Name of the population to get events assignments.
#' @return Dataframe.
#' @keywords flowMagic
#' @export
#' @examples 
#' \donttest{map_to_root()}


map_to_root<-function(gh,pop){
  ######### verify that the input is correct ######
  print("#### veryfing input... ####")
  if(is.null(gh) & is.null(pop)){
    stop("Input Error: gh and pop are mandatory arguments")
  }
  if(class(gh)!="GatingHierarchy"){
    stop("gh must be a GatingHierarchy object")
  }else if(class(gh)=="GatingHierarchy"){
    gh<-gh
  }
  print("Input correct")
  # Name of the manually gated sample, to which the GatingHierarchy is associated.
  sample_name<-sampleNames(gh)
  print(sample_name)
  #################### get the expression matrix of the root population ###################
  print("####### get expression matrix of the root population #########")
  f<-gh_pop_get_data(gh,"root") # It takes automatically the flowFrame associated to the sample of the GatingHierarchy
  exprs_matrix<-exprs(f)
  dataframe_expr<-as.data.frame(exprs_matrix) # expression matrix of the root population.
  ######################## make the binary dataset for the selected node  #######################
  # Now we make the binary dataset for the selected pop (events with 1 belongs to the pop, otherwise 0)
  #  Note: this assignation refers to the root (i.e., all events considered)
  print("###### generating the binary dataset for the selected pop ########")
  # ---- get indices of the selected pop of the GatingHierarchy
  print("get indices selected pop")
  logic<-gh_pop_get_indices(gh,pop)
  inds_selected_pop<-which(logic==T)
  new_column<-rep(0,nrow(dataframe_expr))
  new_column[inds_selected_pop]<-1
  new_dataframe_expr<-cbind(dataframe_expr,new_column)
  #------ check for correctness of the operations
  print("check operations")
  result_check<-identical(which(new_dataframe_expr$new_column==1),inds_selected_pop)
  if(result_check==FALSE){
    stop("Error detected: Wrong operations during training set preparation
         which(new_dataframe_expr$new_column==1) and list_indices[[ind_pop]]$indices not identical objects")
  }
  # rename gate label column
  print("rename gate column")
  ind<-which(colnames(new_dataframe_expr)=="new_column")
  colnames(new_dataframe_expr)[ind]<-pop # expression dataframe with assigned gates with all dimensions
  if(pop!="root"){ # root is not a gate, so the dimensions are not important
    #--------------- we find the dimensions of the selected pop --------------
    print("------ find dimensions (i.e., the markers name) of the selected population ------")
    polgate<-gh_pop_get_gate(gh,pop) # get gate object of the selected pop
    print(sprintf("finding dimensions of %s",polgate@filterId))
    if(class(polgate)=="polygonGate"){
      dims<-colnames(polgate@boundaries)
    }else if(class(polgate)=="rectangleGate"){
      rectGate<-polgate
      vec<-rectGate@min
      dims<-names(vec)
      # in case there is only one dimension indicated in the gh,as second dimension we use by default SSC_A
      if(length(dims)==1){
        dims<-c(dims,"SSC-A")
      }
    }else if(class(polgate)=="booleanFilter"){
      # I use the same dimensions of the parent.
      path_to_parent<-gs_pop_get_parent(gh_sample_gated,pop)
      splitted_path<-strsplit(path_to_parent,"/")[[1]]
      name_parent<-tail(splitted_path,1)
      polgate_parent<-gh_pop_get_gate(gh,name_parent) # get gate object of the selected pop
      if(class(polgate_parent)=="polygonGate"){
        dims<-colnames(polgate_parent@boundaries)
      }else if(class(polgate_parent)=="rectangleGate"){
        rectGate<-polgate_parent
        vec<-rectGate@min
        dims<-names(vec)
        # in case there is only one dimension indicated in the gh,as second dimension we use by default SSC_A
        if(length(dims)==1){
          dims<-c(dims,"SSC-A")
        }
      }else{
        stop("Dimensions of the parent of the boolean gate is still a boolean pop")
      }
    }else{
      print(polgate)
      stop("polgate class not recognized")
    }
    print(dims)
    print(sprintf("Dimensions found: %s",paste0(dims,collapse = ",")))
    new_dataframe_expr<-new_dataframe_expr[,c(dims,pop)]
  }
  print("Done")
  return(new_dataframe_expr)
}

#' map_to_parent
#' 
#' function to get binary dataset with gate assignation based on the mother population of the selected pop (instead of the root).
#' @param gh GatingHierarchy.
#' @param binary_df Dataframe generated by the map_to_root function.
#' @return Dataframe.
#' @keywords flowMagic
#' @export
#' @examples 
#' \donttest{map_to_parent()}

map_to_parent<-function(gh,binary_df){
  current_pop<-colnames(binary_df)[length(colnames(binary_df))]
  print(sprintf("finding mother pop of : %s",current_pop))
  mother_pop<-tail(strsplit(gs_pop_get_parent(gh,current_pop),"/")[[1]],1)
  print(sprintf("mother pop: %s",mother_pop))
  #------  get binary dataset for mother pop
  print("get binary dataset of the mother pop based on the root")
  df_tot_mother<-map_to_root(gh=gh,pop=mother_pop) # root dataset indicating with the integer 1,the mother cells of pop of the binary df
  df_sub_mother<-subset(df_tot_mother,df_tot_mother[,length(colnames(df_tot_mother))]==1) # we select only the mother cells from the root
  indices<-rownames(df_sub_mother) # this command finds the indices of the mother cells of the pop of binary_df from the root.
  # ----- get modified binary dataset
  print("get modified binary dataset with indices based on the mother pop")
  binary_df<-binary_df[indices,1:length(colnames(binary_df))] # from the root of binary_df we select only the cells of the mother pop of binary_df.
  # In this case the mother cells are indicated with 0.
  inds<-which(binary_df[,ncol(binary_df)]==1)
  print(sprintf("n_gate_events: %s",length(inds)))
  print(sprintf("n_mother_events: %s",nrow(binary_df)))
  print("Done")
  return(binary_df) # modified binary_df
  
}

#' get_pop_multiclass
#' 
#' function to get the children of the selected population, 
#' useful for multiclass classification. The output is only the pops with same dimensions 
#' @param gh GatingHierarchy.
#' @param pop Name of population to get events assignments.
#' @return Vector of characters.
#' @keywords flowMagic
#' @export
#' @examples 
#' \donttest{get_pop_multiclass()}


get_pop_multiclass<-function(gh,pop){
  child_pops_x_list<-gs_pop_get_children(gh,pop) # get children selected pop
  # if the selected pop has more than one child means that we need to classify more than one gate for this gating step
  # However even if the pop has more than one child the multiclass situation will occur only when there are
  # more than one child with the same dimensions
  ########################### case in which the selected pop has only one child ######################
  if(length(child_pops_x_list)==1){ # the input pop has only one child we don't need to do anything else.
    child_name<-tail(strsplit(child_pops_x_list[1],"/")[[1]],1)
    return(sprintf("mother_%s:%s",pop,child_name))
  }
  ########################### case in which the selected pop has more than one child ######################
  if(length(child_pops_x_list)>=1){
    ########### we get several info about the children populations and their dimensions
    child_pop_X_vector<-sapply(1:length(child_pops_x_list), function(i){
      pop<-tail(strsplit(child_pops_x_list[i],"/")[[1]],1)
      return(pop)
    })
    print(sprintf("All children pops of selected pop: %s",paste0(child_pop_X_vector,collapse = ",")))
    s<-sapply(1:length(child_pop_X_vector),function(i){
      polgate<-gh_pop_get_gate(gh,child_pop_X_vector[i])
      if(class(polgate)=="polygonGate"){
        dims<-colnames(polgate@boundaries)
        return(sprintf("%s_%s:%s",dims[1],dims[2],child_pop_X_vector[i]))
      }else if(class(polgate)=="rectangleGate"){
        rectGate<-polgate
        vec<-rectGate@min
        dims<-names(vec)
        # in case there is only one dimenion indicated in the gh,as second dimension we use by default SSC_A
        if(length(dims)==1){
          dims<-c(dims,"SSC-A")
        }
        return(sprintf("%s_%s:%s",dims[1],dims[2],child_pop_X_vector[i]))
      }else if(class(polgate)=="booleanFilter"){
        # we use the dimenions of the parent (so the current pop under analysis)
        polgate_parent<-gh_pop_get_gate(gh,pop)
        if(class(polgate_parent)=="polygonGate"){
          dims<-colnames(polgate_parent@boundaries)
          return(sprintf("%s_%s:%s",dims[1],dims[2],child_pop_X_vector[i]))
        }else if(class(polgate_parent)=="rectangleGate"){
          rectGate<-polgate_parent
          vec<-rectGate@min
          dims<-names(vec)
          # in case there is only one dimension indicated in the gh,as second dimension we use by default SSC_A
          if(length(dims)==1){
            dims<-c(dims,"SSC-A")
          }
        }else if(class(polgate_parent)=="booleanFilter"){
          stop( "class of the parent is also boolean")
        }
      }else{
        print(polgate)
        stop("class polgate not recognized")
      }
    })
    print("---- all children pops with their dimensions")
    print(s)
    # the s variable reports the dimensions of the children populations and the name of the children populations
    c<-strsplit(s,":")
    s2<-sapply(1:length(c),function(i){
      return(c[[i]][1])
    })
    # s2= only the dimensions without the name
    dims_set<-unique(s2) # unique dimensions set (vectors that contains all the dimensions of the children pops)
    #################### Based on the dimensions set we divide the children in two categories
    # we look for pops that share the same dimensions (same_dims category) 
    # and pops that have distinct dimensions from the rest (diff dims category)
    s3<-lapply(1:length(dims_set),function(i){
      inds_dims_set_i<-which(s2==dims_set[i]) # ith element of the dimensions set
      pops_current_dims_set<-child_pop_X_vector[inds_dims_set_i] # pops that have these dimensions
      if(length(pops_current_dims_set)>=2){ # 2 or more pops share the same dimensions (same entry in the list)
        pops<-paste0(pops_current_dims_set,collapse = ",")
        pops<-paste0(pops,"#same_dims")
        return(pops)
      }else if(length(pops_current_dims_set)==1){ # only 1 pop has these dimensions
        pop<-paste0(pops_current_dims_set,"#diff_dims")
        return(pop)
      }
    })
    # s3 reports the children populations of the input pop that have the same dimensions,as a list.
    print("---- all children pops with same dimensions")
    print(s3)
    # we generate the final categorization of the children pops of the selected pop
    final_s<-sapply(1:length(s3), function(i){
      entry_i<-s3[[i]]
      return(entry_i)
    })
    final_s<-paste0(final_s,collapse = ";") # contains the final categorization of the children pops
    print("------ final categorization in a unique vector")
    return(paste0(sprintf("mother_%s:%s",pop,final_s)))
  }else{
    ########################### case in which the selected pop has no children ######################
    # we just report this info
    return(paste0(sprintf("mother_%s:no_children",pop)))
  }
}

#' get_hierarchy_all_pops
#' 
#' function to get the hierarchy of all pops from the sample manually gated (the input gating hierarchy).
#' Based on its output the function magicTrain will perform the training step using 
#' the sample manually gated (local training set) and the project discovery data (global training set).
#' @param gh GatingHierarchy.
#' @param export_visnet If true, it export visnetwork object.
#' @param path.output Path to save visnetwork object. Default to None.
#' @return List of Dataframes.
#' @keywords flowMagic
#' @export
#' @examples 
#' \donttest{get_hierarchy_all_pops()}


get_hierarchy_all_pops<-function(gh,export_visnet=F,path.output="None"){
  start<-Sys.time()
  ##### check loading of necessary libraries #############
  if (("package:stringr" %in% search()) == F) {
    warning("stringr library not loaded, attempt to loading...")
    library(stringr)
  }
  if (("package:data.tree" %in% search()) == F) {
    warning("data.tree library not loaded, attempt to loading...")
    library(data.tree)
  }
  if (("package:visNetwork" %in% search()) == F) {
    warning("visNetwork library not loaded, attempt to loading...")
    library(visNetwork)
  }
    if (("package:randomcoloR" %in% search()) == F) {
    warning("randomcoloR library not loaded, attempt to loading...")
    library(randomcoloR)
  }
  ################################ get info about the hierarchy #####################
  all_pops<-name_pop_gating(gh)
  print("------- list all pops in gh--------")
  print(all_pops)
  print("-----------------------------------")
  print(paste0("total number of pops:",length(all_pops)))
  ################## calculating the hierarchy of the input gh #######################
  # make df that contains the hierarchy info
  print("########### make df with hierarchy info ########")
  df<-data.frame(Children=character(),Mother=character(),stringsAsFactors=FALSE)
  for(i in 1:length(all_pops)){
    pop<-all_pops[i]
    # get children info current pop
    pops_multiclass<-get_pop_multiclass(gh,pop) 
    s<-strsplit(pops_multiclass,":")
    # analyze info children current pop
    for (t in 1:length(s)){
      mother<-s[[t]][1]
      mother<-strsplit(mother,"mother_")[[1]][2]
      all_children_info<-s[[t]][2]
      all_children_info<-strsplit(all_children_info,";")[[1]]
      for (children_info in all_children_info){
        new_row<-cbind(mother,children_info)
        df<-rbind(df,new_row)
      }
    }
  }
  # reset indices dataframe
  rownames(df)<-NULL
  df <- data.frame(lapply(df, as.character), stringsAsFactors=FALSE)
  colnames(df)<-c("Mother","Children")
  # remove no children entry
  inds<-which(df$Children=="no_children")
  df<-df[-inds,]
  # rearrange df
  df_v2<-data.frame(Children=character(),Mother=character(),stringsAsFactors=FALSE)
  for (i in 1:nrow(df)){
    row_i<-df$Children[i]
    ind<-grep("#",row_i)
    if(length(ind)>0){
      str<-strsplit(row_i,"#")[[1]][1]
      str<-strsplit(str,",")[[1]]
      for (pop in str){
        new_row<-data.frame(df$Mother[i],pop)
        names(new_row)<-c("Mother","Children")
        df_v2<-rbind(df_v2,new_row)
      }
    }else{
      new_row<-data.frame(df$Mother[i],df$Children[i])
      names(new_row)<-c("Mother","Children")
      df_v2<-rbind(df_v2,new_row)
    }
  }
  #------------------- make string hierarchical path column ----------------------------------------
  print("##### make path column #######")
  #define pathstring for each row of df
  pathString_col<-data.frame(pathString=character(),stringsAsFactors=FALSE)
  for (i in 1:nrow(df_v2)){
    row_i_child<-as.character(df_v2$Children[i])
    row_i_mother<-df_v2$Mother[i]
    path_i<-paste(gs_pop_get_parent(gh,row_i_child,path="full"),row_i_child,sep="/")
    if(i>1){
      path_i<-paste0("root",path_i)
    }
    path_i<-data.frame(path_i,stringsAsFactors=FALSE)
    names(path_i)<-c("pathString")
    pathString_col<-rbind(pathString_col,path_i)
  }
  if(nrow(df_v2)!=nrow(pathString_col)){
    stop("error extraction gating hierarchy:nrow(df_v2)!=nrow(pathString_col)")
  }
  # ----------------------- make dimension info column ------------------
  print("##### make Dimensions column #######")
  inds<-grep("same_dims",df$Children)
  dim_col<-character(nrow(df_v2))
  counter<-0
  for (ind in inds){
    counter<-counter+1
    child_ind<-df$Children[ind]
    str<-strsplit(child_ind,"#")[[1]][1]
    str<-strsplit(str,",")[[1]]
    logic<-df_v2$Children %in% str
    inds<-which(logic==T)
    dim_col[inds]<-sprintf("same_dims_%s",counter)
  }
  df_v2$Dimensions<-dim_col
  #--------------------- generate hierarchical tree df ----------------------------
  print("##### build final hierarchical pop tree ######")
  df_tree<-cbind(df_v2,pathString_col)
  df_tree$Mother<-as.character(df_tree$Mother)
  df_tree$Children<-as.character(df_tree$Children)
  df_tree$Dimensions<-as.character(df_tree$Dimensions)
  df_tree_v2<-df_tree
  df_tree_v2$Mother[1]<-"root2"
  df_tree_v2$pathString<-str_replace(df_tree_v2$pathString,"root","root2")
  hierarchical_tree <- as.Node(df_tree_v2,df_tree_v2$Mother,df_tree_v2$Children)
  df_tree_levels<-ToDataFrameTree(hierarchical_tree, "level")
  #------------------ make visNetwork hierarchy plot ---------------
  print("######## make hierarchical visNetwork ######")
  # nodes colored in the same way share the same dims (stay in same plot)
  # convert df_tree in the two visnet dfs
  # ---- make edges df
  df_tree_edges<-df_tree[,c(1,2)]
  colnames(df_tree_edges)<-c("from","to")
  # ---- make nodes df
  colors_v<-rep("white",nrow(df_tree)) # white is the default color for the nodes 
  # the populations single children or with different dims (separate plots)
  unique_dim<-unique(df_tree$Dimensions) # "" entry refers to situation with pops with different dims,
  # so separate plots (we color them with default color)
  ind<-which(unique_dim=="")
  if(all(unique_dim=="")==F){
    if(length(ind)!=0){
      unique_dim<-unique_dim[-ind]
    }
    n <- length(unique_dim)
    palette <- distinctColorPalette(n)
    set.seed(123)
    inds_col<-sample(1:length(palette),length(unique_dim)) # generate random indices
    i<-0
    for(dim in unique_dim){
      i<-i+1
      inds<-which(df_tree$Dimensions==dim)
      colors_v[inds]<-palette[inds_col[i]]
    }
  }
  df_tree_nodes<-cbind(df_tree$Children,df_tree$Children,colors_v)
  df_tree_nodes<-as.data.frame(df_tree_nodes)
  colnames(df_tree_nodes)<-c("id","label","color.background")
  # add root row
  root_row<-data.frame("root","root","white")
  colnames(root_row)<-c("id","label","color.background")
  df_tree_nodes<-rbind(root_row,df_tree_nodes)
  # add border color
  colorboder<-rep("black",nrow(df_tree_nodes))
  df_tree_nodes<-cbind(df_tree_nodes,colorboder)
  colnames(df_tree_nodes)<-c("id","label","color.background","color.border")
  #------ make visnetwork plot
  visnet<-visNetwork(edges=df_tree_edges,nodes=df_tree_nodes)
  visnet<- visnet %>% visEdges(arrows = "to") %>% visHierarchicalLayout() %>% visNodes(borderWidth = 2)
  visnet<-visnet %>% visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE)
  #-------------------- export visnetwork -------------------------
  if(export_visnet==T){
    dir.create(paste0(path.output,"/Gatinghierarchy_Plots/"),recursive = T)
  }
  end<-Sys.time()
  time_taken<-end-start
  print("Time of execution:")
  print(time_taken)
  print("Done")
  return(list(df_tree=df_tree,df_tree_levels=df_tree_levels,visnet=visnet,hierarchical_tree=hierarchical_tree))
}


#' get_local_train_sets
#' 
#' Based on the hierarchy calculated using the get_hierarchy function, we generate the local training sets gated by the biologists.
#' @param gh GatingHierarchy.
#' @param hierarchical_tree Dataframe of hierarchy information generated by the get_hierarchy function.
#' @param info_hierarchy List of hierarchy information generated by the get_hierarchy function.
#' @return List of Dataframes.
#' @keywords flowMagic
#' @export
#' @examples 
#' \donttest{get_local_train_sets()}

# -------- function to generate the list of the training sets hierarchically distributed
#
get_local_train_sets<-function(gh,hierarchical_tree,info_hierarchy){
  start<-Sys.time()
  
  #sink('/home/rstudio/results/logs/analysis-output.txt')
  ######### get info levels of input hierarchy ########
  df_tree_levels<-info_hierarchy$df_tree_levels
  df_tree<-info_hierarchy$df_tree
  pops<-names(hierarchical_tree$Get('level'))
  levels_for_pop<-as.vector(hierarchical_tree$Get('level'))
  unique_levels<-unique(levels_for_pop)
  max_level<-max(unique_levels) # max_level is the lowest level of the hierarchy
  print(max_level)
  list_train_sets_all_levels<-list() # list that reports the train sets for each parent node of each level of the hierarchy
  ######## generate training sets for all levels
  for(level in 2:max_level){ # we start from the second level, because the root has no mother.
    print(sprintf("##################################### train sets of level %s #######################################",level))
    # we select the pops of the current level
    inds<-which(levels_for_pop==level)
    selected_pops_current_level<-pops[inds]
    print(sprintf("all pops current level:%s",paste0(selected_pops_current_level,collapse = " ")))
    # we select the mothers of all pops of the current level
    string<-paste0(selected_pops_current_level,collapse = "|")
    logic<-df_tree$Children %in% selected_pops_current_level
    inds<-which(logic==T)
    if(length(inds)==0){
      stop(sprintf("fatal error occured. No pops found for %s",selected_pops_current_level))
    }
    mothers_all_pops_current_level<-as.character(df_tree$Mother[inds])
    print(sprintf("mothers pops current level:%s",paste0(mothers_all_pops_current_level,collapse = " ")))
    # dimensions set of the pops of the current level
    dims_info_pops_current_level<-as.character(df_tree$Dimensions[inds])
    print(sprintf("info_dims_current_level:%s",paste0(dims_info_pops_current_level,collapse = " ")))
    # we test if some of the pops of the current level have same dims
    inds_same_dims_pops<-grep("same_dims",dims_info_pops_current_level)
    if(length(inds_same_dims_pops)>0){ # some populations with same dims in the current level
      print("some populations with same dims in the current level")
      selected_pops_current_level_same_dims<-selected_pops_current_level[inds_same_dims_pops]
      selected_pops_current_level_different_dims<-selected_pops_current_level[-inds_same_dims_pops]
    }else{ # all populations have different dims in the current level
      print("all populations have different dims in the current level")
      selected_pops_current_level_different_dims<-selected_pops_current_level
    }
    list_df_pops_current_level<-list()
    # -------------------- train set of pops with different dims --------
    # we generate the train sets of pops with different dimensions for the current level (binary dataset where cells associated to the pop are 1)
    if(length(selected_pops_current_level_different_dims)!=0){
      print("---- Generate train set of pops with different dims -----")
      vec_labels_pops_diff_dims<-c()
      for(pop in selected_pops_current_level_different_dims){
        binary_df_root_pop<-map_to_root(gh=gh,pop=pop)
        binary_df_mother_pop<-map_to_parent(gh=gh,binary_df=binary_df_root_pop) # train set current pop
        list_df_pops_current_level[[sprintf("%s",pop)]]<-binary_df_mother_pop
        vec_labels_pops_diff_dims<-append(vec_labels_pops_diff_dims,sprintf("%s:%s",pop,1))
      }
      list_df_pops_current_level[[sprintf("labels_%s",paste0(selected_pops_current_level_different_dims,collapse = "_"))]]<-vec_labels_pops_diff_dims
    }
    # -------------------- train set of pops with same dims --------
    if(length(inds_same_dims_pops)>0){
      print("---- Generate train set of pops with same dims----")
      set_same_dims<-unique(dims_info_pops_current_level)
      inds_remove_empty<-which(set_same_dims=="")
      if(length(inds_remove_empty)>0){
        set_same_dims<-set_same_dims[-inds_remove_empty] 
      }
      print(sprintf("set_same_dims:%s",paste0(set_same_dims,collapse = " ")))
      # print(dims_info_pops_current_level)
      # print(selected_pops_current_level)
      for(element in set_same_dims){
        print(sprintf("Current set in analysis: %s",element))
        inds<-which((df_tree$Dimensions %in% element)==T)
        pops_same_plot<-df_tree$Children[inds]
        print(sprintf("current pops same plot:%s",paste0(pops_same_plot,collapse = " ")))
        pop_temp<-pops_same_plot[1]
        binary_df_root_pop<-map_to_root(gh=gh,pop=pop_temp)
        multiclass_df_mother_pop_temp<-map_to_parent(gh=gh,binary_df=binary_df_root_pop) # needs to be modified with the correct labels
        dims_element<-colnames(multiclass_df_mother_pop_temp)[c(1,2)]
        vec_rownames_all_pop<-c() # vector that contains the rownames of all pops of the current element
        for (pop in pops_same_plot){
          binary_df_root_pop<-map_to_root(gh=gh,pop=pop)
          binary_df_mother_pop<-map_to_parent(gh=gh,binary_df=binary_df_root_pop) # train set current pop. Binary,but we need it multiclass
          inds<-which(binary_df_mother_pop[,ncol(binary_df_mother_pop)]==1)
          binary_df_mother_pop<-binary_df_mother_pop[inds,]
          rownames_pop<-row.names(binary_df_mother_pop)
          vec_rownames_all_pop<-append(vec_rownames_all_pop,sprintf("%s:%s",pop,paste0(rownames_pop,collapse = "_")))
        }
        # rename multiclass mother pop with correct multiclass labels 
        l<-0 # counter of the labels
        vec_labels_same_dims_current_element<-c()
        for(p in vec_rownames_all_pop){
          l<-l+1
          splitstring<-strsplit(p,":")[[1]]
          pop<-splitstring[1]
          rownames<-splitstring[2]
          rownames<-strsplit(rownames,"_")[[1]]
          multiclass_df_mother_pop_temp[rownames,ncol(multiclass_df_mother_pop_temp)]<-l
          vec_labels_same_dims_current_element<-append(vec_labels_same_dims_current_element,sprintf("%s:%s",pop,l))
        }
        colnames(multiclass_df_mother_pop_temp)[ncol(multiclass_df_mother_pop_temp)]<-paste0(pops_same_plot,collapse="_")
        multiclass_df_mother_pop<-multiclass_df_mother_pop_temp # now the df is truly multiclass
        list_df_pops_current_level[[sprintf("%s",paste0(pops_same_plot,collapse="_"))]]<-multiclass_df_mother_pop
        list_df_pops_current_level[[sprintf("labels_%s",paste0(pops_same_plot,collapse="_"))]]<-vec_labels_same_dims_current_element
      }
    }
    list_train_sets_all_levels[[sprintf("level:%s",level)]]<-list_df_pops_current_level
  }
  #sink('/home/rstudio/results/logs/analysis-output.txt', append=TRUE)
  #closeAllConnections()
  end<-Sys.time()
  time_taken<-end-start
  print("training time:")
  print(time_taken)
  print("Done")
  return(list_train_sets_all_levels)
  
}



