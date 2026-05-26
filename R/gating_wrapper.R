
#' magic_manual_gating
#'
#' Perform interactive manual gating on selected samples and add or update the
#' resulting polygon gate in a GatingSet.
#'
#' This function wraps `magicGating()` for one manual polygon gate and integrates
#' the resulting `polygonGate` into a `GatingSet`. If the target child node does
#' not already exist, it is created under `parent_node`. If the node already
#' exists, the gate is updated only for the samples present in the gating output.
#'
#' @param gs_input A `GatingSet` object.
#' @param parent_node Character. Name of the parent node where the new gate should
#'   be added or updated.
#' @param new_child_node Character. Name of the child node to create or update.
#' @param samples_id Character vector. Sample names to manually gate. These names
#'   must match `sampleNames(gs_input)`.
#' @param channel_x Character. Name of the channel to plot on the x-axis.
#' @param channel_y Character. Name of the channel to plot on the y-axis.
#' @param size_points Numeric. Point size used in the interactive gating plot.
#'   Default is `0.5`.
#' @param return_list_out Logical. If `TRUE`, return only the output from
#'   `magicGating()` without modifying the `GatingSet`. Default is `FALSE`.
#' @param input_list_out Optional list. Existing output from `magicGating()`
#'   containing a `list_poly_gates` slot. If provided, interactive gating is
#'   skipped and these polygon gates are added or updated in the `GatingSet`.
#' @param return_gs_only Logical. If `TRUE`, return only the modified `GatingSet`.
#'   If `FALSE`, return a list containing the modified `GatingSet` and the
#'   gating output. Default is `TRUE`.
#' @param ... Additional arguments passed to `magicGating()`.
#'
#' @return If `return_gs_only = TRUE`, a modified `GatingSet`. Otherwise, a list
#'   with two elements: `gs`, the modified `GatingSet`, and `list_out`, the
#'   gating output used to add or update the gate.
#'
#' @keywords flowMagic 
#' @export
#'
#' @examples
#' \donttest{
#' gs <- magic_manual_gating(
#'   gs_input = gs,
#'   parent_node = "root",
#'   new_child_node = "Live_cells",
#'   samples_id = c("sample_01.fcs"),
#'   channel_x = "FSC-A",
#'   channel_y = "LIVE DEAD Blue-A"
#' )
#' }

magic_manual_gating<-function(gs_input,parent_node,new_child_node,samples_id,
                              channel_x,channel_y,size_points=0.5,return_list_out=F,
                              input_list_out=NULL,return_gs_only=T,...){
  
  message("$$$ checking input $$$")

  #------ check correct input format -------
  if(is.character(parent_node)==F){
    stop("parent_node must be a character")
  }
  if(is.character(new_child_node)==F){
    stop("new_child_node must be a character")
  }

  # ----- Perform interactive gating ------
  if(is.null(input_list_out)==T){
    message("$$$ Perform interactive gating $$$ ")
  
    list_out<-magicGating(fs=gs_input, gs_node = parent_node,sample_id = samples_id,
                          channel_x = channel_x,channel_y = channel_y,size_points=size_points,
                          label_pol = "1",...)
    
    if(return_list_out==T){
      message("only list_out returned")
      return(list_out)
    }
    
  }else{
    message("$$$ no interactive gating, external list of polygon gates provided $$$ ")
    check_list <- "list_poly_gates" %in% names(input_list_out)
    if(check_list==F){
      stop("external list provided does not have a slot named list_poly_gates")
    }
    list_out<-input_list_out
  }
  
  # ----- Integrate polygonGate into the GatingSet -----
  message("$$$ Integrate polygonGate into the GatingSet $$$ ")
  
  current_nodes_names<-basename(flowWorkspace::gs_get_pop_paths(gs_input, path = "full"))
  
  # check if child node already in GatingSet
  if((new_child_node %in% current_nodes_names)==F){
    warning("child node name selected not in hierarchy, making new node with first gate 
            and apply other gates based on samples names")
    flowWorkspace::gs_pop_add(gs_input, list_out$list_poly_gates[[1]], parent = parent_node, name = new_child_node)
  }
  
  # ---- update GatingSet ------
  names_samples_gated<-names(list_out$list_poly_gates)
  names_samples_in_gs<-sampleNames(gs_input)
  
  # check for error in names_samples_gated
  missing_samples <- setdiff(names_samples_gated, names_samples_in_gs)
  
  if (length(missing_samples) > 0) {
    stop(
      "These gated samples are not present in sampleNames(gs_input): ",
      paste(missing_samples, collapse = ", ")
    )
  }
  # It will still allow partial updates, for example 2 out of 10 samples. 
  # It only stops when one of the provided sample names does not exist in the GatingSet.
  
  inds_samples_to_update<-which(names_samples_in_gs %in% names_samples_gated)
  flowWorkspace::gs_pop_set_gate(gs_input[inds_samples_to_update],value = list_out$list_poly_gates,y = new_child_node)
  
  # ------- Recompute the gating -----
  message("recompute GatingSet")
  recompute(gs_input)
  
  message("Done")
  if(return_gs_only==T){
    message("return only GatingSet")
    return(gs_input)
  }else{
    message("return both GatingSet and list gated data")
    return(list(gs=gs_input,list_out=list_out))
  }
}

#' magic_template_train
#'
#' Build a template training set from manually gated samples and train a
#' flowMagic classification model.
#'
#' This function performs one or more rounds of interactive manual gating on
#' selected samples. Each gate is assigned a different label. If more than one
#' gate is drawn, the labeled outputs are merged into a single template data
#' object. The resulting templates are converted into training data using
#' `get_train_data()` and used to train a model with `magicTrain()`.
#'
#' Existing template data can be supplied through `input_list_out_final`, allowing
#' users to add new manually gated templates without redrawing previous ones.
#'
#' @param gs_input A `GatingSet` object.
#' @param parent_node Character. Name of the parent node from which expression
#'   data are extracted for template gating.
#' @param samples_id Character vector. Sample names to manually gate as template
#'   samples.
#' @param channel_x Character. Name of the channel to plot on the x-axis.
#' @param channel_y Character. Name of the channel to plot on the y-axis.
#' @param size_points Numeric. Point size used in the interactive gating plot.
#'   Default is `0.5`.
#' @param n_gates Integer. Number of template gates to draw per selected sample.
#'   Each gate is assigned a label from `1` to `n_gates`. Default is `1`.
#' @param input_list_out_final Optional list. Existing template data object,
#'   usually from a previous call to `magic_template_train()$list_out_final`.
#'   If supplied, newly gated templates are added to this existing template set.
#'   Samples with matching names are replaced.
#' @param train_model Character. Model type passed to `magicTrain()`. Default is
#'   `"rf"`.
#' @param ... Additional arguments passed to `magicGating()`.
#'
#' @return A list with five elements: `ref_model_info`, the trained model
#'   information returned by `magicTrain()`; `ref_train`, the training data
#'   returned by `get_train_data()`; `list_out_final`, the complete merged
#'   template set; `list_out_new`, the newly generated template data; and
#'   `all_list_out_obj`, the raw outputs from each call to `magicGating()`.
#'
#' @keywords flowMagic template training machine-learning gating
#' @export
#'
#' @examples
#' \donttest{
#' train_out <- magic_template_train(
#'   gs_input = gs,
#'   parent_node = "root",
#'   samples_id = c("sample_01.fcs", "sample_02.fcs"),
#'   channel_x = "FSC-A",
#'   channel_y = "LIVE DEAD Blue-A",
#'   n_gates = 2,
#'   train_model = "rf"
#' )
#'
#' # Add additional template samples later
#' train_out2 <- magic_template_train(
#'   gs_input = gs,
#'   parent_node = "root",
#'   samples_id = c("sample_03.fcs"),
#'   channel_x = "FSC-A",
#'   channel_y = "LIVE DEAD Blue-A",
#'   n_gates = 2,
#'   input_list_out_final = train_out$list_out_final
#' )
#' }

magic_template_train <- function(gs_input, parent_node, samples_id,
                                 channel_x, channel_y,
                                 size_points = 0.5,
                                 n_gates = 1,
                                 input_list_out_final = NULL,
                                 train_model = "rf",
                                 ...) {
  
  message("$$$ Starting template training workflow $$$")
  message(sprintf("Parent node: %s", parent_node))
  message(sprintf("Channels: %s vs %s", channel_x, channel_y))
  message(sprintf("Number of template gates to draw per sample: %s", n_gates))
  message(sprintf("Number of samples selected for template gating: %s", length(samples_id)))
  
  if (!is.null(input_list_out_final)) {
    message("Existing template object detected.")
    message(sprintf("Existing templates: %s", length(input_list_out_final)))
    message("New templates will be added to the existing template set.")
  } else {
    message("No existing template object provided. A new template set will be created.")
  }
  
  message("$$$ Perform interactive template gating $$$")
  
  all_list_out_obj <- list()
  
  for (n in 1:n_gates) {
    
    message(sprintf("---- Gating template/gate %s of %s ----", n, n_gates))
    message(sprintf("Current label assigned to this gate: %s", as.character(n)))
    
    list_out_n <- magicGating(
      fs = gs_input,
      gs_node = parent_node,
      sample_id = samples_id,
      channel_x = channel_x,
      channel_y = channel_y,
      size_points = size_points,
      label_pol = as.character(n),
      ...
    )
    
    all_list_out_obj[[as.character(n)]] <- list_out_n
    
    message(sprintf("Finished gating label %s.", n))
  }
  
  message("$$$ Combining template labels $$$")
  
  if (n_gates == 1) {
    
    message("Only one gate was drawn. No label merging needed.")
    list_out_new <- all_list_out_obj[[1]]$list_gated_data
    
  } else {
    
    message(sprintf("Merging labels from %s gates.", n_gates))
    
    all_gated_data <- lapply(all_list_out_obj, function(x) x$list_gated_data)
    
    list_out_new <- Reduce(
      f = function(x, y) {
        merge_magicGating_labels(
          list_out_1 = x,
          list_out_2 = y,
          gated_data_only  = TRUE
        )
      },
      x = all_gated_data
    )
    
    message("Finished merging template labels.")
  }
  
  message("$$$ Updating template set $$$")
  
  if (is.null(input_list_out_final)) {
    
    message("Creating new template set from current gated samples.")
    list_out_final <- list_out_new
    
  } else {
    
    duplicated_samples <- intersect(
      names(input_list_out_final),
      names(list_out_new)
    )
    
    if (length(duplicated_samples) > 0) {
      warning(
        "These samples already exist in the template set and will be replaced: ",
        paste(duplicated_samples, collapse = ", ")
      )
    }
    
    new_samples <- setdiff(
      names(list_out_new),
      names(input_list_out_final)
    )
    
    if (length(new_samples) > 0) {
      message(
        "Adding new template samples: ",
        paste(new_samples, collapse = ", ")
      )
    }
    
    list_out_final <- input_list_out_final
    list_out_final[names(list_out_new)] <- list_out_new
    
    message(sprintf("Updated template set now contains %s samples.", length(list_out_final)))
  }
  
  message("$$$ Convert templates to training data $$$")
  
  ref_train <- get_train_data(paths_file = list_out_final)
  
  message(sprintf("Training data created with %s events.", nrow(ref_train)))
  message(sprintf("Training classes detected: %s", paste(sort(unique(ref_train$classes)), collapse = ", ")))
  
  message("$$$ Train model $$$")
  message(sprintf("Training model type: %s", train_model))
  
  ref_model_info <- magicTrain(
    df_train = ref_train,
    train_model = train_model
  )
  
  message("$$$ Template training workflow completed $$$")
  
  return(
    list(
      ref_model_info = ref_model_info,
      ref_train = ref_train,
      list_out_final = list_out_final,
      list_out_new = list_out_new,
      all_list_out_obj = all_list_out_obj
    )
  )
}


#' magic_template_predict
#'
#' Predict gated populations from a trained flowMagic template model and
#' optionally add the predicted gates to a GatingSet.
#'
#' This function exports bivariate expression data from a selected GatingSet node,
#' predicts event labels using a previously trained flowMagic model, converts
#' selected predicted labels into polygon gates, and optionally adds or updates
#' these gates in the input `GatingSet`.
#'
#' The `label_to_node` argument controls how predicted labels are mapped to
#' GatingSet node names. For example, `c("1" = "Live_cells", "2" = "Dead_cells")`
#' converts predicted label `"1"` into a node named `"Live_cells"` and predicted
#' label `"2"` into a node named `"Dead_cells"`.
#'
#' @param gs_input A `GatingSet` object.
#' @param parent_node Character. Name of the parent node from which raw data are
#'   exported and under which predicted gates are added.
#' @param channel_x Character. Name of the channel to use on the x-axis.
#' @param channel_y Character. Name of the channel to use on the y-axis.
#' @param ref_train Data frame. Training data generated by `get_train_data()`,
#'   typically `train_out$ref_train` from `magic_template_train()`.
#' @param ref_model_info List. Trained model information returned by
#'   `magicTrain()`, typically `train_out$ref_model_info` from
#'   `magic_template_train()`.
#' @param label_to_node Named character vector. Names are predicted labels and
#'   values are the corresponding GatingSet node names. For example,
#'   `c("1" = "Live_cells", "2" = "Dead_cells")`.
#' @param n_cores Integer. Number of cores used for prediction and polygon
#'   extraction. Default is `1`.
#' @param concavity_val Numeric. Concavity parameter passed to
#'   `flowmagic_pred_to_poly_gates()`. Higher values generally produce smoother
#'   polygon gates. Default is `50`.
#' @param add_to_gs Logical. If `TRUE`, convert predicted labels to polygon gates
#'   and add or update them in the `GatingSet`. If `FALSE`, only prediction
#'   results are returned. Default is `TRUE`.
#' @param return_gs_only Logical. If `TRUE`, return only the modified `GatingSet`.
#'   Default is `FALSE`.
#' @param ... Additional arguments passed to `export_raw_gs_plots()`.
#'
#' @return If `return_gs_only = TRUE`, a modified `GatingSet`. Otherwise, a list
#'   with four elements: `gs`, the modified or original `GatingSet`;
#'   `list_test_data`, the exported test data; `list_dfs_pred`, the predicted
#'   labeled data frames; and `list_poly_gates`, the polygon gates generated for
#'   each requested node.
#'
#' @keywords flowMagic template prediction automated gating GatingSet
#' @export
#'
#' @examples
#' \donttest{
#' pred_out <- magic_template_predict(
#'   gs_input = gs,
#'   parent_node = "root",
#'   channel_x = "FSC-A",
#'   channel_y = "LIVE DEAD Blue-A",
#'   ref_train = train_out$ref_train,
#'   ref_model_info = train_out$ref_model_info,
#'   label_to_node = c("1" = "Live_cells", "2" = "Dead_cells"),
#'   n_cores = 1,
#'   add_to_gs = TRUE
#' )
#'
#' gs <- pred_out$gs
#' }

magic_template_predict <- function(gs_input,
                                   parent_node,
                                   channel_x,
                                   channel_y,
                                   ref_train,
                                   ref_model_info,
                                   label_to_node,
                                   n_cores = 1,
                                   concavity_val = 50,
                                   add_to_gs = TRUE,
                                   return_gs_only = FALSE,
                                   ...) {
  
  message("$$$ Starting template prediction workflow $$$")
  
  # check inputs
  if (!is.character(parent_node)) {
    stop("parent_node must be a character")
  }
  
  if (!is.character(channel_x) || !is.character(channel_y)) {
    stop("channel_x and channel_y must be characters")
  }
  
  if (is.null(names(label_to_node))) {
    stop("label_to_node must be a named character vector, e.g. c('1' = 'Live_cells', '2' = 'Dead_cells')")
  }
  
  if (is.null(ref_train)) {
    stop("ref_train is NULL.")
  }
  
  if (is.null(ref_model_info)) {
    stop("ref_model_info is NULL.")
  }
  
  message(sprintf("Parent node: %s", parent_node))
  message(sprintf("Channels: %s vs %s", channel_x, channel_y))
  message(sprintf("Labels to export: %s", paste(names(label_to_node), collapse = ", ")))
  message(sprintf("Target nodes: %s", paste(label_to_node, collapse = ", ")))
  
  # get test data
  message("$$$ Export raw data from GatingSet $$$")
  
  list_test_data <- flowMagic::export_raw_gs_plots(
    gs = gs_input,
    node_name = parent_node,
    channel_x = channel_x,
    channel_y = channel_y,
    return_data = TRUE,
    ...
  )
  
  message(sprintf("Exported data for %s samples.", length(list_test_data)))
  
  # perform automated gating
  message("$$$ Predict labels using template model $$$")
  
  list_dfs_pred <- magicPred_all(
    list_test_data = list_test_data,
    magic_model = NULL,
    ref_data_train = ref_train,
    ref_model_info = ref_model_info,
    n_cores = n_cores
  )

  message("Prediction completed.")
  
  # optionally add/update GatingSet
  list_poly_gates_all <- list()
  
  if (add_to_gs) {
    
    message("$$$ Convert predicted labels to polygon gates and update GatingSet $$$")
    
    for (pred_label in names(label_to_node)) {
      
      node_name <- label_to_node[[pred_label]]
      
      message(sprintf("---- Processing predicted label %s -> node %s ----", pred_label, node_name))
      
      
      list_poly_gates <- flowMagic::flowmagic_pred_to_poly_gates(
        list_df = list_dfs_pred,
        gate_label = node_name,
        pred_label = pred_label,
        n_cores = n_cores,
        concavity_val = concavity_val
      )
      
      list_poly_gates_all[[node_name]] <- list_poly_gates
      
      current_nodes_names <- basename(
        flowWorkspace::gs_get_pop_paths(gs_input, path = "full")
      )
      
      if (!(node_name %in% current_nodes_names)) {
        
        message(sprintf("Node %s does not exist. Adding it under %s.", node_name, parent_node))
        
        flowWorkspace::gs_pop_add(
          gs_input,
          list_poly_gates[[1]],
          parent = parent_node,
          name = node_name
        )
        
      } else {
        
        message(sprintf("Node %s already exists. Updating existing gates.", node_name))
      }
      
      names_samples_gated <- names(list_poly_gates)
      names_samples_in_gs <- sampleNames(gs_input)
      
      missing_samples <- setdiff(names_samples_gated, names_samples_in_gs)
      
      if (length(missing_samples) > 0) {
        stop(
          "These predicted samples are not present in sampleNames(gs_input): ",
          paste(missing_samples, collapse = ", ")
        )
      }
      
      inds_samples_to_update <- which(names_samples_in_gs %in% names_samples_gated)
      
      flowWorkspace::gs_pop_set_gate(
        gs_input[inds_samples_to_update],
        value = list_poly_gates,
        y = node_name
      )
    }
    
    message("recompute GatingSet")
    recompute(gs_input)
  }
  
  message("$$$ Template prediction workflow completed $$$")
  
  if (return_gs_only) {
    return(gs_input)
  }
  
  return(
    list(
      gs = gs_input,
      list_test_data = list_test_data,
      list_dfs_pred = list_dfs_pred,
      list_poly_gates = list_poly_gates_all
    )
  )
}


