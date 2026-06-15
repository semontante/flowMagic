
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
#' Existing template objects can be supplied through `input_template_obj`. This
#' allows users either to reuse the same manually gated templates for
#' reproducible training, or to use an existing template object as the starting
#' point before drawing new interactive templates.
#'
#' @param gs_input A `GatingSet` object.
#'
#' @param parent_node Character. Name of the parent node from which expression
#'   data are extracted for template gating.
#'
#' @param samples_id Character vector. Sample names to manually gate as template
#'   samples. All values must be present in `sampleNames(gs_input)`. Duplicate
#'   sample names are not allowed.
#'
#' @param channel_x Character. Name of the channel to plot on the x-axis.
#'
#' @param channel_y Character. Name of the channel to plot on the y-axis.
#'
#' @param n_gates Integer. Number of template gates to draw per selected sample.
#'   Each gate is assigned a label from `1` to `n_gates`. Default is `1`.
#'
#' @param input_template_obj Optional list. Existing template-training object,
#'   usually returned by a previous call to `magic_template_train()`. The object
#'   must contain `all_list_out_obj` if `use_existing_templates = TRUE`.
#'
#' @param use_existing_templates Logical. If `FALSE`, interactive gating is
#'   performed and new templates are drawn. If `TRUE`, `input_template_obj` must
#'   be supplied and the stored templates in
#'   `input_template_obj$all_list_out_obj` are reused without redrawing gates.
#'   This is useful for reproducing training from previously saved templates.
#'   Default is `FALSE`.
#'
#' @param train_model Character. Model type passed to `magicTrain()`. Default is
#'   `"rf"`.
#'
#' @param prop_down Numeric or `NULL`. Proportion of events to keep when building
#'   the training dataframe with `get_train_data()`. If `NULL`, no proportion is
#'   specified explicitly. Default is `NULL`.
#'
#' @param n_points_per_plot Numeric, integer, or `NULL`. Approximate number of
#'   points to keep per gated dataset when building the training dataframe with
#'   `get_train_data()`. This is useful to avoid training on all events. If
#'   `NULL`, no fixed number of points is specified. Default is `NULL`.
#'
#' @param seed Integer. Random seed set immediately before calling
#'   `get_train_data()`. This makes downsampling reproducible when `prop_down`
#'   or `n_points_per_plot` is used. Default is `123`.
#'
#' @param return_train_data Logical. If `TRUE`, return the training dataframe in
#'   the `ref_train` element. If `FALSE`, `ref_train` is set to `NULL` before
#'   returning the output. Default is `TRUE`.
#'
#' @param ... Additional arguments passed to internal functions. Arguments are
#'   automatically matched to the formal arguments of `magicGating()`,
#'   `get_train_data()`, or `magicTrain()`. This allows users to pass valid
#'   options for those functions without adding every possible argument to this
#'   wrapper.
#'
#' @return A list with four elements:
#'   \describe{
#'     \item{ref_model_info}{The trained model information returned by
#'       `magicTrain()`.}
#'     \item{ref_train}{The training dataframe returned by `get_train_data()`, or
#'       `NULL` if `return_train_data = FALSE`.}
#'     \item{list_out_final}{The complete merged template-gated data used to
#'       build the training dataframe.}
#'     \item{all_list_out_obj}{The raw template-gating outputs from each call to
#'       `magicGating()`. This object can be saved and later reused with
#'       `input_template_obj` and `use_existing_templates = TRUE`.}
#'   }
#'
#' @details
#' The function first checks that all selected samples exist in the input
#' `GatingSet` and that `samples_id` does not contain duplicated sample names.
#' Duplicate sample names are not allowed because sample names are used as list
#' names in downstream template objects.
#'
#' If `use_existing_templates = FALSE`, the function calls `magicGating()`
#' interactively once for each template gate. The labels assigned to the gates
#' are `"1"`, `"2"`, ..., up to `n_gates`.
#'
#' If `use_existing_templates = TRUE`, no interactive gating is performed. The
#' function reuses `input_template_obj$all_list_out_obj`, allowing the same
#' manually drawn templates to be used again for reproducible model training.
#'
#' After template generation or reuse, template-gated data are merged with
#' `merge_magicGating_labels()` when more than one gate is present. The merged
#' object is then passed to `get_train_data()`. The arguments `prop_down` and
#' `n_points_per_plot` control how many events are used for training. The seed is
#' set immediately before `get_train_data()` to make this downsampling
#' reproducible.
#'
#' The resulting training dataframe is passed to `magicTrain()` to build the
#' classification model.
#'
#' @keywords flowMagic template training machine-learning gating
#' @export
#'
#' @examples
#' \donttest{
#' # Train a template model by drawing two manual gates
#' train_out <- magic_template_train(
#'   gs_input = gs,
#'   parent_node = "root",
#'   samples_id = c("sample_01.fcs", "sample_02.fcs"),
#'   channel_x = "FSC-A",
#'   channel_y = "LIVE DEAD Blue-A",
#'   n_gates = 2,
#'   train_model = "rf",
#'   n_points_per_plot = 5000,
#'   seed = 123
#' )
#'
#' # Save the full template-training object for reproducibility
#' saveRDS(train_out, "template_training_root.rds")
#'
#' # Reuse the same templates later without redrawing gates
#' old_template <- readRDS("template_training_root.rds")
#'
#' train_out_reproduced <- magic_template_train(
#'   gs_input = gs,
#'   parent_node = "root",
#'   samples_id = c("sample_01.fcs", "sample_02.fcs"),
#'   channel_x = "FSC-A",
#'   channel_y = "LIVE DEAD Blue-A",
#'   input_template_obj = old_template,
#'   use_existing_templates = TRUE,
#'   train_model = "rf",
#'   n_points_per_plot = 5000,
#'   seed = 123
#' )
#'
#' # Use a proportion instead of a fixed number of points
#' train_out_prop <- magic_template_train(
#'   gs_input = gs,
#'   parent_node = "root",
#'   samples_id = c("sample_01.fcs", "sample_02.fcs"),
#'   channel_x = "FSC-A",
#'   channel_y = "LIVE DEAD Blue-A",
#'   n_gates = 2,
#'   prop_down = 0.25,
#'   seed = 123
#' )
#' }

magic_template_train <- function(gs_input,
                                 parent_node,
                                 samples_id,
                                 channel_x,
                                 channel_y,
                                 n_gates = 1,
                                 input_template_obj = NULL,
                                 use_existing_templates = FALSE,
                                 train_model = "rf",
                                 prop_down = NULL,
                                 n_points_per_plot = NULL,
                                 seed = 123,
                                 return_train_data = TRUE,
                                 ...) {
  
  # =========================================================================
  # 1. Start workflow
  # =========================================================================
  # This function creates or reuses manual template gates, converts the gated
  # events into a training dataframe, and trains a classifier.
  #
  # Reproducibility modes:
  #
  #   input_template_obj = NULL
  #     New templates are drawn interactively.
  #
  #   input_template_obj != NULL and use_existing_templates = FALSE
  #     Existing templates are loaded, and new templates are drawn interactively.
  #     Newly gated samples are added to the existing template object.
  #     Samples with matching names are replaced.
  #
  #   input_template_obj != NULL and use_existing_templates = TRUE
  #     No new interactive gating is performed. The saved templates in
  #     input_template_obj$all_list_out_obj are reused.
  
  message("$$$ Starting template training workflow $$$")
  message(sprintf("Parent node: %s", parent_node))
  message(sprintf("Channels: %s vs %s", channel_x, channel_y))
  message(sprintf("Number of template gates to draw per sample: %s", n_gates))
  message(sprintf("Number of samples selected for template gating: %s", length(samples_id)))
  
  # Check that selected samples exist in the GatingSet.
  missing_samples <- setdiff(samples_id, sampleNames(gs_input))
  
  if (length(missing_samples) > 0) {
    stop(
      "These samples are not present in sampleNames(gs_input): ",
      paste(missing_samples, collapse = ", ")
    )
  }
  
  # Check duplicated samples.
  # Duplicated sample names can break downstream list handling because sample
  # names are used as list names.
  if (anyDuplicated(samples_id)) {
    duplicated_samples <- unique(samples_id[duplicated(samples_id)])
    stop(
      "samples_id contains duplicated sample names: ",
      paste(duplicated_samples, collapse = ", "),
      ". Please provide each sample only once."
    )
  }
  
  
  # =========================================================================
  # 2. Collect optional arguments from ...
  # =========================================================================
  # Users may pass additional arguments for the internal functions through ...
  #
  # These arguments may belong to:
  #   - magicGating()
  #   - magicPlot_template(), through magicGating()
  #   - get_train_data()
  #   - magicTrain()
  #
  # We split them automatically by checking the formal arguments of each
  # function. This avoids passing arguments to functions that do not use them.
  
  args_all <- list(...)
  
  get_matching_args <- function(args_list, fun) {
    
    fun_args <- names(formals(fun))
    
    args_list[names(args_list) %in% fun_args]
  }
  
  # Arguments for magicGating().
  #
  # Important:
  # magicGating() has ... and internally forwards some of those arguments to
  # magicPlot_template().
  #
  # Therefore, for magicGating(), we keep arguments that match either:
  #   1. magicGating() formal arguments
  #   2. magicPlot_template() formal arguments
  #
  # This is needed so arguments like size_points are not removed here.
  args_magicGating_names <- unique(c(
    names(formals(flowMagic::magicGating)),
    names(formals(flowMagic::magicPlot_template))
  ))
  
  args_magicGating <- args_all[
    names(args_all) %in% args_magicGating_names
  ]
  
  args_get_train_data <- get_matching_args(
    args_list = args_all,
    fun = flowMagic::get_train_data
  )
  
  args_magicTrain <- get_matching_args(
    args_list = args_all,
    fun = flowMagic::magicTrain
  )
  
  # Warn the user if some arguments in ... were not used by any internal
  # function. This helps catch spelling mistakes.
  used_args <- unique(c(
    names(args_magicGating),
    names(args_get_train_data),
    names(args_magicTrain)
  ))
  
  unused_args <- setdiff(names(args_all), used_args)
  
  if (length(unused_args) > 0) {
    warning(
      "These arguments in ... were not matched to magicGating(), ",
      "get_train_data(), or magicTrain(): ",
      paste(unused_args, collapse = ", ")
    )
  }
  
  
  # =========================================================================
  # 3. Prepare or reuse template object
  # =========================================================================
  # This section decides whether templates are drawn interactively or reused from
  # input_template_obj.
  
  if (!is.null(input_template_obj) && use_existing_templates == TRUE) {
    
    message("Existing template object provided.")
    message("use_existing_templates = TRUE, so interactive gating will be skipped.")
    message("Templates stored in input_template_obj$all_list_out_obj will be reused.")
    
    all_list_out_obj <- input_template_obj$all_list_out_obj
    
    if (is.null(all_list_out_obj)) {
      stop("input_template_obj does not contain all_list_out_obj.")
    }
    
    if (length(all_list_out_obj) == 0) {
      stop("input_template_obj$all_list_out_obj is empty.")
    }
    
    # When reusing existing templates, update n_gates so that downstream
    # messages reflect the actual number of stored gates.
    n_gates <- length(all_list_out_obj)
    
    message(sprintf("Number of existing template gates reused: %s", n_gates))
    
  } else {
    
    if (is.null(input_template_obj)) {
      
      message("No existing template object provided. A new template set will be created.")
      all_list_out_obj <- list()
      
    } else {
      
      message("Existing template object provided.")
      message("use_existing_templates = FALSE, so new interactive templates will be drawn.")
      message("Existing templates will be used as the starting template object.")
      message("Newly gated samples will be added to the existing templates.")
      message("Samples with matching names will be replaced.")
      
      all_list_out_obj <- input_template_obj$all_list_out_obj
      
      if (is.null(all_list_out_obj)) {
        all_list_out_obj <- list()
      }
    }
    
    
    # =========================================================================
    # 4. Perform interactive template gating
    # =========================================================================
    # Each template gate is drawn by calling magicGating().
    #
    # label_pol is set automatically:
    #   gate 1 -> label "1"
    #   gate 2 -> label "2"
    #   etc.
    #
    # Additional arguments that belong to magicGating() are passed through
    # args_magicGating.
    
    message("$$$ Perform interactive template gating $$$")
    
    for (n in seq_len(n_gates)) {
      
      label_current <- as.character(n)
      
      message(sprintf("---- Gating template/gate %s of %s ----", n, n_gates))
      message(sprintf("Current label assigned to this gate: %s", label_current))
      
      args_gating_current <- c(
        list(
          fs = gs_input,
          gs_node = parent_node,
          sample_id = samples_id,
          channel_x = channel_x,
          channel_y = channel_y,
          label_pol = label_current
        ),
        args_magicGating
      )
      
      new_list_out_obj <- do.call(
        flowMagic::magicGating,
        args_gating_current
      )
      
      # ---------------------------------------------------------------------
      # Add new template samples to an existing template object
      # ---------------------------------------------------------------------
      # If all_list_out_obj[[n]] already exists, this means we are adding new
      # manually gated samples to an existing template gate.
      #
      # In that case:
      #   - old samples are kept,
      #   - new samples are added,
      #   - samples with matching names are replaced.
      #
      # If all_list_out_obj[[n]] does not exist yet, store the newly created gate
      # object directly.
      
      if (!is.null(all_list_out_obj[[n]])) {
        
        old_gated_data <- all_list_out_obj[[n]]$list_gated_data
        new_gated_data <- new_list_out_obj$list_gated_data
        
        old_poly_gates <- all_list_out_obj[[n]]$list_poly_gates
        new_poly_gates <- new_list_out_obj$list_poly_gates
        
        old_gated_data[names(new_gated_data)] <- new_gated_data
        old_poly_gates[names(new_poly_gates)] <- new_poly_gates
        
        all_list_out_obj[[n]]$list_gated_data <- old_gated_data
        all_list_out_obj[[n]]$list_poly_gates <- old_poly_gates
        
      } else {
        
        all_list_out_obj[[n]] <- new_list_out_obj
      }
    }
  }
  
  
  # =========================================================================
  # 5. Merge template-gated data
  # =========================================================================
  # magicGating() returns one object per template gate.
  #
  # If there is only one template gate:
  #   use the list_gated_data directly.
  #
  # If there is more than one template gate:
  #   merge labels across template gates using merge_magicGating_labels().
  
  message("$$$ Merge template gated data $$$")
  
  if (length(all_list_out_obj) == 1) {
    
    list_out_final <- all_list_out_obj[[1]]$list_gated_data
    
  } else {
    
    all_gated_data <- lapply(
      all_list_out_obj,
      function(x) x$list_gated_data
    )
    
    list_out_final <- Reduce(
      function(x, y) {
        flowMagic::merge_magicGating_labels(
          list_out_1 = x,
          list_out_2 = y,
          from_gs = TRUE
        )
      },
      all_gated_data
    )
  }
  
  
  # =========================================================================
  # 6. Build training data
  # =========================================================================
  # get_train_data() converts the manually gated data into the final training
  # dataframe used by magicTrain().
  #
  # prop_down and n_points_per_plot are explicit arguments of this wrapper
  # because they are important for controlling training size.
  #
  # seed is set immediately before get_train_data() so that downsampling is
  # reproducible when prop_down or n_points_per_plot is used.
  #
  # Additional arguments that belong to get_train_data() are passed through
  # args_get_train_data.
  
  message("$$$ Build training data $$$")
  message(sprintf("Set seed before get_train_data(): %s", seed))
  set.seed(seed)
  
  args_train_data <- c(
    list(
      paths_file = list_out_final,
      prop_down = prop_down,
      n_points_per_plot = n_points_per_plot
    ),
    args_get_train_data
  )
  
  ref_train <- do.call(
    flowMagic::get_train_data,
    args_train_data
  )
  
  
  # =========================================================================
  # 7. Train model
  # =========================================================================
  # magicTrain() trains the classifier using the dataframe produced by
  # get_train_data().
  #
  # Additional arguments that belong to magicTrain() are passed through
  # args_magicTrain.
  
  message("$$$ Train template model $$$")
  
  args_train_model <- c(
    list(
      df_train = ref_train,
      train_model = train_model
    ),
    args_magicTrain
  )
  
  ref_model_info <- do.call(
    flowMagic::magicTrain,
    args_train_model
  )
  
  
  # =========================================================================
  # 8. Return output
  # =========================================================================
  # The returned object contains:
  #
  #   ref_model_info:
  #     trained model object and model metadata returned by magicTrain().
  #
  #   ref_train:
  #     training dataframe used by magicTrain().
  #
  #   list_out_final:
  #     merged manually gated data used to build the training dataframe.
  #
  #   all_list_out_obj:
  #     original template-gating objects returned by magicGating().
  #
  # all_list_out_obj is important for reproducibility because it allows the same
  # templates to be reused later with:
  #
  #   input_template_obj = saved_object
  #   use_existing_templates = TRUE
  
  message("$$$ Template training workflow completed $$$")
  
  out <- list(
    ref_model_info = ref_model_info,
    ref_train = ref_train,
    list_out_final = list_out_final,
    all_list_out_obj = all_list_out_obj
  )
  
  if (return_train_data == FALSE) {
    out$ref_train <- NULL
  }
  
  return(out)
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


