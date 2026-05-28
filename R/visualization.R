#' magicPlot
#'
#' function to generate the scatter plot with colored density of the events.
#' @param df Dataframe of bivariate markers expression (with labels if gates to plot).
#' @param type Type of plot to generated. "dens"= plot with bivariate density. "ML"= plot with color layer based on machine learning class assignments 
#' @param polygons_coords_list list of gates coordinates. Needed if labels not included in df. Default to NULL.
#' @param show_legend Show legend if type="ML". Default to True.
#' @param size_axis_text Size of axis ticks labels. Default to 18.
#' @param size_title_x Size of x axis label title. Default to 20.
#' @param size_title_y Size of y axis label title. Default to 20.
#' @param treat_0_as_gate Treat 0 label as gate. Defaul to False (0 label is background)
#' @param x_lab Label of x axis. Default to NULL (inherited from input dataframe by default).
#' @param y_lab Label of y axis. Default to NULL (inherited from input dataframe by default).
#' @param gates_to_plot Select labels to plot.
#' @param apply_manual_scale Apply predifined scale of colors. Default to True.
#' @param size_points Size of points in scatter plot.
#' @param concavity_val Concavity value. Default to 20. Higher value, less jagged boundaries
#' @param aspect_ratio Aspect ratio value. If = 1, y and x axis have equal distance between the ticks labels. Default to NULL.
#' @param x_lim1 Minimum limit x axis. Default to NULL.
#' @param x_lim2 Max limit x axis. Default to NULL.
#' @param y_lim1 Minimum limit y axis. Default to NULL.
#' @param y_lim2 Max limit y axis. Default to NULL.
#' @param add_labels add polygon labels. Default to FALSE.
#' @param map_label_polygon map of polygon labels (to assign custom labels). Default to NULL.
#' @param size_pol_name size polygon labels. Default to 6.
#' @param show_marginals show 1D density next to axis. Default to False.
#' @return Object of class ggplot
#' @keywords flowMagic
#' @export
#' @examples 
#' \donttest{magicPlot()}

magicPlot<-function(df, type = "dens", polygons_coords_list = NULL, show_legend = T, 
          size_axis_text = 18, size_title_x = 20, size_title_y = 20, 
          treat_0_as_gate = F, x_lab = NULL, y_lab = NULL, gates_to_plot = NULL, 
          apply_manual_scale = F, hull_only = F, size_points = 1, concavity_val = 20, 
          aspect_ratio = NULL, x_lim1 = NULL, x_lim2 = NULL, y_lim1 = NULL, 
          y_lim2 = NULL,add_labels=F,map_label_polygon=NULL,size_pol_name=6,show_marginals=F,...){
  # copy of attributes
  manual_polygon <- attr(df, "manual_polygon")
  gate_source <- attr(df, "gate_source")
  # capture all extra arguments
  args_list <- list(...)
  # Arguments for get_hull_all_gates
  args_get_hull <- args_list[names(args_list) %in% names(formals(get_hull_all_gates))]

  # save original column names
  original_col_names<-colnames(df)
  # ----------- plot with no gates --------------
  if (type == "no_gate" || ncol(df) == 2) {
    colPalette <- colorRampPalette(c("blue", "turquoise", "green", "yellow", "orange", "red"))
    col <- densCols(df[, c(1, 2)], colramp = colPalette, nbin = 200)
    
    df_plot <- data.frame(x = df[[1]], y = df[[2]], col = col)
    
    p <- ggplot(df_plot, aes(x = x, y = y)) + geom_point(color = df_plot$col, size = size_points, shape = 16)
    
    if(is.null(x_lab)==T && is.null(y_lab)==T){
      p <- p + xlab(original_col_names[1]) + ylab(original_col_names[2])
    }else if(is.null(x_lab)==F || is.null(y_lab)==F){
      p <- p + xlab(x_lab) + ylab(y_lab)
    }
    
    
    p <- p + theme(legend.position = "none")
    p <- p + theme(axis.text = element_text(size = size_axis_text), 
                   axis.title.x = element_text(size = size_title_x, 
                                               face = "bold"), axis.title.y = element_text(size = size_title_y, 
                                                                                           face = "bold"))
    p <- p + theme(panel.grid.major = element_blank(), 
                   panel.grid.minor = element_blank(), panel.background = element_blank(), 
                   axis.line = element_line(colour = "black"))
    # Apply limits if specified
    if (!is.null(x_lim1) && !is.null(x_lim2)) {
      p <- p + xlim(x_lim1, x_lim2)
    }
    if (!is.null(y_lim1) && !is.null(y_lim2)) {
      p <- p + ylim(y_lim1, y_lim2)
    }
    
    # Add marginal densities
    if (show_marginals == TRUE) {
      p <- ggExtra::ggMarginal(p, type = "density", margins = "both", 
                               groupColour = FALSE, groupFill = F, fill="gray")
    }
    
    return(p)
  }
  colnames(df) <- c("x", "y", "classes")
  # performing some modifications on the classes column based on user choice
  df$classes <- as.character(df$classes)
  if (treat_0_as_gate == T) {
    df$classes <- as.numeric(df$classes)
    inds <- which(df$classes == 0)
    df$classes[inds] <- max(unique(df$classes)) + 1
    df$classes <- as.character(df$classes)
  }
  if (is.null(gates_to_plot) == F) {
    check_classes <- df$classes %in% gates_to_plot
    inds <- which(check_classes == T)
    df$classes[-inds] <- "0"
  }
  # ----------- plot with gates --------------
  if(type == "dens"){
    # Plotting with density layer ------

    colPalette <- colorRampPalette(c("blue", "turquoise", 
                                     "green", "yellow", "orange", "red"))
    col <- densCols(df[, c(1, 2)], colramp = colPalette, 
                    nbin = 200)
    df <- cbind(df, col)

    # restore attributes lost with cbind (because cbind make a new object dropping the old attributes)
    attr(df, "manual_polygon") <- manual_polygon
    attr(df, "gate_source") <- gate_source

    df$col <- as.character(df$col)
    if (is.null(polygons_coords_list) == T) {
      list_df_hull<- do.call(get_hull_all_gates,c(list(gated_df=df,concavity_val = concavity_val), args_get_hull))
      vec <- names(list_df_hull)
      inds <- which(vec == "0")
      if (length(inds) != 0) {
        list_df_hull <- list_df_hull[-inds]
      }
      df_hull <- do.call(rbind, list_df_hull)
    }
    else {
      df_hull <- do.call(rbind, polygons_coords_list)
    }
    magicggplot <- ggplot(data = df, aes(x = x, y = y, col = factor(1:nrow(df)))) + 
      geom_point(size = size_points, shape = 16)
    magicggplot <- magicggplot + scale_colour_manual(values = col)
    magicggplot <- magicggplot + theme(legend.position = "none")
    magicggplot <- magicggplot + theme(axis.text = element_text(size = size_axis_text), 
                                       axis.title.x = element_text(size = size_title_x, 
                                                                   face = "bold"), axis.title.y = element_text(size = size_title_y, 
                                                                                                               face = "bold"))
    magicggplot <- magicggplot + theme(panel.grid.major = element_blank(), 
                                       panel.grid.minor = element_blank(), panel.background = element_blank(), 
                                       axis.line = element_line(colour = "black"))
    magicggplot <- magicggplot + geom_polygon(data = df_hull, 
                                              aes(x = x, y = y, group = group_gate), fill = NA, 
                                              color = "black", linewidth = 1)
  }else if (type == "ML") {
    # Plotting with ML class color layer ------

    if(apply_manual_scale==T){
      # check if there are real text labels (like full strings,not just numbers-like strings)
      df<-convert_to_integers_chr(df=df)
      df$classes <- as.factor(df$classes)
      magicggplot <- ggplot(data = df, aes(x = x, y = y)) + 
        geom_point(aes(color = classes), size = size_points, 
                   shape = 16)
      magicggplot <- magicggplot + scale_color_manual(values = c(`0` = "black", 
                                                                 `1` = "red", `2` = "blue", `3` = "green", `4` = "yellow", 
                                                                 `5` = "purple", `6` = "orange", `7` = "pink", 
                                                                 `8` = "brown", `9` = "gray"))
    }else{
      df$classes <- as.factor(df$classes)
      all_levels <- levels(df$classes)
      # Rename "0" to "background" in a copy of the factor for plotting
      df$classes_legend <- df$classes
      levels(df$classes_legend)[levels(df$classes_legend) == "0"] <- "background"
      
      # Now assign colors, Separating "0|background" from other classes still keyed to the original factor levels, but color vector keys to new labels
      other_levels <- setdiff(levels(df$classes_legend), "background")
      n_other <- length(other_levels)
      if(n_other<3){
        n_min_colors<-3
        # Assign colors
        colors_other_all <- RColorBrewer::brewer.pal(n_min_colors, "Set2")
        colors_other <- colors_other_all[1:n_other]
        names(colors_other) <- other_levels
      }else if(n_other >= 3 & n_other <=8){
        # Assign colors
        colors_other <- RColorBrewer::brewer.pal(n_other, "Set2")
        names(colors_other) <- other_levels
      }else if(n_other > 8 ){
        stop("Set2 supports only up to 8 colors for classes other than '0'")
        
      }
      
      # Assign black for "0|background"
      colors_all <- c("background" = "black", colors_other)
      magicggplot <- ggplot(data = df, aes(x = x, y = y)) + 
        geom_point(aes(color = classes_legend), size = size_points, 
                   shape = 16)
      magicggplot <- magicggplot + scale_color_manual(values = colors_all)
    }
    
    magicggplot <- magicggplot + theme(axis.text = element_text(size = size_axis_text), 
                                       axis.title.x = element_text(size = size_title_x, 
                                                                   face = "bold"), axis.title.y = element_text(size = size_title_y, 
                                                                                                               face = "bold"))
    magicggplot <- magicggplot + theme(panel.grid.major = element_blank(), 
                                       panel.grid.minor = element_blank(), panel.background = element_blank(), 
                                       axis.line = element_line(colour = "black"))
    magicggplot <- magicggplot + theme(legend.key.size = unit(1, 
                                                              "cm"), legend.title = element_text(size = 20), legend.text = element_text(size = 20))
    magicggplot <- magicggplot + guides(color = guide_legend(override.aes = list(size = 10))) + 
      labs(color = "Assignment")
    if (show_legend == F) {
      magicggplot <- magicggplot + theme(legend.position = "none")
    }
  }
  #----- additional options ------------
  if(is.null(x_lab)==T && is.null(y_lab)==T){
    magicggplot <- magicggplot + xlab(original_col_names[1]) + ylab(original_col_names[2])
  }else if(is.null(x_lab)==F || is.null(y_lab)==F){
    magicggplot <- magicggplot + xlab(x_lab) + ylab(y_lab)
  }
  if (is.null(aspect_ratio) == F) {
    magicggplot <- magicggplot + coord_fixed(ratio = aspect_ratio) + 
      theme(plot.margin = unit(c(0, 0, 0, 0), "pt"))
  }
  # Apply limits if specified
  if (!is.null(x_lim1) && !is.null(x_lim2)) {
    magicggplot <- magicggplot + xlim(x_lim1, x_lim2)
  }
  if (!is.null(y_lim1) && !is.null(y_lim2)) {
    magicggplot <- magicggplot + ylim(y_lim1, y_lim2)
  }
  if(add_labels==T & type=="dens"){
    if(is.null(map_label_polygon)==F){
      df_hull$group_gate<-map_label_polygon[df_hull$group_gate]
    }
    label_coords <- do.call(rbind, lapply(split(df_hull, df_hull$group_gate), function(sub) {
      data.frame(
        group = unique(sub$group),
        label_x = mean(sub$x),
        label_y = mean(sub$y)
      )
    }))
    magicggplot<- magicggplot + geom_text(data = label_coords, aes(x = label_x, y = label_y, label = group), 
                                          fontface = "bold",inherit.aes = FALSE,size = size_pol_name)
  }
  
  # Add marginal densities
  if (show_marginals == TRUE) {
    magicggplot <- ggExtra::ggMarginal(magicggplot, type = "density", margins = "both", 
                                       groupColour = FALSE, groupFill = F, fill="gray")
  }
  return(magicggplot)
}




#' magic_plot_wrap
#'
#' function to wrap all plots in one plot
#' @param list_gated_data List of dataframes composed by three columns (could be generated by the get_list_df_gated_plots() function)
#' @param n_col_wrap number of columns in the wrapped plot. Default to 3.
#' @param size_title size of title of each plot. Default to 10.
#' @return wrapped plot
#' @keywords flowMagic
#' @export
#' @examples 
#' \donttest{magic_plot_wrap()}

magic_plot_wrap<-function(list_gated_data,n_col_wrap=3,size_title=10,...){
  all_names<-names(list_gated_data)
  plot_list <- lapply(seq_along(list_gated_data),function(i){
    df <- list_gated_data[[i]]
    if("df_test_original" %in% names(df)){
      df<-df$df_test_original
    }
    p <- magicPlot(df,...) + ggtitle(all_names[i]) + theme(plot.title = element_text(size = size_title))
    return(p)
  })
  
  # Combine into a grid layout, e.g., 3 columns
  combined_plot <- patchwork::wrap_plots(plot_list, ncol = n_col_wrap)
  
  return(combined_plot)
}


#' magicplot_3D
#'
#' function to make a 3D plot of data.
#' @param df input dataframe composed of at least three columns. Four columns if gate needs to be plot.
#' @param class_col Gate to plot? Default to False.
#' @param x_lab label of the x axis
#' @param y_lab Label of the y axis
#' @param z_lab Label of the z axis
#' @param type Assuming class_col==T, If type==ML, Generate a plot colored based on the gate assignment. If type=="mesh", generate a 3d polygon gate.
#' @param size_p size of scatter plot points. Default to 1.
#' @return plotly plot
#' @keywords flowMagic
#' @export
#' @examples 
#' \donttest{magicplot_3D()}

magicplot_3D<-function(df,class_col=F,x_lab="x",y_lab="y",z_lab="z",type="ML",size_p=1){
  if (!requireNamespace("plotly", quietly = TRUE)) {
    stop("The 'plotly' package is required for this function. Please install it.")
  }
  if (!requireNamespace("magrittr", quietly = TRUE)) {
    stop("The 'magrittr' package is required for piping (%>%). Please install it.")
  }

  if (!class_col && ncol(df) != 3) {
    stop("The input df must have 3 columns if class_col is FALSE")
  }
  if (class_col && ncol(df) != 4) {
    stop("The input df must have 4 columns if class_col is TRUE")
  }

  if (!class_col) {
    plotly_plot <- plotly::plot_ly(
      x = df[,1], y = df[,2], z = df[,3],
      type = "scatter3d", mode = "markers",
      marker = list(size = size_p)
    ) %>%
      plotly::layout(
        scene = list(
          xaxis = list(title = x_lab),
          yaxis = list(title = y_lab),
          zaxis = list(title = z_lab)
        ),
        showlegend = FALSE
      )
  } else if (class_col && type == "ML") {
    df[,4] <- as.factor(df[,4])
    plotly_plot <- plotly::plot_ly(
      x = df[,1], y = df[,2], z = df[,3],
      type = "scatter3d", mode = "markers",
      marker = list(size = size_p),
      color = df[,4],
      colors = c("blue", "red")
    ) %>%
      plotly::layout(
        scene = list(
          xaxis = list(title = x_lab),
          yaxis = list(title = y_lab),
          zaxis = list(title = z_lab)
        ),
        showlegend = FALSE
      )
  } else if (class_col && type == "mesh") {
    if (!requireNamespace("geometry", quietly = TRUE)) {
      stop("The 'geometry' package is required for mesh plotting. Please install it.")
    }

    df_class <- df[df[,4] == 1, ]
    coords <- cbind(df_class[,1], df_class[,2], df_class[,3])
    hull <- geometry::convhulln(coords)
    red_rgb <- rep("rgb(255,0,0)", nrow(coords))

    plotly_plot <- plotly::plot_ly() %>%
      plotly::add_trace(
        type = "scatter3d",
        mode = "markers",
        x = df[,1], y = df[,2], z = df[,3],
        color = df[,4],
        marker = list(size = size_p),
        colors = c("blue", "black"),
        name = "Parent"
      ) %>%
      plotly::add_trace(
        type = "mesh3d",
        x = coords[,1], y = coords[,2], z = coords[,3],
        i = hull[,1] - 1, j = hull[,2] - 1, k = hull[,3] - 1,
        vertexcolor = red_rgb,
        opacity = 0.3,
        color = I("red"),
        flatshading = TRUE,
        lighting = list(ambient = 1, diffuse = 0, specular = 0, roughness = 0, fresnel = 0),
        lightposition = list(x = 0, y = 0, z = 100),
        name = "Gated pop"
      ) %>%
      plotly::layout(
        scene = list(
          xaxis = list(title = x_lab),
          yaxis = list(title = y_lab),
          zaxis = list(title = z_lab)
        ),
        showlegend = FALSE
      )
  }

  return(plotly_plot)
}

#' magicPlot_template
#'
#' function to interactively generate a template based on bivariate ungated data.
#' @param df input dataframe composed of two columns reporting markers expression (first column = marker 1, second column = marker 2).
#' @param size_points Size of points in scatter plot. Default to 1.
#' @return Dataframe
#' @keywords flowMagic
#' @export
#' @examples 
#' \donttest{magicPlot_template()}

magicPlot_template<-function(df,size_points=1){
  colPalette <- colorRampPalette(c("blue", "turquoise", "green", "yellow", "orange", "red"))
  col <- densCols(df[,1], df[,2], colramp = colPalette, nbin = 200)
  
  par(mar = c(5, 5, 2, 2))  # Set margins
  names<-colnames(df)
  plot(df[,1], df[,2],
       col = col,
       pch = 16,
       cex = size_points,
       xlab = names[1],
       ylab = names[2],
       bty = "n")
  
  # Use locator() to draw polygon interactively
  cat("Click to draw a polygon point-by-point. Press ESC when done.\n")
  pts <- graphics::locator(type = "o")  # type = "n" means collect points only, not draw
  
  # Convert to dataframe
  polygon_df <- data.frame(x = pts$x, y = pts$y)
  
  return(polygon_df)
  
}

#' magicPlot_fs
#'
#' function to plot FCM data directly from objects of class flowSet or flowFrame. The data is shown in a density scatter plot.
#' @param fs object of class flowSet or flowFrame.
#' @param sample_id Name of the sample data to plot. No effect if fs is an object of class flowFrame.
#' @param channel_x Name of the channel to plot (x-axis).
#' @param channel_y Name of the channel to plot (y-axis).
#' @return Object of class ggplot
#' @keywords flowMagic
#' @export
#' @examples 
#' \donttest{magicPlot_fs()}

magicPlot_fs<-function(fs,sample_id,channel_x,channel_y,...){
  if(class(fs)=="flowFrame"){
    ff<-fs
  }else if(class(fs)=="flowSet"){
    ff<-fs[[sample_id]]
  }else{
    stop("Error input: fs needs to be an object of either flowSet or flowFrame class")
  }
  expr_matrix_ff <- exprs(ff) # get expression matrix of selected flowFrame
  df_exprs<-as.data.frame(expr_matrix_ff) # convert to dataframes.
  df_exprs_selected_channels<-df_exprs[,c(channel_x,channel_y)]
  out<-flowMagic::magicPlot(df = df_exprs_selected_channels,...)
  return(out)
}

#' magicPlot_gs_node
#'
#' Plot selected channels from one GatingSet node.
#'
#' This function extracts one sample and one population node from a `GatingSet`
#' as a `flowFrame`, then plots two selected channels using `magicPlot_fs()`.
#' Unlike `magic_plot_gs_gates()`, this function does not overlay multiple gate
#' labels. It simply plots the events contained in the selected node.
#'
#' @param gs_input A `GatingSet` object.
#' @param node_name Character. Name of the population node to extract, such as
#'   `"root"`, `"Singlets"`, or `"Live_cells"`.
#' @param sample_id Sample identifier used to select one sample from the
#'   `GatingSet`. This can be a sample name or an index, depending on
#'   `get_flowframe_from_gs()`.
#' @param channel_x Character. Name of the channel to plot on the x-axis.
#' @param channel_y Character. Name of the channel to plot on the y-axis.
#' @param ... Additional arguments passed to `magicPlot_fs()`.
#'
#' @return A plot object returned by `magicPlot_fs()`.
#'
#' @keywords flowMagic plotting GatingSet
#' @export
#'
#' @examples
#' \donttest{
#' magicPlot_gs_node(
#'   gs_input = gs,
#'   node_name = "Live_cells",
#'   sample_id = "sample_01.fcs",
#'   channel_x = "FSC-A",
#'   channel_y = "LIVE DEAD Blue-A"
#' )
#' }

magicPlot_gs_node <- function(gs_input, node_name, sample_id, channel_x, channel_y, ...) {
  
  # Extract the selected sample and node from the GatingSet.
  # The returned object is a flowFrame containing the events in node_name.
  ff <- flowMagic::get_flowframe_from_gs(
    gs = gs_input,
    node_name = node_name,
    sample_id = sample_id
  )
  
  # Plot the selected x/y channels from the extracted flowFrame.
  # Additional plotting options are passed to magicPlot_fs() through ... .
  out <- flowMagic::magicPlot_fs(
    fs = ff,
    channel_x = channel_x,
    channel_y = channel_y,
    ...
  )
  
  return(out)
}

#' magicPlot_gs_gates
#'
#' Plot one sample from a GatingSet with one or more gated populations overlaid.
#'
#' This function extracts one or more gated populations from a single sample in a
#' `GatingSet` and visualizes them using `magicPlot()`. If multiple gate names are
#' supplied, the extracted gated data are merged so the selected populations are
#' displayed together on the same bivariate plot. Each population is labeled using
#' its node name.
#'
#' @param gs A `GatingSet` object.
#' @param sample_id Character. Name of the sample to plot. Must be present in
#'   `sampleNames(gs)`.
#' @param gate_names Character vector. One or more GatingSet node names to
#'   extract and display.
#' @param size_points Numeric. Point size passed to `magicPlot()`. Default is
#'   `0.5`.
#' @param concavity_val Numeric. Concavity parameter available for downstream
#'   polygon plotting. Default is `50`.
#' @param ... Additional arguments passed to `magicPlot()`, such as `type`,
#'   `size_axis_text`, `size_title_x`, or `size_title_y`.
#'
#' @return A `ggplot` object showing the selected sample and gated populations.
#'
#' @keywords flowMagic plotting GatingSet
#' @export
#'
#' @examples
#' \donttest{
#' # Plot one gate
#' p1 <- magicPlot_gs_gates(
#'   gs = gs,
#'   sample_id = "sample_01.fcs",
#'   gate_names = "Live_cells",
#'   size_points = 0.5
#' )
#'
#' # Plot multiple gates together
#' p2 <- magicPlot_gs_gates(
#'   gs = gs,
#'   sample_id = "sample_01.fcs",
#'   gate_names = c("Live_cells", "Dead_cells"),
#'   type = "ML",
#'   size_points = 0.5
#' )
#' }


magicPlot_gs_gates <- function(gs,
                                 sample_id,
                                 gate_names,
                                 size_points = 0.5,
                                 concavity_val=50,
                                 ...) {
  
  message("$$$ Plot GatingSet sample $$$")
  message(sprintf("Sample: %s", sample_id))
  message(sprintf("Gates: %s", paste(gate_names, collapse = ", ")))
  
  if (!(sample_id %in% sampleNames(gs))) {
    stop("sample_id is not present in sampleNames(gs): ", sample_id)
  }
  
  list_gates <- list()
  
  for (i in seq_along(gate_names)) {
    
    current_gate <- gate_names[i]
    current_label <- current_gate
    
    message(sprintf(
      "Extracting gate %s and assigning label %s",
      current_gate,
      current_label
    ))
    
    list_gated_i <- get_list_df_gated_plots(
      gs = gs[sample_id],
      gate_name = current_gate,
      label_pop = current_label
    )
    
    list_gates[[current_gate]] <- list_gated_i
  }
  
  if (length(list_gates) == 1) {
    
    list_gated_data <- list_gates[[1]]
    
  } else {
    
    message("Merging selected gates into one plot.")
    
    list_gated_data <- Reduce(
      f = function(x, y) {
        merge_magicGating_labels(
          list_out_1 = x,
          list_out_2 = y,
          gated_data_only  = TRUE
        )
      },
      x = list_gates
    )
  }
  
  df_plot <- list_gated_data[[sample_id]]
  
  message(sprintf("Plot selected sample: %s", sample_id))
  
  p <- magicPlot(
    df = df_plot,
    size_points = size_points,
    concavity_val=concavity_val,
    ...
  )
  
  p <- p + ggplot2::ggtitle(sample_id)
  
  return(p)
}


#' magicPlot_gs_hierarchy
#'
#' Plot all gates from one GatingSet sample as a hierarchy-style figure.
#'
#' This function generates a wrapped figure containing one bivariate plot for
#' each gating step in a `GatingSet` hierarchy. For each subplot, the title
#' indicates the parent population and the gate labels indicate the child
#' population or populations shown on that plot.
#'
#' Gates that share the same bivariate dimensions are automatically grouped into
#' the same subplot. For example, if `NK_cells`, `B_cells`, and `T_cells` are
#' sibling gates drawn on the same two channels, they are displayed together in
#' one plot. Gates that do not share dimensions with sibling gates are plotted
#' individually.
#'
#' The hierarchy order is determined from `get_hierarchy_all_pops()`, and the
#' final plots are combined using `patchwork::wrap_plots()`.
#'
#' @param gs A `GatingSet` object.
#'
#' @param sample_id Character. Name of the sample to plot. The value must be
#'   present in `sampleNames(gs)`.
#'
#' @param n_col_wrap Integer. Number of columns used to arrange the wrapped
#'   hierarchy figure. Default is `4`.
#'
#' @param size_points Numeric. Point size passed to `magicPlot_gs_gates()` and
#'   ultimately to `magicPlot()`. Default is `0.5`.
#'
#' @param size_title Numeric or `NULL`. Size of the title for each subplot. Each
#'   subplot title indicates the parent population, for example
#'   `"Parent pop: CD45+"`. If `NULL` and `auto_size = TRUE`, the value is
#'   computed automatically from the grid size. If `NULL` and
#'   `auto_size = FALSE`, the default value is `10`.
#'
#' @param concavity_val Numeric. Concavity value used when polygon boundaries
#'   are estimated from gated events. This value is passed to
#'   `magicPlot_gs_gates()` and then to `magicPlot()`. Larger values generally
#'   produce smoother polygon boundaries. Default is `50`.
#'
#' @param add_labels Logical. If `TRUE`, gate names are printed inside the
#'   polygon gates. This is passed to `magicPlot()`. Labels are only shown for
#'   plotting modes that support polygon labels, such as `type = "dens"`.
#'   Default is `TRUE`.
#'
#' @param size_pol_name Numeric or `NULL`. Size of the gate labels printed inside
#'   polygon gates. If `NULL` and `auto_size = TRUE`, the value is computed
#'   automatically from the grid size. If `NULL` and `auto_size = FALSE`, the
#'   default value is `4`.
#'
#' @param auto_size Logical. If `TRUE`, automatically estimates text sizes for
#'   axis tick labels, axis titles, polygon labels, and subplot titles based on
#'   the number of rows and columns in the wrapped figure. User-supplied values
#'   passed through `...`, such as `size_axis_text`, `size_title_x`, or
#'   `size_title_y`, are respected and not overwritten. Default is `TRUE`.
#'
#' @param return_plot_list Logical. If `FALSE`, return only the wrapped plot. If
#'   `TRUE`, return a list containing the wrapped plot and intermediate objects
#'   used to build it. Default is `FALSE`.
#'
#' @param path_output Character or `NULL`. Full file path where the wrapped
#'   hierarchy plot should be saved. If `NULL`, no file is exported. The file
#'   format is inferred by `ggplot2::ggsave()` from the file extension, for
#'   example `.tiff`, `.png`, or `.pdf`. Default is `NULL`.
#'
#' @param plot_width Numeric. Width, in inches, assigned to each individual
#'   subplot when exporting the full wrapped figure. The final exported width is
#'   calculated as `plot_width * number_of_columns`. Default is `4`.
#'
#' @param plot_height Numeric. Height, in inches, assigned to each individual
#'   subplot when exporting the full wrapped figure. The final exported height is
#'   calculated as `plot_height * number_of_rows`. Default is `4`.
#'
#' @param dpi Numeric. Resolution used when exporting raster image formats such
#'   as `.png` or `.tiff`. Default is `300`.
#'
#' @param ... Additional plotting arguments passed to `magicPlot_gs_gates()` and
#'   then to `magicPlot()`. Common examples include `type`, `size_axis_text`,
#'   `size_title_x`, `size_title_y`, `show_legend`, `x_lim1`, `x_lim2`,
#'   `y_lim1`, and `y_lim2`.
#'
#' @return If `return_plot_list = FALSE`, a wrapped `patchwork` plot object. If
#'   `return_plot_list = TRUE`, a list with the following elements:
#'   \describe{
#'     \item{plot}{The final wrapped hierarchy plot.}
#'     \item{plot_list}{A named list of individual `ggplot` objects before wrapping.}
#'     \item{plot_groups}{A named list describing which gate or gates were plotted in each subplot.}
#'     \item{plot_titles}{A named list containing the parent population title for each subplot.}
#'     \item{df_tree}{The hierarchy table returned by `get_hierarchy_all_pops()` and used to build the figure.}
#'   }
#'
#' @details
#' The function uses `get_hierarchy_all_pops()` to obtain a hierarchy table. The
#' table is walked row by row so that plots follow the gating hierarchy order.
#' If a child population has an empty `Dimensions` value, it is plotted alone.
#' If multiple child populations share the same `Dimensions` value, they are
#' grouped into the same subplot and shown together.
#'
#' The plotting itself is delegated to `magicPlot_gs_gates()`. This wrapper
#' determines which gates should appear in each subplot, then calls
#' `magicPlot_gs_gates()` repeatedly and combines the resulting plots.
#'
#' @keywords flowMagic plotting GatingSet hierarchy
#' @export
#'
#' @examples
#' \donttest{
#' # Plot the full hierarchy for one sample
#' p <- magicPlot_gs_hierarchy(
#'   gs = gs,
#'   sample_id = "sample_01.fcs",
#'   n_col_wrap = 3,
#'   type = "dens"
#' )
#'
#' # Plot and export the hierarchy figure as a TIFF file
#' p <- magicPlot_gs_hierarchy(
#'   gs = gs,
#'   sample_id = "sample_01.fcs",
#'   n_col_wrap = 3,
#'   type = "dens",
#'   path_output = "~/main/Results/sample_01_hierarchy.tiff",
#'   plot_width = 4,
#'   plot_height = 4,
#'   dpi = 300
#' )
#'
#' # Manually control text sizes instead of using automatic sizing
#' p <- magicPlot_gs_hierarchy(
#'   gs = gs,
#'   sample_id = "sample_01.fcs",
#'   n_col_wrap = 3,
#'   auto_size = FALSE,
#'   size_title = 10,
#'   size_pol_name = 4,
#'   size_axis_text = 8,
#'   size_title_x = 10,
#'   size_title_y = 10,
#'   type = "dens"
#' )
#' }

magicPlot_gs_hierarchy <- function(gs,
                                   sample_id,
                                   n_col_wrap = 4,
                                   size_points = 0.5,
                                   size_title = NULL,
                                   concavity_val = 50,
                                   add_labels = TRUE,
                                   size_pol_name = NULL,
                                   auto_size = TRUE,
                                   return_plot_list = FALSE,
                                   path_output = NULL,
                                   plot_width = 4,
                                   plot_height = 4,
                                   dpi = 300,
                                   ...) {
  
  # =========================================================================
  # 1. Start workflow and validate the selected sample
  # =========================================================================
  # This function creates a hierarchy-style summary figure from a GatingSet.
  #
  # For one selected sample, it:
  #   1. reads the gating hierarchy,
  #   2. identifies which gates should be plotted alone,
  #   3. identifies which sibling gates share the same bivariate dimensions
  #      and should therefore be plotted together,
  #   4. generates one plot per gating step,
  #   5. combines all plots into one wrapped figure,
  #   6. optionally exports the wrapped figure to disk.
  #
  # In each subplot:
  #   - the title is the parent population,
  #   - the polygon labels are the child gate names.
  #
  # Example:
  #   title: "Parent pop: CD45+"
  #   labels inside plot: "NK_cells", "B_cells", "T_cells"
  
  message("$$$ Plot GatingSet hierarchy $$$")
  message(sprintf("Sample: %s", sample_id))
  
  # Check that the selected sample exists in the GatingSet.
  # This is important because later we use gs[[sample_id]].
  if (!(sample_id %in% sampleNames(gs))) {
    stop("sample_id is not present in sampleNames(gs): ", sample_id)
  }
  
  
  # =========================================================================
  # 2. Extract hierarchy information from the selected GatingHierarchy
  # =========================================================================
  # get_hierarchy_all_pops() is expected to return a list containing df_tree.
  #
  # df_tree should contain at least these columns:
  #   Mother     : parent population
  #   Children   : child gate/population
  #   Dimensions : empty if the child gate is plotted alone,
  #                or a shared ID such as "same_dims_1" when multiple sibling
  #                gates should be plotted together.
  #
  # Example rows:
  #
  #   Mother   Children   Dimensions
  #   root     Live       ""
  #   CD45+    NK_cells   same_dims_1
  #   CD45+    B_cells    same_dims_1
  #   CD45+    T_cells    same_dims_1
  #
  # In this example, Live is plotted alone, while NK_cells, B_cells, and
  # T_cells are plotted together because they share the same Dimensions ID.
  
  message("Get hierarchy information")
  
  hierarchy_info <- flowMagic::get_hierarchy_all_pops(gs[[sample_id]])
  df_tree <- hierarchy_info$df_tree
  
  # The root node is not itself a gate drawn inside another population.
  # If it appears as a child, remove it from the plotting table.
  df_tree <- df_tree[df_tree$Children != "root", ]
  
  # Stop if no actual gated populations are available.
  if (nrow(df_tree) == 0) {
    stop("No non-root populations found in the GatingSet hierarchy.")
  }
  
  
  # =========================================================================
  # 3. Build plotting groups by walking through df_tree in hierarchy order
  # =========================================================================
  # We now convert df_tree into two parallel lists:
  #
  #   plot_groups:
  #     tells the function which gate or gates should appear in each subplot.
  #
  #   plot_titles:
  #     tells the function what title each subplot should have.
  #
  # The names of plot_groups are internal plot IDs. They are mostly useful for
  # storing and returning the plots.
  #
  # Example:
  #
  #   plot_groups[["Live_cells"]] <- "Live_cells"
  #   plot_titles[["Live_cells"]] <- "root"
  #
  #   plot_groups[["NK_cells_B_cells_T_cells"]] <- c("NK_cells", "B_cells", "T_cells")
  #   plot_titles[["NK_cells_B_cells_T_cells"]] <- "CD45+"
  #
  # The function iterates row-by-row through df_tree to preserve the biological
  # hierarchy order: root-to-leaf.
  #
  # used_dims keeps track of same-dimension groups that were already added.
  # Without used_dims, grouped gates would be plotted repeatedly, once per row.
  
  plot_groups <- list()
  plot_titles <- list()
  used_dims <- character(0)
  
  message("Build plot groups following df_tree order")
  
  for (i in seq_len(nrow(df_tree))) {
    
    # Extract parent, child, and dimension-group information for this row.
    mother_name <- df_tree$Mother[i]
    child_name <- df_tree$Children[i]
    dim_id <- df_tree$Dimensions[i]
    
    
    # -----------------------------------------------------------------------
    # 3A. Single-gate case
    # -----------------------------------------------------------------------
    # If dim_id is empty or NA, this gate does not belong to a group of sibling
    # gates drawn on the same bivariate dimensions. Therefore, it gets its own
    # subplot.
    
    if (is.na(dim_id) || dim_id == "") {
      
      # Store the single gate as its own plotting group.
      plot_groups[[child_name]] <- child_name
      
      # The plot title should be the parent population, not the child gate.
      plot_titles[[child_name]] <- mother_name
      
      message(sprintf(
        "Added single plot: mother = %s, gate = %s",
        mother_name,
        child_name
      ))
      
      
      # -----------------------------------------------------------------------
      # 3B. Grouped-gate case
      # -----------------------------------------------------------------------
      # If dim_id is not empty, this child belongs to a group of sibling gates
      # that share the same bivariate dimensions.
      #
      # Example:
      #   CD45+ -> NK_cells, B_cells, T_cells
      #
      # These gates should be overlaid together in the same subplot.
      
    } else {
      
      # Add this same-dimension group only once.
      # The first time this dim_id appears, we collect all gates with this dim_id.
      # Later rows with the same dim_id are skipped.
      if (!(dim_id %in% used_dims)) {
        
        # Select all gates belonging to this same-dimension group.
        df_dim <- df_tree[df_tree$Dimensions == dim_id, ]
        
        # These are the child gates that will be plotted together.
        gate_names_current <- df_dim$Children
        
        # All rows in the same-dimension group should share the same parent.
        # Use the first parent name as the subplot title.
        mother_name_current <- df_dim$Mother[1]
        
        # Build an internal name for this grouped subplot.
        # This name is used as the list key in plot_groups and plot_list.
        group_name <- paste(gate_names_current, collapse = "_")
        
        # Store the grouped gates and their parent title.
        plot_groups[[group_name]] <- gate_names_current
        plot_titles[[group_name]] <- mother_name_current
        
        # Mark this dimension group as used.
        used_dims <- c(used_dims, dim_id)
        
        message(sprintf(
          "Added grouped plot: mother = %s, gates = %s",
          mother_name_current,
          paste(gate_names_current, collapse = ", ")
        ))
      }
    }
  }
  
  
  # =========================================================================
  # 4. Compute wrapped-grid dimensions
  # =========================================================================
  # The number of plots determines the number of rows needed in the final figure.
  #
  # n_col_wrap controls the maximum number of columns.
  # If there are fewer plots than n_col_wrap, the number of columns is reduced.
  #
  # These values are used for:
  #   - patchwork layout,
  #   - automatic text-size scaling,
  #   - exported figure size.
  
  n_plots <- length(plot_groups)
  n_cols <- min(n_col_wrap, n_plots)
  n_rows <- ceiling(n_plots / n_cols)
  
  message(sprintf("Number of hierarchy plots: %s", n_plots))
  message(sprintf("Grid: %s columns x %s rows", n_cols, n_rows))
  
  
  # =========================================================================
  # 5. Collect optional user arguments from ...
  # =========================================================================
  # The ... argument allows users to pass additional plotting options without
  # adding every possible magicPlot() option to this wrapper.
  #
  # For example, the user may call:
  #
  #   magicPlot_gs_hierarchy(
  #     gs = gs,
  #     sample_id = "sample_01.fcs",
  #     type = "dens",
  #     size_axis_text = 10,
  #     size_title_x = 12,
  #     size_title_y = 12,
  #     show_legend = FALSE
  #   )
  #
  # In this function, those extra arguments are captured by ...
  # This line converts them into a normal named list:
  #
  #   args_list <- list(...)
  #
  # If the user supplied:
  #
  #   type = "dens", size_axis_text = 10
  #
  # then args_list becomes:
  #
  #   list(type = "dens", size_axis_text = 10)
  #
  # We use a list because it can be inspected, modified, and then passed to
  # another function later.
  
  args_list <- list(...)
  
  
  # =========================================================================
  # 6. Automatically scale plot text sizes, if requested
  # =========================================================================
  # When many plots are shown in one grid, text sizes need to be adjusted.
  #
  # auto_size = TRUE fills in reasonable defaults for:
  #   - axis tick labels,
  #   - x-axis title,
  #   - y-axis title,
  #   - polygon gate labels,
  #   - subplot title.
  #
  # Important:
  #   User-provided values are not overwritten.
  #
  # Example:
  #   If the user passes size_axis_text = 12 through ...,
  #   this function keeps that value and does not replace it.
  
  if (auto_size == TRUE) {
    
    message("Apply automatic plot-size scaling")
    
    # Use the larger grid dimension as a simple measure of layout complexity.
    # A 4 x 4 grid has scale_factor = 4.
    # A larger scale_factor leads to slightly smaller text.
    scale_factor <- max(n_cols, n_rows)
    
    # Axis tick label size.
    # Only set it automatically if the user did not provide size_axis_text.
    if (!("size_axis_text" %in% names(args_list))) {
      args_list$size_axis_text <- max(7, round(20 / sqrt(scale_factor), 1))
    }
    
    # X-axis title size.
    # Only set it automatically if the user did not provide size_title_x.
    if (!("size_title_x" %in% names(args_list))) {
      args_list$size_title_x <- max(9, round(24 / sqrt(scale_factor), 1))
    }
    
    # Y-axis title size.
    # Only set it automatically if the user did not provide size_title_y.
    if (!("size_title_y" %in% names(args_list))) {
      args_list$size_title_y <- max(9, round(24 / sqrt(scale_factor), 1))
    }
    
    # Gate label size printed inside the polygon.
    # This is a formal argument of magicPlot_gs_hierarchy(), so it is checked
    # directly rather than inside args_list.
    if (is.null(size_pol_name)) {
      size_pol_name <- max(4.5, round(8 / sqrt(scale_factor), 1))
    }
    
    # Subplot title size, for example "Parent pop: CD45+".
    if (is.null(size_title)) {
      size_title <- max(10, round(17 / sqrt(scale_factor), 1))
    }
    
    message(sprintf("Auto-size axis text: %s", args_list$size_axis_text))
    message(sprintf("Auto-size x title: %s", args_list$size_title_x))
    message(sprintf("Auto-size y title: %s", args_list$size_title_y))
    message(sprintf("Auto-size gate label: %s", size_pol_name))
    message(sprintf("Auto-size panel title: %s", size_title))
    
    
    # -------------------------------------------------------------------------
    # If auto_size is disabled, use simple defaults for the values that were not
    # supplied directly.
    # -------------------------------------------------------------------------
    
  } else {
    
    if (is.null(size_pol_name)) {
      size_pol_name <- 4
    }
    
    if (is.null(size_title)) {
      size_title <- 10
    }
  }
  
  
  # =========================================================================
  # 7. Generate one ggplot object per hierarchy step
  # =========================================================================
  # This section loops through plot_groups and creates each subplot.
  #
  # Each element of plot_groups is passed to magicPlot_gs_gates().
  #
  # Example 1: single gate
  #   current_gate_names = "Live_cells"
  #
  # Example 2: grouped gates
  #   current_gate_names = c("NK_cells", "B_cells", "T_cells")
  #
  # magicPlot_gs_gates() handles the extraction and plotting of these gates.
  
  plot_list <- list()
  
  for (i in seq_along(plot_groups)) {
    
    current_plot_name <- names(plot_groups)[i]
    current_gate_names <- plot_groups[[i]]
    current_title <- plot_titles[[current_plot_name]]
    
    message(sprintf(
      "---- Plotting mother %s with gates: %s ----",
      current_title,
      paste(current_gate_names, collapse = ", ")
    ))
    
    
    # -----------------------------------------------------------------------
    # 7A. Build the argument list for magicPlot_gs_gates()
    # -----------------------------------------------------------------------
    # We need to call magicPlot_gs_gates() many times, once per subplot.
    #
    # Some arguments are controlled by this wrapper:
    #   - gs
    #   - sample_id
    #   - current_gate_names
    #   - size_points
    #   - concavity_val
    #   - add_labels
    #   - size_pol_name
    #
    # Other arguments come from the user through ... and are stored in args_list.
    #
    # Example:
    #   args_list might contain:
    #     list(type = "dens", size_axis_text = 10)
    #
    # The c() function combines both lists into one complete argument list.
    #
    # After this block, args_gates may look like:
    #
    #   list(
    #     gs = gs,
    #     sample_id = "sample_01.fcs",
    #     gate_names = c("NK_cells", "B_cells", "T_cells"),
    #     size_points = 0.5,
    #     concavity_val = 50,
    #     add_labels = TRUE,
    #     size_pol_name = 4,
    #     type = "dens",
    #     size_axis_text = 10
    #   )
    
    args_gates <- c(
      list(
        gs = gs,
        sample_id = sample_id,
        gate_names = current_gate_names,
        size_points = size_points,
        concavity_val = concavity_val,
        add_labels = add_labels,
        size_pol_name = size_pol_name
      ),
      args_list
    )
    
    
    # -----------------------------------------------------------------------
    # 7B. Call magicPlot_gs_gates() using do.call()
    # -----------------------------------------------------------------------
    # Because the function arguments are stored inside the list args_gates,
    # we use do.call().
    #
    # do.call(function_name, argument_list) means:
    #   "Call this function using the named values stored in this list."
    #
    # In other words:
    #
    #   do.call(flowMagic::magicPlot_gs_gates, args_gates)
    #
    # is equivalent to manually writing:
    #
    #   flowMagic::magicPlot_gs_gates(
    #     gs = gs,
    #     sample_id = sample_id,
    #     gate_names = current_gate_names,
    #     size_points = size_points,
    #     concavity_val = concavity_val,
    #     add_labels = add_labels,
    #     size_pol_name = size_pol_name,
    #     type = "dens",
    #     size_axis_text = 10
    #   )
    #
    # The advantage is that args_list can be modified earlier by auto_size
    # before we pass it onward.
    
    # args_gates now contains the complete set of arguments that will be passed
    # to magicPlot_gs_gates().
    #
    # It combines:
    #   1. arguments controlled directly by magicPlot_gs_hierarchy(), such as
    #      gs, sample_id, gate_names, size_points, concavity_val, add_labels,
    #      and size_pol_name;
    #   2. optional plotting arguments passed by the user through ..., stored in
    #      args_list, such as type, size_axis_text, size_title_x, size_title_y,
    #      show_legend, x/y limits, etc.
    #
    # Not every argument of magicPlot_gs_hierarchy() is included in args_gates.
    # For example, plot_width, plot_height, dpi, and path_output are used only for
    # exporting the final wrapped figure, so they are not passed to
    # magicPlot_gs_gates().
    
    p <- do.call(flowMagic::magicPlot_gs_gates, args_gates)
    
    
    # -----------------------------------------------------------------------
    # 7C. Add subplot title
    # -----------------------------------------------------------------------
    # The title should describe the parent population.
    # The gate names themselves are printed inside the polygons when
    # add_labels = TRUE.
    
    p <- p +
      ggplot2::ggtitle(sprintf("Parent pop: %s", current_title)) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(size = size_title)
      )
    
    plot_list[[current_plot_name]] <- p
  }
  
  
  # =========================================================================
  # 8. Combine all subplots into one wrapped figure
  # =========================================================================
  # patchwork::wrap_plots() arranges the list of ggplot objects into a grid.
  # The order of plot_list follows plot_groups, which was built in df_tree order.
  # Therefore, the final plot follows the gating hierarchy order.
  
  out <- patchwork::wrap_plots(
    plot_list,
    ncol = n_col_wrap
  )
  
  
  # =========================================================================
  # 9. Optionally export the wrapped figure to disk
  # =========================================================================
  # path_output is treated as the full output file path.
  #
  # Example:
  #   path_output = "~/main/Results/hierarchy_plot.tiff"
  #
  # ggsave() infers the export format from the file extension.
  #
  # The exported image size is computed from:
  #   plot_width  * number of columns
  #   plot_height * number of rows
  #
  # Example:
  #   If n_cols = 3, n_rows = 4, plot_width = 4, plot_height = 4:
  #   final export size = 12 x 16 inches.
  
  if (!is.null(path_output)) {
    
    message("$$$ Export hierarchy plot to disk $$$")
    
    export_width <- plot_width * n_cols
    export_height <- plot_height * n_rows
    
    message(sprintf("Export size: %.2f x %.2f inches", export_width, export_height))
    message(sprintf("Export file: %s", path_output))
    
    ggplot2::ggsave(
      filename = path_output,
      plot = out,
      width = export_width,
      height = export_height,
      dpi = dpi,
      units = "in"
    )
    
    message("Export completed.")
  }
  
  
  # =========================================================================
  # 10. Return output
  # =========================================================================
  # By default, return only the wrapped plot.
  #
  # If return_plot_list = TRUE, return a list with useful intermediate objects:
  #
  #   plot        : final wrapped patchwork plot
  #   plot_list   : individual ggplot objects before wrapping
  #   plot_groups : gate names used in each subplot
  #   plot_titles : parent population title for each subplot
  #   df_tree     : hierarchy table used to build the plot
  #
  # These intermediate objects are useful for debugging or for developing a
  # more advanced tree-layout version later.
  
  if (return_plot_list == TRUE) {
    return(
      list(
        plot = out,
        plot_list = plot_list,
        plot_groups = plot_groups,
        plot_titles = plot_titles,
        df_tree = df_tree
      )
    )
  }
  
  return(out)
}


#' get_hierarchy_all_pops
#'
#' function to plot hierarchy of all nodes/cell populations from a GatingHierarchy object.
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

  if (("package:flowWorkspace" %in% search()) == F) {
    warning("flowWorkspace library not loaded, attempt to loading...")
    library(flowWorkspace)
  }

  if (("package:flowCore" %in% search()) == F) {
    warning("flowCore library not loaded, attempt to loading...")
    library(flowCore)
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
  # ------------------------------------------------------------
  # Add channel information as node hover text
  # ------------------------------------------------------------
  
  # Get one flowFrame from the gating hierarchy / gating set.
  # This is only used to access the parameter metadata:
  # ff@parameters@data
  ff <- flowMagic::get_flowframe_from_gs(
    gs = gh,
    node_name = "root",
    sample_id = 1
  )
  
  # Extract metadata for the current flowFrame.
  # This contains columns such as:
  # name, desc, range, minRange, maxRange
  df_metadata <- ff@parameters@data
  
  # Add a title column to df_tree_nodes.
  #
  # In visNetwork, a column named `title` automatically becomes
  # the hover tooltip for each node.
  #
  # For each population:
  # 1. get the gate channels used to define that population
  # 2. convert those channels into readable labels using df_metadata
  # 3. store the result as HTML tooltip text

  df_tree_nodes$title <- sapply(
    df_tree_nodes$label,
    function(pop) {
      
      channels <- get_gate_channels(
        gs = gh,
        pop = pop
      )
      
      channel_info <- format_channel_info(
        channels = channels,
        df_metadata = df_metadata
      )
      
      paste0(
        "<b>Population:</b> ",
        pop,
        "<br><br><b>Channels used for this gate:</b><br>",
        gsub("\n", "<br>", channel_info)
      )
    }
  )
  
  # ------------------------------------------------------------
  # Build visNetwork
  # ------------------------------------------------------------
  
  visnet <- visNetwork(
    edges = df_tree_edges,
    nodes = df_tree_nodes
  )
  
  visnet <- visnet %>%
    visEdges(arrows = "to") %>%
    visHierarchicalLayout() %>%
    visNodes(borderWidth = 2)
  
  visnet <- visnet %>%
    visOptions(
      highlightNearest = TRUE,
      nodesIdSelection = TRUE
    )
  
  # ------------------------------------------------------------
  # Copy the node label when clicking a node
  # ------------------------------------------------------------

  visnet <- visnet %>%
    visEvents(
      selectNode = "
        function(params) {

          // This function runs whenever the user selects/clicks a node
          // in the visNetwork plot.

          // params contains information about what was selected.
          // params.nodes is a list of selected node IDs.
          // If the user clicked somewhere empty, params.nodes may be empty.

          if (!params.nodes || params.nodes.length === 0) {
            // If there are no selected nodes, stop here.
            // This prevents errors when the user clicks empty space.
            return;
          }

          // Get the ID of the first selected node.
          // Usually only one node is selected, so we use params.nodes[0].
          var nodeId = params.nodes[0];

          // Use the node ID to get the full node information
          // from the visNetwork internal node dataset.
          //
          // this = the visNetwork object
          // this.body.data.nodes = the table of nodes used by the plot
          // get(nodeId) = retrieve the node with this ID
          var node = this.body.data.nodes.get(nodeId);

          // Continue only if:
          // 1. the node was found
          // 2. the node has a label
          //
          // This protects against errors if something unexpected happens.
          if (node && node.label) {

            // Create a temporary invisible text box.
            // JavaScript can copy text from a selected text box.
            var textarea = document.createElement('textarea');

            // Put the node label into the temporary text box.
            // This is the text we want to copy.
            textarea.value = node.label;

            // Add the temporary text box to the web page.
            // It has to be part of the page before we can select/copy from it.
            document.body.appendChild(textarea);

            // Select the text inside the temporary text box.
            // This is like highlighting the text manually with the mouse.
            textarea.select();

            try {
              // Copy the selected text to the clipboard.
              // This is the actual copy step.
              document.execCommand('copy');

            } catch (err) {
              // If automatic copying fails, show a small backup box.
              // The user can manually copy the node label from this box.
              window.prompt('Copy this text:', node.label);
            }

            // Remove the temporary text box from the page.
            // The user never needs to see it.
            document.body.removeChild(textarea);
          }
        }
        "
    )
  
  #-------------------- export visnetwork -------------------------
  if(export_visnet==T){
    dir.create(paste0(path.output,"/visnet_plot/"),recursive = F)
  }
  end<-Sys.time()
  time_taken<-end-start
  print("Time of execution:")
  print(time_taken)
  print("Done")
  return(list(df_tree=df_tree,df_tree_levels=df_tree_levels,visnet=visnet,hierarchical_tree=hierarchical_tree))
}

