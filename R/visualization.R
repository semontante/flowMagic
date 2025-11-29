#' magicPlot
#'
#' function to generate the scatter plot with colored density of the events.
#' @param df Dataframe of bivariate markers expression (with labels if gates to plot).
#' @param type Type of plot to generated. "dens"= plot with bivariate density. "ML"= plot with color layer based on machine learning class assignments 
#' @param polygons_coords_list list of gates coordinates. Needed if labels not included in df. Default to NULL.
#' @param show_legend Show legend if type="ML". Default to True.
#' @param size_axis_text Size of axis ticks labels. Default to 18.
#' @param size_title_x Size of x axis label title.
#' @param size_title_y Size of y axis label title.
#' @param treat_0_as_gate Treat 0 label as gate. Defaul to False (0 label is background)
#' @param x_lab Label of x axis. Default to NULL.
#' @param y_lab Label of y axis. Default to NULL.
#' @param gates_to_plot Select labels to plot.
#' @param apply_manual_scale Apply predifined scale of colors. Default to True.
#' @param size_points Size of points in scatter plot.
#' @param concavity_val Concavity value. Default to 5. Higher value, less jagged boundaries
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
    if("df_test_original" %in% colnames(df)){
      df<-df$df_test_original
    }
    p <- magicPlot(df,...) + ggtitle(all_names[i]) + theme(plot.title = element_text(size = size_title))
    return(p)
  })
  
  # Combine into a grid layout, e.g., 3 columns
  combined_plot <- wrap_plots(plot_list, ncol = n_col_wrap)
  
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

