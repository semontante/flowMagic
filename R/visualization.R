#' magicPlot
#'
#' function to generate the scatter plot with colored density of the events.
#' @param df Dataframe of bivariate markers expression (with labels if gates to plot).
#' @param type Type of plot to generated. "dens"=bivariate density plot. "ML"=events assignments plot
#' @param polygons_coords_list list of gates coordinates. Needed if labels not included in df. Default to NULL.
#' @param show_legend Show legend if type="ML". Default to True.
#' @param size_axis_text Size of axis ticks labels. Default to 18.
#' @param size_title_x Size of x axis label title.
#' @param size_title_y Size of y axis label title.
#' @param treat_0_as_gate Treat 0 label as gate. Defaul to False (0 label is background)
#' @param x_lab Label of x axis.
#' @param y_lab Label off y axis.
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
#' @return ggplot.
#' @export
#' @examples 
#' \donttest{magicPlot()}

magicPlot<-function(df, type = "dens", polygons_coords_list = NULL, show_legend = T, 
          size_axis_text = 18, size_title_x = 20, size_title_y = 20, 
          treat_0_as_gate = F, x_lab = "x", y_lab = "y", gates_to_plot = NULL, 
          apply_manual_scale = T, hull_only = F, size_points = 1, concavity_val = 20, 
          aspect_ratio = NULL, x_lim1 = NULL, x_lim2 = NULL, y_lim1 = NULL, 
          y_lim2 = NULL,add_labels=F,map_label_polygon=NULL,size_pol_name=6,show_marginals=F){

  if (type == "no_gate" || ncol(df) == 2) {
    colPalette <- colorRampPalette(c("blue", "turquoise", "green", "yellow", "orange", "red"))
    col <- densCols(df[, c(1, 2)], colramp = colPalette, nbin = 200)
    
    df_plot <- data.frame(x = df[[1]], y = df[[2]], col = col)
    
    p <- ggplot(df_plot, aes(x = x, y = y)) + geom_point(color = df_plot$col, size = size_points, shape = 16) 
    p <- p + xlab(x_lab) + ylab(y_lab)
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
  if (type == "dens") {
    colPalette <- colorRampPalette(c("blue", "turquoise", 
                                     "green", "yellow", "orange", "red"))
    col <- densCols(df[, c(1, 2)], colramp = colPalette, 
                    nbin = 200)
    df <- cbind(df, col)
    df$col <- as.character(df$col)
    if (is.null(polygons_coords_list) == T) {
      list_df_hull <- get_hull_all_gates(df, concavity_val = concavity_val)
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
    magicggplot <- magicggplot + xlab(x_lab) + ylab(y_lab)
  }
  else if (type == "ML") {
    df$classes <- as.factor(df$classes)
    magicggplot <- ggplot(data = df, aes(x = x, y = y)) + 
      geom_point(aes(color = classes), size = size_points, 
                 shape = 16)
    magicggplot <- magicggplot + theme(axis.text = element_text(size = size_axis_text), 
                                       axis.title.x = element_text(size = size_title_x, 
                                                                   face = "bold"), axis.title.y = element_text(size = size_title_y, 
                                                                                                               face = "bold"))
    magicggplot <- magicggplot + theme(panel.grid.major = element_blank(), 
                                       panel.grid.minor = element_blank(), panel.background = element_blank(), 
                                       axis.line = element_line(colour = "black"))
    if (apply_manual_scale == T) {
      magicggplot <- magicggplot + scale_color_manual(values = c(`0` = "black", 
                                                                 `1` = "red", `2` = "blue", `3` = "green", `4` = "yellow", 
                                                                 `5` = "purple", `6` = "orange", `7` = "pink", 
                                                                 `8` = "brown", `9` = "gray"))
    }
    magicggplot <- magicggplot + theme(legend.key.size = unit(1, 
                                                              "cm"), legend.title = element_text(size = 20), legend.text = element_text(size = 20))
    magicggplot <- magicggplot + guides(color = guide_legend(override.aes = list(size = 10))) + 
      labs(color = "Assignment")
    if (show_legend == F) {
      magicggplot <- magicggplot + theme(legend.position = "none")
    }
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
  if(add_labels==T){
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
#' @export
#' @examples 
#' \donttest{magic_plot_wrap()}

magic_plot_wrap<-function(list_gated_data,n_col_wrap=3,size_title=10,...){
  all_names<-names(list_gated_data)
  plot_list <- lapply(seq_along(list_gated_data),function(i){
    df <- list_gated_data[[i]]
    p <- magicPlot(df, type = "dens",...) + ggtitle(all_names[i]) + theme(plot.title = element_text(size = size_title))
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
#' @export
#' @examples 
#' \donttest{magicplot_3D()}

magicplot_3D<-function(df,class_col=F,x_lab="x",y_lab="y",z_lab="z",type="ML",size_p=1){
  if(class_col==F){
    if(ncol(df)!=3){
      stop("the input df is expected to have 3 columns if class_col==F")
    }
    plotly_plot<-plot_ly(x = df[,1], y = df[,2], 
                         z = df[,3], type = "scatter3d", mode = "markers",
                         marker = list(size = size_p))
    
    plotly_plot<- plotly_plot %>% layout(scene = list(
      xaxis = list(title = x_lab),
      yaxis = list(title = y_lab),
      zaxis = list(title = z_lab)
    ))
  }else if(class_col==T){
    if(ncol(df)!=4){
      stop("the input df is expected to have 4 columns if class_col==T")
    }
    if(type=="ML"){
      df[,4]<-as.factor(df[,4])
      
      plotly_plot<-plot_ly(x = df[,1], y = df[,2], 
                           z = df[,3], type = "scatter3d", mode = "markers",marker = list(size = size_p),
                           color =df[,4],colors = c("blue", "red"))
      
      plotly_plot<- plotly_plot %>% layout(scene = list(
        xaxis = list(title = x_lab),
        yaxis = list(title = y_lab),
        zaxis = list(title = z_lab)
      ))
      plotly_plot<- plotly_plot %>% layout(showlegend = FALSE)
      
    }else if(type=="mesh"){
      df_class<-df[df[,4]==1,]
      
      coords<-cbind(df_class[,1], df_class[,2],
                    df_class[,3])
      
      hull <- convhulln(coords)
      
      red_rgb <- rep("rgb(255,0,0)", nrow(coords))
      
      plotly_plot<-plot_ly() %>%
        # Full data as scatter
        add_trace(
          type = "scatter3d",
          mode = "markers",
          x = df[,1],
          y = df[,2],
          z = df[,3],
          color = df[,4],
          marker = list(size = size_p),
          colors = c("blue", "black"),
          name = "Parent"
        ) %>%
        # Add convex hull as simulated 3D gate
        add_trace(
          type = "mesh3d",
          x = coords[, 1],
          y = coords[, 2],
          z = coords[, 3],
          i = hull[, 1] - 1,
          j = hull[, 2] - 1,
          k = hull[, 3] - 1,
          vertexcolor = red_rgb,
          opacity = 0.3,
          color = I("red"),
          flatshading = TRUE,
          lighting = list(ambient = 1, diffuse = 0, specular = 0, roughness = 0, fresnel = 0),
          lightposition = list(x = 0, y = 0, z = 100),
          name = "Gated pop"
        )
      
      plotly_plot<- plotly_plot %>% layout(scene = list(
        xaxis = list(title = x_lab),
        yaxis = list(title = y_lab),
        zaxis = list(title = z_lab)
      ))
      
      plotly_plot<- plotly_plot %>% layout(showlegend = FALSE)
      
    }

  }
  return(plotly_plot)
 
}


