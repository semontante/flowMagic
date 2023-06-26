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
#' @param concavity_val Concavity value. Default to 5.
#' @return ggplot.
#' @export
#' @examples 
#' \donttest{magicPlot()}

magicPlot<-function(df,type="dens",polygons_coords_list=NULL,
                    show_legend=T,size_axis_text=18,size_title_x=20,size_title_y=20,
                    treat_0_as_gate=F,x_lab="x",y_lab="y",gates_to_plot=NULL,apply_manual_scale=T,
                    size_points=1,concavity_val=5){
  
  if(type=="no_gate" || ncol(df)==2){
    colPalette <- colorRampPalette(c("blue", "turquoise", "green", "yellow", "orange", "red"))
    col <- densCols(df[,c(1,2)], colramp = colPalette,nbin = 200) # get colors based on bivariate density
    plot<-graphics::plot(df[,c(1,2)],col=col,pch=".",cex=size_points)
    return(plot)
  }
  colnames(df)<-c("x","y","classes")
  df$classes<-as.character(df$classes)
  if(treat_0_as_gate==T){
    # we treat 0 as a gate
    df$classes<-as.numeric(df$classes)
    inds<-which(df$classes==0)
    df$classes[inds]<-max(unique(df$classes))+1
    df$classes<-as.character(df$classes)
  }
  if(is.null(gates_to_plot)==F){
    check_classes<-df$classes %in% gates_to_plot
    inds<-which(check_classes==T)
    df$classes[-inds]<-"0"
  }
  ##################################### generate plot based on type info ############################
  if(type=="dens"){
    colPalette <- colorRampPalette(c("blue", "turquoise", "green", "yellow", "orange", "red"))
    col <- densCols(df[,c(1,2)], colramp = colPalette,nbin = 200) # get colors based on bivariate density
    df<-cbind(df,col)
    df$col<-as.character(df$col)
    if(is.null(polygons_coords_list)==T){
      list_df_hull<-get_hull_all_gates(df,concavity_val = concavity_val)
      # we don't consider the 0 class (no gate) 
      vec<-names(list_df_hull)
      inds<-which(vec=="0")
      if(length(inds)!=0){
        list_df_hull<-list_df_hull[-inds]
      }
      df_hull<-do.call(rbind,list_df_hull)
    }else{
      df_hull<-do.call(rbind,polygons_coords_list)
    }
    magicggplot<-ggplot(data=df,aes(x=x,y=y,col=factor(1:nrow(df)))) + geom_point(size=size_points,shape=16)
    magicggplot<- magicggplot + scale_colour_manual(values=col)
    magicggplot<- magicggplot + theme(legend.position="none")
    magicggplot<- magicggplot + theme(axis.text=element_text(size=size_axis_text),axis.title.x = element_text(size = size_title_x,face="bold"),
                                      axis.title.y = element_text(size=size_title_y,face="bold"))
    magicggplot <- magicggplot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                       panel.background = element_blank(), axis.line = element_line(colour = "black"))
    magicggplot<- magicggplot + geom_polygon(data = df_hull,aes(x=x,y=y,group=group_gate),fill=NA,color="black",size=1)
    
    magicggplot<-magicggplot + xlab(x_lab) + ylab(y_lab)
    
  }else if(type=="ML"){
    df$classes<-as.factor(df$classes)
    magicggplot<-ggplot(data=df,aes(x=x,y=y)) + geom_point(aes(color=classes),size=size_points,shape=16)
    magicggplot<- magicggplot + theme(axis.text=element_text(size=size_axis_text),axis.title.x = element_text(size = size_title_x,face="bold"),
                                      axis.title.y = element_text(size=size_title_y,face="bold"))
    magicggplot <- magicggplot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                       panel.background = element_blank(), axis.line = element_line(colour = "black"))
    if(apply_manual_scale==T){
      magicggplot<- magicggplot + scale_color_manual(values=c("0"="black","1"="red","2"="blue","3"="green",
                                                              "4"="yellow","5"="purple","6"="orange","7"="pink",
                                                              "8"="brown","9"="gray"))
    }

    magicggplot<-magicggplot + theme(legend.key.size = unit(1, 'cm'),legend.title = element_text(size=20),legend.text = element_text(size=20))
    magicggplot<-magicggplot + guides(color = guide_legend(override.aes = list(size = 10))) + labs(color="Assignment")
    if(show_legend==F){
      magicggplot<- magicggplot + theme(legend.position="none")
    }
    magicggplot<-magicggplot + xlab(x_lab) + ylab(y_lab)
  }
  return(magicggplot)
}

#' .densRange
#'
#' function to get the correct coordinates (correct range) based on the density coordinates x and y generated by the density() function.
#' @param x Coordinates of x axis.
#' @param y Coordinates of y axis
#' @param gate gate threshold.
#' @param pos refere to the pos in deGatePlot() function.
#' @return list of numbers.
#' @export
#' @examples 
#' \donttest{.densRange()}


.densRange <- function(x, y, gate, pos = FALSE){
  
  ##==================================================================
  ## Plots the output of density-estimate gating method
  ##  Args:
  ##   x: 'x' slot of the density returned by the 'density()' function
  ##   y: 'y' slot of the density returned by the 'density()' function
  ##   gate: a threshold given by the '.densityGating()' function
  ##   pos: refer to the 'pos' in 'deGatePlot()' function
  ## Value:
  ##   (x,y) coordinates
  ##------------------------------------------------------------------
  pts <- list()
  if(is.na(pos))
    return(list(x=c(x,tail(x,1),x[1]), y=c(y,min(y),min(y))))
  if(pos){
    x.pts <- c(x[which(x>=gate)], tail(x[which(x>=gate)],1), x[which(x>=gate)][1])
    y.pts <- c( y[which(x>=gate)], min(y[which(x>=gate)]), min(y[which(x>=gate)]))
  }else{
    x.pts <- c(x[which(x<gate)][1], x[which(x<gate)], gate)
    y.pts <- c(min(y[which(x>=gate)]),y[which(x<gate)], min(y[which(x<gate)]))
  }
  pts$x <- x.pts
  pts$y <- y.pts
  return(pts)
}


