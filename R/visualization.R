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
#' @param aspect_ratio Aspect ratio value. If = 1, y and x axis have equal distance between the ticks labels. Default to NULL.
#' @param x_lim1 Minimum limit x axis. Default to NULL.
#' @param x_lim2 Max limit x axis. Default to NULL.
#' @param y_lim1 Minimum limit y axis. Default to NULL.
#' @param y_lim2 Max limit y axis. Default to NULL.
#' @return ggplot.
#' @export
#' @examples 
#' \donttest{magicPlot()}

magicPlot<-function(df,type="dens",polygons_coords_list=NULL,
                    show_legend=T,size_axis_text=18,size_title_x=20,size_title_y=20,
                    treat_0_as_gate=F,x_lab="x",y_lab="y",gates_to_plot=NULL,apply_manual_scale=T,
                    hull_only=F,size_points=1,concavity_val=5,aspect_ratio=NULL,x_lim1=NULL,x_lim2=NULL,
                    y_lim1=NULL,y_lim2=NULL){
  
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
  
  if(is.null(aspect_ratio)==F){
    magicggplot<- magicggplot + coord_fixed(ratio = aspect_ratio) + theme(plot.margin=unit(c(0,0,0,0), "pt"))
  }
  if(is.null(x_lim1)==F){
    magicggplot<- magicggplot+ xlim(x_lim1, x_lim2) + ylim(y_lim1,y_lim2)
  }
  # The ratio represents the number of units on the y-axis equivalent to one unit on the x-axis. 
  # The default, ratio = 1, ensures that one unit on the x-axis is the same length as one unit on the y-axis.
  return(magicggplot)
}




