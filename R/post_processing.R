
#' get_hull_all_gates
#' 
#' function to get the convex hull of all gates.
#' @param gated_df dataframe with labels (third column).
#' @param concavity_val Values of concavity. Default to 1.
#' @return List of dataframes.
#' @export
#' @examples 
#' \donttest{get_hull_all_gates()}

get_hull_all_gates<-function(gated_df,concavity_val=1){
  colnames(gated_df)<-c("x","y","classes")
  all_classes<-unique(gated_df$classes)
  list_df_hull<-list()
  for(classes in all_classes){
    inds<-which(gated_df$classes==classes)
    df_current_classes<-gated_df[inds,]
    df_current_classes_hull_values<-as.data.frame(concaveman::concaveman(as.matrix(df_current_classes[,c(1,2)]),concavity=concavity_val))
    if(concavity_val <= 3){
      df_current_classes_hull_values<-smooth_hull(hull_df,spar=0.7)
    }
    vec_group<-rep(sprintf("%s",classes),nrow(df_current_classes_hull_values))
    df_current_hull<-cbind(df_current_classes_hull_values,vec_group)
    colnames(df_current_hull)<-c("x","y","group_gate")
    list_df_hull[[sprintf("%s",classes)]]<-df_current_hull
  }
  return(list_df_hull)
}

#' extract_polygon_gates
#' 
#' function to extract the polygon gates objects based on the convex hull and classes.
#' @param gated_df dataframe with labels (third column).
#' @param concavity_val Concavity of polygons. Default to 1.
#' @return List of dataframes.
#' @export
#' @examples 
#' \donttest{extract_polygon_gates()}

extract_polygon_gates<-function(gated_df,concavity_val=1){
  row.names(gated_df)<-NULL
  colnames(gated_df)<-c("x","y","classes") 
  gated_df$classes<-as.character(gated_df$classes)
  all_classes<-unique(gated_df$classes)
  ind<-which(all_classes=="0")
  if(length(ind)!=0){
    all_classes<-all_classes[-ind]
  }
  if(length(all_classes)==0){
    # only 0 class, no polygons
    return(NULL)
  }
  message(sprintf("all classes: %s",paste0(all_classes,collapse = ",")))
  ########################## find convex hull for each class ###################
  message(" ############# find convex hull for each class ############")
  list_df_hull<-get_hull_all_gates(gated_df,concavity_val=concavity_val)
  # we don't consider the 0 class (no gate) 
  vec<-names(list_df_hull)
  inds<-which(vec=="0")
  if(length(inds)!=0){
    list_df_hull<-list_df_hull[-inds]
  }
  return(list_df_hull)
}

#' extract_polygon_gates
#' 
#' function to check polygons intersection.
#' @param list_df_hull List of polygons coordinates
#' @return float
#' @export
#' @examples 
#' \donttest{check_polygons_intersection()}

check_polygons_intersection<-function(list_df_hull){
  ######################### convert convex hull in spatial polygon ################
  df_hull<-do.call(rbind,list_df_hull)
  polys <- lapply(unique(df_hull$group_gate), function(i) {
    sp::Polygons(list(sp::Polygon(df_hull[df_hull$group_gate==i, 1:2])), ID=i)
  })
  spa_polys <- sp::SpatialPolygons(polys) # spatial polygons based on a convex hull
  ####################### check final polygons  ntersections ##################
  n_polygons<-length(spa_polys)
  vec_check<-c()
  if(n_polygons==1){
    # There is only ony polygon no intersection
    max_area_intersect<-0
  }else{
    # check all combinations of polygons.
    for(i in 1:n_polygons){
      name_group_poly_i<-spa_polys@polygons[[i]]@ID # we keep track of the class under analysis
      message(sprintf("######## Analysis gate %s ###########",name_group_poly_i))
      poly_i<-spa_polys[i]
      all_new_poly_i_coords<-list() # all new possible coords of poly_i (poly_i intersects with more than one polygon)
      # for each polygon i we look for the situations in which there is an intersection with an other polygon j
      for(j in 1:n_polygons){
        name_group_poly_j<-spa_polys@polygons[[j]]@ID # we keep track of the class under analysis
        message(sprintf("-------- Analysis gate %s vs %s",name_group_poly_i,name_group_poly_j))
        if(name_group_poly_j!=name_group_poly_i){ # avoid comparison with itself
          poly_j<-spa_polys[j]
          poly_i_sf<-sf::st_as_sf(poly_i)
          poly_j_sf<-sf::st_as_sf(poly_j)
          area_intersect<-sf::st_intersection(sf::st_buffer(poly_i_sf, 0), sf::st_buffer(poly_j_sf, 0)) %>% sf::st_area()
          if(length(area_intersect)==0){
            area_intersect<-0
          }
          vec_check<-c(vec_check,area_intersect)
        }
      }
    }
    vec_check<-as.numeric(vec_check)
    max_area_intersect<-max(vec_check)
  }
  return(max_area_intersect)
}


#' compute_gates
#' 
#' function to assign events based on polygon gates.
#' @param gated_df dataframe with labels (third column).
#' @param list_final_polygons_coords List of dataframes containing polygon coordinates.
#' @param no_classes  Generate third column of labels. Default to False.
#' @return Dataframe.
#' @export
#' @examples 
#' \donttest{compute_gates()}
  
compute_gates<-function(gated_df,list_final_polygons_coords,no_classes=F){
  if(no_classes==F){
    row.names(gated_df)<-NULL # we have already saved the original root indices
    colnames(gated_df)<-c("x","y","classes")
    gated_df$classes<-as.character(gated_df$classes)
  }else if(no_classes==T){
    colnames(gated_df)<-c("x","y") 
    gated_df$classes<-rep(0,nrow(gated_df))
    gated_df$classes<-as.character(gated_df$classes)
  }
  n_polygons<-length(list_final_polygons_coords)
  all_classes_name<-names(list_final_polygons_coords)
  for(i in 1:n_polygons){
    class_name<-all_classes_name[i]
    coords_poly_current_class<-list_final_polygons_coords[[i]]
    vec_out<-sp::point.in.polygon(point.x=gated_df[,1], point.y=gated_df[,2], pol.x=coords_poly_current_class[,1],
                              pol.y=coords_poly_current_class[,2], mode.checked=FALSE)
    inds<-which(vec_out!=0)
    gated_df[inds,"classes"]<-class_name
  }
  return(gated_df)
}

#' post_process_gates
#' 
#' function to post process the events after model prediction.
#' @param gated_df dataframe with labels (third column).
#' @param n_cores Number of cores. Default to 1.
#' @param thr_dist  Distance threshold for centroids calculation. Default to 0.15.
#' @param include_zero  Consider centroid of label 0. Default to False.
#' @param remove_centroids  Remove centroids too near each other based on thr_dist value.
#' @param type  Type of post-processing.
#' @param concavity_val  Concavity of polygons for the "polygon" type of post-processing
#' @return Dataframe.
#' @export
#' @examples 
#' \donttest{post_process_gates()}

post_process_gates<-function(gated_df,n_cores=1,thr_dist=0.15,include_zero=F,remove_centroids=T,type="dist",
                             concavity_val=5,normalize_data=T){
  colnames(gated_df)<-c("x","y","classes")
  gated_df$classes<-as.character(gated_df$classes)
  if(type=="dist"){
    message("post-process based on events distance")
    new_df<-assign_events_to_nearest_centroids(gated_df = gated_df,
                                               n_cores = n_cores,thr_dist=thr_dist,
                                               include_zero = include_zero,remove_centroids = remove_centroids)
  }else if(type=="polygon"){
    message("post-process based on events distance checking polygons intersection")
    list_df_hull<-extract_polygon_gates(gated_df = gated_df,concavity_val=concavity_val)
    if(normalize_data==T){
      # check polygons intersections
      max_area_intersect<-check_polygons_intersection(list_df_hull = list_df_hull)
      message(sprintf("max_area_intersect:%f",max_area_intersect))
      # If 0: no polygon intersecion
      if(max_area_intersect>0.08){
        message("There is a relevant intersection")
        new_df<-assign_events_to_nearest_centroids(gated_df = gated_df,
                                                   n_cores = n_cores,thr_dist=0.05,
                                                   include_zero = F,remove_centroids = T)
      }else{
        new_df<-compute_gates(gated_df=gated_df,list_final_polygons_coords =  list_df_hull)
      }
    }else{
      new_df<-compute_gates(gated_df=gated_df,list_final_polygons_coords =  list_df_hull)
    }

  }
  return(new_df)
}

#' get_centroids
#' 
#' function to get centroids for each label
#' @param df dataframe with labels (third column).
#' @param low_thr Lower threshold for quantile calculation.
#' @param up_thr  Upper threshold for quantile calculation.
#' @param thr_dist  Distance threshold for centroids calculation. Default to 0.15.
#' @param include_zero  Consider centroid of label 0. Default to False.
#' @param remove_centroids  Remove centroids too near each other based on thr_dist value.
#' @return Dataframe.
#' @export
#' @examples 
#' \donttest{get_centroids()} 

get_centroids<-function(df,low_thr=0.10,up_thr=0.90,thr_dist=0.15,include_zero=F,remove_centroids=T){
  all_labels<-unique(df[,3])
  if(include_zero==F){
    all_labels<-all_labels[all_labels!="0"]
  }
  list_centroids_all_labels<-list()
  for(l in all_labels){
    #print(l)
    df_l<-df[which(df[,3]==l),]
    # get quantile lower and upper bound
    low_1<-as.numeric(quantile(df_l[,1], low_thr))
    up_1<-as.numeric(quantile(df_l[,1], up_thr))
    inds_1<-which(df_l[,1]>low_1 & df_l[,1]<up_1)
    
    low_2<-as.numeric(quantile(df_l[,2], low_thr))
    up_2<-as.numeric(quantile(df_l[,2], up_thr))
    inds_2<-which(df_l[,2]>low_2 & df_l[,2]<up_2)
    # get centroid of selected events
    mean_1<-mean(df_l[inds_1,1])
    mean_2<-mean(df_l[inds_2,2])
    vec_centroid_coord<-c(l,mean_1,mean_2)
    list_centroids_all_labels[[l]]<-vec_centroid_coord
  }
  df_centroids<-as.data.frame(do.call(rbind,list_centroids_all_labels),stringsAsFactors=T)
  colnames(df_centroids)<-c("label","mean_1","mean_2")
  df_centroids$label<-as.character(df_centroids$label)
  df_centroids$mean_1<-as.numeric(as.character(df_centroids$mean_1))
  df_centroids$mean_2<-as.numeric(as.character(df_centroids$mean_2))
  df_centroids$mean_1<-round(df_centroids$mean_1,2)
  df_centroids$mean_2<-round(df_centroids$mean_2,2)
  #---- remove NA values
  check_na_1<-is.nan(df_centroids$mean_1)
  inds<-which(check_na_1==T)
  if(length(inds)!=0){
    df_centroids<-df_centroids[-inds,]
  }
  check_na_2<-is.nan(df_centroids$mean_2)
  inds<-which(check_na_2==T)
  if(length(inds)!=0){
    df_centroids<-df_centroids[-inds,]
  }
  if(remove_centroids==T){
    # Remove centroids too near each other
    list_vec_dist<-list()
    for(c in 1:nrow(df_centroids)){
      current_label<-df_centroids$label[c]
      ref_centroid<-df_centroids[c,c("mean_1","mean_2")]
      vec_dist<-c()
      labels_test<-c()
      for(c2 in 1:nrow(df_centroids)){
        test_label<-df_centroids$label[c2]
        test_centroid<-df_centroids[c2,c("mean_1","mean_2")]
        m_coords<-rbind(ref_centroid,test_centroid)
        dist_value<-dist(m_coords,method = "euclidean")
        if(current_label==test_label){
          vec_dist<-c(vec_dist,1)
        }else{
          vec_dist<-c(vec_dist,dist_value)
        }
        labels_test<-c(labels_test,test_label)
      }
      names(vec_dist)<-labels_test
      list_vec_dist[[current_label]]<-vec_dist
    }
    df_dist<-do.call(rbind,list_vec_dist)
    all_labels_row<-row.names(df_dist)
    all_labels_col<-colnames(df_dist)
    label_to_remove<-c()
    label_to_remain<-c()
    for(c in 1:nrow(df_dist)){
      current_label<-all_labels_row[c]
      #print(current_label)
      vec_dist<-df_dist[c,]
      inds<-which(vec_dist<thr_dist)
      if(length(inds)!=0){
        labels_too_near<-all_labels_col[inds]
        inds_already_checked<-which((labels_too_near %in% label_to_remain)==T)
        if(length(inds_already_checked)!=0){
          labels_too_near<-labels_too_near[-inds_already_checked]
        }
        if(length(labels_too_near)!=0){
          label_to_remove<-c(label_to_remove,labels_too_near)
        }
      }
      final_check<-current_label %in% label_to_remove
      if(final_check==F){
        label_to_remain<-c(label_to_remain,current_label)
      }
    }
    label_to_remove<-unique(label_to_remove)
    inds_to_remove<-which((df_centroids$label %in% label_to_remove)==T)
    if(length(inds_to_remove)!=0){
      df_centroids<-df_centroids[-inds_to_remove,]
    }
  }
  return(df_centroids)
}

#' assign_events_to_nearest_centroids
#' 
#' function to assign events to class with nearest centroid.
#' @param gated_df dataframe with labels (third column).
#' @param n_cores Number of cores. Default to 1.
#' @param method_dist  Distance method calculation. Default to euclidean.
#' @param thr_dist  Distance threshold for centroids calculation. Default to 0.15.
#' @param include_zero  Consider centroid of label 0. Default to False.
#' @param remove_centroids  Remove centroids too near each other based on thr_dist value.
#' @return Dataframe.
#' @export
#' @examples 
#' \donttest{assign_events_to_nearest_centroids()} 


assign_events_to_nearest_centroids<-function(gated_df,n_cores=1,method_dist="euclidean",thr_dist=0.15,include_zero=F,remove_centroids=T){
  start<-Sys.time()
  df_centroids<-get_centroids(df = gated_df,thr_dist = thr_dist,include_zero = include_zero,remove_centroids = remove_centroids)
  list_new_classes<-parallel::mclapply(1:nrow(gated_df),function(i){
    message(i)
    coords_i<-gated_df[i,c(1,2)]
    colnames(coords_i)<-c("coord_1","coord_2")
    dist_vec<-c()
    labels<-c()
    for(c in 1:nrow(df_centroids)){
      current_l<-df_centroids$label[c]
      coord_centroid<-df_centroids[c,c("mean_1","mean_2")]
      colnames(coord_centroid)<-c("coord_1","coord_2")
      m_coords<-rbind(coords_i,coord_centroid)
      dist_value<-dist(m_coords,method = method_dist)
      dist_vec<-c(dist_vec,dist_value)
      labels<-c(labels,current_l)
    }
    ind_min<-which.min(dist_vec)
    if(length(ind_min)>1){
      ind_min<-ind_min[1]
    }
    selected_l<-labels[ind_min]
    return(selected_l)
  },mc.cores = n_cores)
  new_classes<-unlist(list_new_classes)
  gated_df$classes<-new_classes
  gated_df$classes<-as.character(gated_df$classes)
  end<-Sys.time()
  time_taken<-end-start
  message("Execution time:")
  message(time_taken)
  message("Done")
  return(gated_df)
}


#' smooth_hull
#' 
#' function to smooth concave polygons
#' @param hull_df Dataframe generate by  concaveman functions inside get_huget_hull_all_gates function
#' @param spar Spar value to regulate smoothing process: higher value (max 1) higher smoothing.
#' @return Dataframe.
#' @export
#' @examples 
#' \donttest{smooth_hull()} 


smooth_hull <- function(hull_df, spar = 0.7) {
  colnames(hull_df) <- c("x", "y")
  # Ensure it's closed loop
  hull_df <- rbind(hull_df, hull_df[1, ])
  smoothed <- data.frame(
    x = smooth.spline(hull_df$x, spar = spar)$y,
    y = smooth.spline(hull_df$y, spar = spar)$y
  )
  return(smoothed)
}