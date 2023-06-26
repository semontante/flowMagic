


#' get_test_sets
#' 
#' Function to obtain the data ready to be gated (validation set or test set).
#' It takes as input a flowSet and generates a list of dataframes with a ML structure for each ungated fcs file.
#' Each dataframe is the expression matrix of the root.
#' @param fs flowSet to gate.
#' @param gh Gating hierarchy.
#' @return List of dataframes.
#' @export
#' @examples 
#' \donttest{get_test_sets()} 

get_test_sets<-function(fs,gh){
  start<-Sys.time()
  n_samples<-length(fs)
  list_test_sets<-list()
  # extract names dimensions train set
  f_train<-gh_pop_get_data(gh)
  dimensions_train<-colnames(exprs(f_train))
  print("----- features train set -------")
  print(dimensions_train)
  # importing test samples
  print("------ processing test samples -------")
  for(i in 1:n_samples){
    f<-fs[[i]]
    sample_name<-identifier(f)
    root_exprs<-as.data.frame(f@exprs)
    # check dimensions of the test data (they must be equal to the gh dimensions)
    dimensions_test<-colnames(root_exprs)
    check_length <-length(dimensions_train) == length(dimensions_test)
    if(check_length==F){
      stop("number of test dimensions != number of train dimensions")
    }
    check_names <-dimensions_train == dimensions_test
    if(all(check_names)==F){
      print("name dimensions train samples != name dimensions test data. Fixing in process..." )
      print(dimensions_test)
      inds<-which(check_names==F)
      inds_correct_dimensions<-which(check_names==T)
      
      test_names_to_fix<-dimensions_test[inds]
      # check if test names are different
      for(test_name in test_names_to_fix){
        ind<-grep(test_name,dimensions_train,fixed = T)
        if(length(ind)>1){
          stop(sprint("test dimension %s multiple matches in train dimensions.Check dimensions sample %s",test_name,sample_name))
        }else if(length(ind)==0){
          stop(sprint("test dimension %s not found in train dimensions.Check dimensions sample %s",test_name,sample_name))
        }
      }
      root_exprs<-root_exprs[,dimensions_train]
    }
    dimensions_test<-colnames(root_exprs)
    check_names <-dimensions_train == dimensions_test
    if(all(check_names)==F){
      print(dimensions_test)
      stop(sprintf("dimensions still different. Fixing Failed. Check sample: %s",sample_name))
    }
    # generate list entry
    list_test_sets[[sample_name]]<-root_exprs
  }
  gc()
  end<-Sys.time()
  time_taken<-end-start
  print("Time of execution:")
  print(time_taken)
  print("Done")
  return(list_test_sets)
}
