% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/post_processing.R
\name{get_centroids}
\alias{get_centroids}
\title{get_centroids}
\usage{
get_centroids(
  df,
  low_thr = 0.1,
  up_thr = 0.9,
  thr_dist = 0.15,
  include_zero = F,
  remove_centroids = T
)
}
\arguments{
\item{df}{dataframe with labels (third column).}

\item{low_thr}{Lower threshold for quantile calculation.}

\item{up_thr}{Upper threshold for quantile calculation.}

\item{thr_dist}{Distance threshold for centroids calculation. Default to 0.15.}

\item{include_zero}{Consider centroid of label 0. Default to False.}

\item{remove_centroids}{Remove centroids too near each other based on thr_dist value.}
}
\value{
Dataframe.
}
\description{
function to get centroids for each label
}
\examples{
\donttest{get_centroids()} 
}
