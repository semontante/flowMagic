% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/importing.R
\name{get_train_data}
\alias{get_train_data}
\title{get_train_data}
\usage{
get_train_data(
  paths_file = NULL,
  df_paths = NULL,
  n_cores = 1,
  prop_down = NULL,
  remove_class = NULL,
  n_points_per_plot = NULL,
  normalize_data = T,
  vec_col = NULL
)
}
\arguments{
\item{paths_file}{Vector of paths. Each path points toward a single csv file containin training info (labels and bivariate expression). paths_file can be also directly the list of dataframes containing labels and bivariate expression.}

\item{df_paths}{Dataframe containing the paths of file to read. The paths to data must be in the first column. The associated paths to classes are in second column.}

\item{n_cores}{Number of cores. Default to 1.}

\item{prop_down}{Proportion of events (downsampling). Default to NULL (downsampling using number of points).}

\item{remove_class}{Vector of classes to ignore. Default to NULL.}

\item{n_points_per_plot}{Number of points for downsampling.}

\item{normalize_data}{If True, data is normalized to 0-1 range. Default to True.}

\item{vec_col}{vector of columns names if the input dataframes have more than 3 columns. The third column name must always refer to the column with the gate label of each event. Default to NULL.}
}
\value{
Dataframe.
}
\description{
function to import training data based on paths to files.
}
\examples{
\donttest{get_train_data()}
}
