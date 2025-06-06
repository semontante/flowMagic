% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/training.R
\name{magicTrain_knn}
\alias{magicTrain_knn}
\title{magicTrain_knn}
\usage{
magicTrain_knn(
  Xtrain,
  Ytrain,
  k_cv = 10,
  list_index_train = NULL,
  list_index_val = NULL,
  tune_lenght = 5
)
}
\arguments{
\item{Xtrain}{Dataframe of training features.}

\item{Ytrain}{Dataframe of labels (one column).}

\item{k_cv}{Number of k for cross-validation.}

\item{list_index_train}{List of vector of indices to use in training for each fold.}

\item{list_index_val}{List of vector of indices to use as held out data for each fold.}

\item{tune_lenght}{Number of hyper parameters to test. Default to 5.}
}
\value{
model object.
}
\description{
function to generate a random forest training model.
}
\examples{
\donttest{magicTrain_knn()}
}
