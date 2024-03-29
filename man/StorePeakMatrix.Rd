% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/LoadStorePeakMat.R
\name{StorePeakMatrix}
\alias{StorePeakMatrix}
\title{StorePeakMatrix.}
\usage{
StorePeakMatrix(data_path, data)
}
\arguments{
\item{data_path}{full path including filename where data must be stored.}

\item{data}{a rMSI2 peak list object of the rMSIprocPeakMatrix class.}
}
\description{
Stores a binned peaks matrix to HDD.
Data is stored compressed using RData format with .pkmat extension.
}
