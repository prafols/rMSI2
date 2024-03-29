% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Process.R
\name{RunPreProcessing}
\alias{RunPreProcessing}
\title{RunPreProcessing}
\usage{
RunPreProcessing(
  proc_params,
  output_data_path,
  img_lst,
  data_is_peaklist,
  numOfThreads = min(parallel::detectCores()/2, 6),
  memoryPerThreadMB = 200
)
}
\arguments{
\item{proc_params}{a ProcParams object containing the processing parameters.}

\item{output_data_path}{output path to store the processing results.}

\item{img_lst}{a rMSI objects lis to process when processing spectral data or reference to peaks list when processing peak lists.}

\item{data_is_peaklist}{a boolean indicating wheter the imzML data contains peak lists instead of spectral data.}

\item{numOfThreads}{the number number of threads used to process the data.}

\item{memoryPerThreadMB}{maximum allowed memory by each thread. The total number of trehad will be two times numOfThreads, so the total memory usage will be: 2*numOfThreads*memoryPerThreadMB.}
}
\description{
Process a single image or multiple images with the complete processing workflow.
}
