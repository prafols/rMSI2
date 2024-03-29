% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/librMSIdata.R
\name{Trim_imzMLData}
\alias{Trim_imzMLData}
\title{Trim_imzMLData.}
\usage{
Trim_imzMLData(
  input_data_file,
  output_data_file,
  keepPixelIDs = NULL,
  mass_min = NULL,
  mass_max = NULL,
  intensity_trim = 0
)
}
\arguments{
\item{input_data_file}{complete path to an input imzML file.}

\item{output_data_file}{complete path to an output imzML file (will be overwritten).}

\item{keepPixelIDs}{the pixel ID's to kept in the output imzML file.}

\item{mass_min}{the minimum value of the output mass range.}

\item{mass_max}{the maximum value of the output mass range.}

\item{intensity_trim}{trim intensities below this value.}
}
\description{
Create a partial copy of an imZML file by coping only the data that matches the mass and pixel ID filters.
}
