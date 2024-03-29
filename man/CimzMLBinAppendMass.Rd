% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{CimzMLBinAppendMass}
\alias{CimzMLBinAppendMass}
\title{CimzMLBinAppendMass.
This function appends a new mass axis to a given ibd file.
The last added offset is returned.}
\usage{
CimzMLBinAppendMass(ibdFname, mz_dataTypeString, mzNew)
}
\arguments{
\item{ibdFname:}{full path to the ibd file.}

\item{mz_dataTypeString:}{String to specify the data format used to encode m/z values.}

\item{mzNew:}{The mass axis to append.}
}
\description{
CimzMLBinAppendMass.
This function appends a new mass axis to a given ibd file.
The last added offset is returned.
}
