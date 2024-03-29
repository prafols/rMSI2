% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Process.R
\name{FormatPeakMatrix}
\alias{FormatPeakMatrix}
\title{FormatPeakMatrix.
Formats a C style peak matrix generated by MTPeakPicking::BinPeaks() to a rMSIprocPeakMatrix.}
\usage{
FormatPeakMatrix(cPeakMatrix, posMat, numPixels, names, uuid, posMotors)
}
\arguments{
\item{cPeakMatrix}{a peak matrix with the same format as retured by MTPeakPicking::BinPeaks().}

\item{posMat}{a rMSI image pos matrix.}

\item{numPixels}{a vector including the number of pixels of each sample.}

\item{names}{a vector of strings with the name of each sample.}

\item{uuid}{a vector of img UUID to be also stored in peak matrices}

\item{posMotors}{a rMSI image original motros coordinates matrix.}
}
\value{
the formated matrix.
}
\description{
FormatPeakMatrix.
Formats a C style peak matrix generated by MTPeakPicking::BinPeaks() to a rMSIprocPeakMatrix.
}
