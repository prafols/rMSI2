% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/libimgramdisk.R
\name{getIdsFromCoords}
\alias{getIdsFromCoords}
\title{Obtain the image Identifiers from a given set of images coords}
\usage{
getIdsFromCoords(Img, Coords)
}
\arguments{
\item{Img}{the rMSI object where the data is stored (ramdisk).}

\item{Coords}{a coords vector of spectra to load represented as complex numbers where real part corresponds to X and imaginary to Y.}
}
\value{
a vector of identifiers corresponding to given coords.

Calculate the pixel identifiers from the pixel coordinates.
}
\description{
Obtain the image Identifiers from a given set of images coords
}
