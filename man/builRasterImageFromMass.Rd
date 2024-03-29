% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/libimgramdisk.R
\name{builRasterImageFromMass}
\alias{builRasterImageFromMass}
\title{Build a image slice from specified datacube masses.}
\usage{
builRasterImageFromMass(Img, Mass, Tolerance, Normalization = NULL)
}
\arguments{
\item{Img}{the rMSI object where the data is stored (ramdisk).}

\item{Mass}{the slected mass.}

\item{Tolerance}{a tolerance expressed in daltons around the slected mass.}

\item{Normalization}{optionally a vector of the normalization factor for each pixel.}
}
\value{
list with a matrix with the same size as image size containing the pixel values, used mass and tolerance.
}
\description{
Builds a image from the selected masses in the rMSI object. The image is returned arranged in a matrix containing each pixel value.
}
