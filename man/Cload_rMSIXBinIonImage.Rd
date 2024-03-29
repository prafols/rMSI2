% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{Cload_rMSIXBinIonImage}
\alias{Cload_rMSIXBinIonImage}
\title{Cload_rMSIXBinIonImage.}
\usage{
Cload_rMSIXBinIonImage(
  rMSIobj,
  ionIndex,
  ionCount,
  normalization_coefs,
  number_of_threads
)
}
\arguments{
\item{normalization_coefs}{a vector containing the intensy normalization coeficients.}

\item{rMSIobj:}{an rMSI object prefilled with a parsed imzML.}

\item{ionIndex:}{the first mass channel at which the image starts.}

\item{ionCount:}{the numer of mass channels used to construct the ion image (a.k.a. image tolerance window).}

\item{number_of_threads:}{number of threads used for imgStream encoding.}
}
\value{
the ion image as a NumericMatrix using max operator with all the ion images of the mass channels.
}
\description{
loads a ion image from the .BrNSI img stream.
}
