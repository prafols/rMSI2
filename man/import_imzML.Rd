% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/imzMLreader.R
\name{import_imzML}
\alias{import_imzML}
\title{import_imzML.}
\usage{
import_imzML(
  imzML_File,
  ibd_File = paste(sub("\\\\.[^.]*$", "", imzML_File), ".ibd", sep = ""),
  fun_progress = NULL,
  fun_text = NULL,
  close_signal = NULL,
  verifyChecksum = F,
  convertProcessed2Continuous = T,
  subImg_rename = NULL,
  subImg_Coords = NULL,
  fixBrokenUUID = F
)
}
\arguments{
\item{imzML_File}{full path to .imzML file}

\item{ibd_File}{path to the binary file (default the same as imzML file but with .ibd extension)}

\item{fun_progress}{This is a callback function to update the progress of loading data. See details for more information.}

\item{fun_text}{This is a callback function to update the label widget of loading data. See details for more information.}

\item{close_signal}{function to be called if loading process is abored.}

\item{verifyChecksum}{if the binary file checksum must be verified, it is disabled by default for convenice with really big files.}

\item{convertProcessed2Continuous}{if true (the default) an imzML file in processed mode will be converted to a continuous mode.}

\item{subImg_rename}{alternative image name, new rMSI files will be created with the given name.}

\item{subImg_Coords}{a Complex vector with the motors coordinates to be included in the rMSI data.}

\item{fixBrokenUUID}{set to FALSE by default to automatically fix an uuid mismatch between the ibd and the imzML files (a warning message will be raised).

 Imports an imzML image to an rMSI data object.
 It is recomanded to use rMSI2::LoadMsiData directly instead of this function.}
}
\value{
an rMSI data object.
}
\description{
import_imzML.
}
