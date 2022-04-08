#This script automatizes the rMSI and rMSIproc installation process on Windows machines where it is difficult to get RGtk2 binaries

if(version$major != "4")
{
  stop("You need R version 4.1 to run this!")
}

installPackage_gdrive <- function(package_name, package_URL)
{
  zipfile <- file.path(tempdir(), paste0(package_name,".zip"))
  download.file(package_URL, zipfile, method = "libcurl", mode = "wb")
  install.packages( zipfile , repos = NULL)
}

#Install RGtk2
installPackage_gdrive("RGtk2", "https://drive.google.com/uc?export=download&id=1yM1lAJHtpTKylgRLBjtvl__-zaX7HR4C")
gtkzipfile <- tempfile(fileext = ".zip")
download.file("https://drive.google.com/uc?export=download&id=17H96joCZ1wHg_T2rbkTMTq4_AnJMrQ6-", gtkzipfile, method = "libcurl", mode = "wb")
gtkPath <- file.path(find.package("RGtk2"), "gtk", "x64")
dir.create(gtkPath, recursive = T)
unzip(gtkzipfile, exdir = gtkPath)

#Install gWidgets2RGtk2
installPackage_gdrive("gWidgets2RGtk2","https://drive.google.com/uc?export=download&id=1SO2lTLeMge52FmFGM3x5KsCQqoS0eT65")

#Install cairoDevice
installPackage_gdrive("cairoDevice","https://drive.google.com/uc?export=download&id=1OdVZcRAU90EycpEs88KH0YCUuPoDc0AM")

#Install the devtools
install.packages("devtools")

#Install rMSI
devtools::install_github("prafols/rMSI", ref = "0.9.1")

#Install rMSIproc
devtools::install_github("prafols/rMSIproc", ref = "0.3.1")
