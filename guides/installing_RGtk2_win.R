#This script automatizes the RGtk2 installation process on Windows machines.
if(version$major != "4" && verison.minor != "0.5")
{
  stop("You need R version 4.0.5 to run this!")
}

installPackage_gdrive <- function(package_name, package_URL)
{
  zipfile <- file.path(tempdir(), paste0(package_name,".zip"))
  download.file(package_URL, zipfile, method = "libcurl", mode = "wb")
  install.packages( zipfile , repos = NULL)
}

#Install RGtk2
installPackage_gdrive("RGtk2", "https://drive.google.com/uc?export=download&id=184mBmTeLEOOcKM2LAX0YCSDIHkHt2WXn")
gtkzipfile <- tempfile(fileext = ".zip")
download.file("https://drive.google.com/uc?export=download&id=1OfUoXPi1Gbzgl8dPjp8kVWXYM7rM5Tlq", gtkzipfile, method = "libcurl", mode = "wb")
gtkPath <- file.path(find.package("RGtk2"), "gtk", "x64")
dir.create(gtkPath, recursive = T)
unzip(gtkzipfile, exdir = gtkPath)

#Install gWidgets2RGtk2
installPackage_gdrive("gWidgets2RGtk2","https://drive.google.com/uc?export=download&id=13K5TaqH1tKurDpEPYDo5XAYY0chPVbfX")

#Install cairoDevice
installPackage_gdrive("cairoDevice","https://drive.google.com/uc?export=download&id=1N_BIKyidiEc1CbdoGOoEdwLs-nL4eREW")
