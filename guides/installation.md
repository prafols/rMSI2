### Installation

*rMSI2* provides a quite complex data model together with a graphical user interface (GUI), consequently rMSI2 depends on some other R packages that must be also installed. RGtk2 is one of these packages and is known to be problematic to install on Windows systems since it is not maintained anymore. Thus, we provide an installation guides for each operating system (MAC OSX is not available yet).

We are very aware of the issues related to the GUI installation for rMSI2 and the downsides of it and therefore, we are working to develop a new GUI model to overcome this issue. But, in the meanwhile, we provide this guide to facilitate a successful installation of rMSI2. 

Before installing rMSI2 ensure you have the RGtk2 related dependencies installed properly by following the guide below corresponding to your operating system.

* ##### Instructions for Windows users

RGtk2 is not available anymore as a binary package easy to install from CRAN to Windows systems. Building RGtk2 from source on Windows platform is rather complex and requires to have the Gtk2 SDK properly installed. In order to simplify the installation of rMSI2 we provide and R script to automate the RGtk2 installation on Windows based on archive binaries of the required packages and libraries. So, the steps below describe the recommended procedure to get rMSI2 working on a Windows machine:

1. The provided solution for Windows is only valid for a specific R version: 4.0.5. So, if you have another version of R installed we recommend to completely uninstall it before continuing this guide.
2. Get the R 4.0.5 for Windows from this link: [get R 4.0.5 for Windows](https://cloud.r-project.org/bin/windows/base/old/4.0.5/R-4.0.5-win.exe "Download R 4.0.5").
3. Start the installation of R 4.0.5 but be very careful to disable the installation of 32 bits version of R. The provided rMSI2 solution only works for 64 bits systems and having a mixed environment only claims for troubles. Thus, check out the following screen shoot of the R installer and ensure yours looks exactly the same (32 bits unchecked).

![alt text](https://github.com/prafols/rMSI2/blob/main/guides/images/Screenshot_installingR.png "R 4.0.5 Installation")

4. In Windows it is needed to install Rtools in order to be able to complie rMSI2 package from source. You can get it from here: TODO.. While installing Rtools ensure you select the option to set the system path (see the screenshot below) otherwise the installtion of rMSI2 may not work.

5. Finally we are ready to execute the script that is going to automatically pull the RGtk2 dependencies for Windows. Just run the following command in your R console:
 
 Then... source!

* ##### Instructions for Linux users
Gtk2 is a graphical library that comes pre-installed in almost any Linux distribution and therefore it is really straight forward to get rMSI2 working on a Linux machine. If you already have a working R environment on Linux, nothing special has to be made before installing rMSI2. No special version of R is needed on Linux, just continure reading to the next step.

* ##### Common installation of rMSI2 for all platforms
The simplest way to install rMSI2 and keep it updated is using devtools package. Install devtools from CRAN into your R session:
```R
>  install.packages("devtools")
```
Then simply tell devtools to install rMSI2 from github latest release:
```R
> devtools::install_github("prafols/rMSI2")
```

Please note that is a development version and no release has been made yet. So, keep looking at this page for future updates and releases.

