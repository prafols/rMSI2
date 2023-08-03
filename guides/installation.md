### Installation

*rMSI2* provides a quite complex data model together with an optional graphical user interface (GUI), consequently rMSI2 depends on some other R packages that must be also installed. Previous versions of rMSI2 relied on the RGtk2 package to drawn its GUI. This GTK dependency was known to be problematic to install on Windows and MAC systems. Currently, we are proud to provide a rewritten GUI that do not depend on GTK anymore. Instead, the renewed GUI uses the R built-in TCL/TK library to provide the required graphics functions. This allows a much easier installation on any R platform. 
Follow the instructions below corresponding to your operating system to install the dependencies needed for the GUI .

* ##### Instructions for Windows users
To install rMSI2 from sources, we recommend a recent version of R (from 4.2 and above) and the corresponding Rtools. You can obtain the Rtools matching your R version here: [get Rtools](https://cran.r-project.org/bin/windows/Rtools/ "Rtools"). Next, you can proceed with the installation instructions below.

* ##### Common installation of rMSI2 for all platforms
The simplest way to install rMSI2 and keep it updated is using devtools package. Install devtools from CRAN into your R session:
```R
>  install.packages("devtools")
```
Then tell devtools to install rMSI2 from github latest release:
```R
> devtools::install_github("prafols/rMSI2")
```

Please note that is a development version and no release has been made yet. So, keep looking at this page for future updates and releases.

