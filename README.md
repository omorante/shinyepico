shinyÉPICo
================

<!-- README.md is generated from README.Rmd. Please edit that file -->

<!-- badges: start -->

[![License: AGPL
v3](https://img.shields.io/badge/License-AGPL%20v3-blue.svg)](https://www.gnu.org/licenses/agpl-3.0)
![Lines of code](https://img.shields.io/tokei/lines/github/omorante/shinyepico)
<!-- badges: end -->

<img src="https://github.com/omorante/shinyepico/blob/master/inst/images/logo.png" width="200px" />

# Description

ShinyÉPICo is a web interface based on Shiny that makes it easy to do
differentially methylated CpGs analysis from Illumina EPIC or 450k DNA
methylation arrays. This program allows following a standard pipeline of
normalization (with minfi package), model creation and statistical
analysis (with limma package), with different options in each step and
plots to be able to choose properly. Moreover, you can select different
options in the final heatmap and download an RMarkdown report with all
the steps chosen.

# System Requirements

ShinyÉPICo can run in GNU/Linux, Windows or macOS. The package
dependencies are automatically tried to install when you install the
package.

**R 4.0** or higher with updated packages is required.

Since the application allows to follow interactively all the analysis
process, many objects have to be stored in RAM memory. Therefore, the
application can be quite memory demanding, especially when trying to
analyze a large number of samples. We recommend at least 12GB of RAM for
a smooth use of the application, but depending on the number of samples
analyzed and whether they are EPIC or 450k, the needs may be lower or
higher.

# Installation and use

Our package is now included in the 3.13 Release of [Bioconductor](https://bioconductor.org/packages/release/bioc/html/shinyepico.html).

To install the package through Bioconductor you have to first install R 4.1 and enter:

``` r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("shinyepico")
```

To run the package:

``` r
library("shinyepico")
run_shinyepico()
```

# Docker container

In addition to the shinyÉPICo package, now an official shinyÉPICo container is available in [Docker](https://hub.docker.com/repository/docker/omorante/shinyepico). This is very useful to run shinyÉPICo in a webserver. [Shinyproxy](https://www.shinyproxy.io/) is an excellent option to provide ShinyApps to several users from a private or public server.

# Shiny App

The recommended way to run shinyÉPIco locally is through Bioconductor, as previously explained. However, if you would like to configure a Shiny Server with shinyÉPICo and Docker/Shinyproxy is not an option, it is possible to use shinyÉPICo as a shinyApp. It can be downloaded from the [shiny_app branch](https://github.com/omorante/shinyepico/tree/shiny_app). However, the Shiny Server should be properly configured, including all the packages needed for shinyÉPIco. An easy way to do that is install the shinyepico package from Bioconductor.

# Parameters

The function run\_shinyepico has 4 parameters that can be modified:

  - **n\_cores:** The application is partially compatible with parallel
    computing. This numeric parameter controls the number of cores and,
    by default, it is half of the detected cores. The parallelization of
    this application does not have a great impact on the RAM memory
    consumption and it is recommended even in low-end computers.
  - **max\_upload\_size:** By default, shiny applications have an upload
    limit (in MB), useful when the application is running in a web
    server. By default, this parameter is 2000MB.
  - **host:** IP used to deploy the application. By default, this
    parameter is your local IP (127.0.0.1) which means that only you,
    from your computer, will have access to the application. However, it
    is possible to make the app reachable to other computers in the same
    LAN changing the IP to 0.0.0.0.
  - **port:** Port used to deploy the application. By default, a random
    free port.

# Input data

The first step in the shinyÉPICo workflow is prepare the data in the
properly format. iDAT files should be compressed in a .zip file. The
name of the files should follow the standard convention:
**XXXXXXXXXXXX\_YYYYYY\_ZZZ.idat** being XXXXXXXXXXXX the Sentrix\_ID,
YYYYYY the Sentrix\_Position and ZZZ Grn or Red (corresponding,
respectively, to the Red and Green signal file).

Moreover, a CSV (comma-separated) file with the annotation of the
experiment should be included. It is mandatory to include the
Sentrix\_ID and Sentrix\_Position columns that allows the software to
find their respective iDAT files. Moreover, other columns should be
added to reflect the different variables (e.g. sample name,
health/disease, treatment/control, age, sex, hybridization day, etc.).
Usually, iDATs and sample sheet are obtained meeting these features by
default, and no additional work is required.

You can find an example of the final .zip to
upload to the application in [the example\_data folder of the doc branch](https://github.com/omorante/shinyepico/blob/doc/example_data/Li_NAR_2019.zip).

# Package documentation

The package includes a detailed vignette explaining the steps and options of the application. In addition, you can access the updated version of the vignette in pdf at this [link](https://omorante.github.io/shinyepico/shiny_epico.pdf)


