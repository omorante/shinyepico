shinyÉPICo
================

<!-- README.md is generated from README.Rmd. Please edit that file -->

<!-- badges: start -->

[![License: AGPL
v3](https://img.shields.io/badge/License-AGPL%20v3-blue.svg)](https://www.gnu.org/licenses/agpl-3.0)
<!-- badges: end -->

<img src="/home/omorante/Documentos/R_projects/master/inst/images/logo.png" width="200px" />

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

Since the application allows to follow interactively all the analysis
process, many objects have to be stored in RAM memory. Therefore, the
application can be quite memory demanding, especially when trying to
analyze a large number of samples. We recommend at least 16GB of RAM for
a smooth use of the application, but depending on the number of samples
analyzed and whether they are EPIC or 450k, the needs may be lower or
higher.

# Installation and use

To install shinyÉPICo, you have to use the GitHub repository. It is easy
to install it directly in R using the install\_github function from the
remotes package. It is highly recommended to update all packages when
install\_github asks:

``` r
install.packages("remotes")
library("remotes")
install_github("omorante/shinyepico")
```

To run the package:

``` r
library("shinyepico")
run_shinyepico()
```

You can assign the number of cores or change the upload memory limit
with the arguments of that function. The parallelization of this
application does not have a great impact on the RAM memory consumption
and it is recommended even in low-end computers.
