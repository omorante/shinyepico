
<!-- README.md is generated from README.Rmd. Please edit that file -->
ShinyÉPICo
==========

<!-- badges: start -->
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental) <!-- badges: end -->

ShinyÉPICO is a web interface based on Shiny that makes it easy to do differentially methylated CpGs analysis from Illumina EPIC methylation arrays. This program allows following a standard pipeline of normalization (with minfi package), model creation and statistical analysis (with limma package), with different options in each step and plots to be able to choose properly. Moreover, you can select different options in the final heatmap and download an RMarkdown report with all the steps chosen.

Installation and use
--------------------

To install ShinyEPICO, you have to use the GitHub repository. It is easy to install it directly in R using the install\_github function from the remotes package:

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

You can assign the number of cores or the upload memory limit with the arguments of that function. The parallelization makes the application faster, but it requires more RAM available. By default, n\_cores is 1.
