shinyÉPICo (documentation)
================

<!-- README.md is generated from README.Rmd. Please edit that file -->

<!-- badges: start -->

[![License: AGPL
v3](https://img.shields.io/badge/License-AGPL%20v3-blue.svg)](https://www.gnu.org/licenses/agpl-3.0)
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

R 3.6 or higher with updated packages is required.

Since the application allows to follow interactively all the analysis
process, many objects have to be stored in RAM memory. Therefore, the
application can be quite memory demanding, especially when trying to
analyze a large number of samples. We recommend at least 12GB of RAM for
a smooth use of the application, but depending on the number of samples
analyzed and whether they are EPIC or 450k, the needs may be lower or
higher.

# Installation and use

To install shinyÉPICo, you have to use the GitHub repository. It is easy
to install it directly in R using the install\_github function from the
remotes package.

``` r
install.packages("remotes")
library("remotes")
install_github("omorante/shinyepico", upgrade="always")
```

To run the package:

``` r
library("shinyepico")
run_shinyepico()
```

It is highly recommended to update all packages when install\_github
asks to avoid problems during the use of the application. The minimum
versions required are:

    DT (>= 0.15.0),
    IlluminaHumanMethylation450kanno.ilmn12.hg19,
    IlluminaHumanMethylation450kmanifest,
    IlluminaHumanMethylationEPICanno.ilm10b4.hg19,
    IlluminaHumanMethylationEPICmanifest,
    data.table (>= 1.13.0),
    doParallel (>= 1.0.0),
    dplyr (>= 1.0.0),
    foreach (>= 1.5.0),
    ggplot2 (>= 3.3.0),
    gplots (>= 3.0.0),
    heatmaply (>= 1.1.0),
    limma (>= 3.44.0),
    minfi (>= 1.34.0),
    plotly (>= 4.9.2),
    reshape2 (>= 1.4.0),
    rlang (>= 0.4.0),
    rmarkdown (>= 2.3.0),
    shiny (>= 1.5.0),
    shinyWidgets (>= 0.5.0),
    shinycssloaders (>= 0.3.0),
    shinyjs (>= 1.1.0),
    shinythemes (>= 1.1.0),
    statmod (>= 1.4.0),
    tidyr (>= 1.1.0)

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

You can find examples of both the sample sheet and the final .zip to
upload to the application in the example\_data folder.
