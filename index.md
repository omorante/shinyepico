# shinyÉPICo

![](https://github.com/omorante/shinyepico/blob/master/inst/images/logo.png)

ShinyÉPICo is a web interface based on Shiny that makes it easy to do differentially methylated CpGs analysis from Illumina EPIC or 450k DNA methylation arrays. This program allows following a standard pipeline of normalization (with minfi package), model creation and statistical analysis (with limma package), with different options in each step and plots to be able to choose properly. Moreover, you can select different options in the final heatmap and download an RMarkdown report with all the steps chosen.

## System Requirements

ShinyÉPICo can run in GNU/Linux, Windows or macOS. The package dependencies are automatically tried to install when you install the package.

Since the application allows to follow interactively all the analysis process, many objects have to be stored in RAM memory. Therefore, the application can be quite memory demanding, especially when trying to analyze a large number of samples. We recommend at least 12GB of RAM for a smooth use of the application, but depending on the number of samples analyzed and whether they are EPIC or 450k, the needs may be lower or higher.

## Installation and use

To install shinyÉPICo, you have to use the GitHub repository. It is easy to install it directly in R using the install_github function from the remotes package. It is highly recommended to update all packages when install_github asks:

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

You can assign the number of cores or change the upload memory limit with the arguments of that function. The parallelization of this application does not have a great impact on the RAM memory consumption and it is recommended even in low-end computers.

## Input data

The first step in the shinyÉPICo workflow is prepare the data in the properly format. iDAT files should be compressed in a .zip file. The name of the files should follow the standard convention: **XXXXXXXXXXXX_YYYYYY_ZZZ.idat** being XXXXXXXXXXXX the Sentrix_ID, YYYYYY the Sentrix_Position and ZZZ Grn or Red (corresponding, respectively, to the Red and Green signal file).

Moreover, a CSV (comma-separated) file with the annotation of the experiment should be included. It is mandatory to include the Sentrix_ID and Sentrix_Position columns that allows the software to find their respective iDAT files. Moreover, other columns should be added to reflect the different variables (e.g. sample name, health/disease, treatment/control, age, sex, hybridization day, etc.). Usually, iDATs and sample sheet are obtained meeting these features by default, and no additional work is required.

You can find examples of both the sample sheet and the final .zip to upload to the application in the example_data folder.
