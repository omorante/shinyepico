# This branch of shinyepico is intended to use with ShinyServer. The server must be properly configured and included all the packages needed by shinyepico. 
# An easy way to do it is install the shinyepico package from Bioconductor

library(shiny)
library(data.table)
library(shinyWidgets)
library(statmod)
library(rlang)

source('R/app_ui.R')
source('R/app_server.R')
source('R/utils_analysis.R')
source('R/utils_download.R')
source('R/utils_graphs.R')

shinyOptions(shiny.maxRequestSize = 10000 * 1024^2, 
             n_cores = 1)

shinyApp(app_ui, app_server)

