## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----eval=FALSE---------------------------------------------------------------
#  install.packages("remotes")
#  library("remotes")
#  install_github("omorante/shinyepico", upgrade="always")

## ----eval=FALSE---------------------------------------------------------------
#  library("shinyepico")
#  set.seed(123)
#  run_shinyepico()

## ---- echo=FALSE, eval=TRUE---------------------------------------------------
knitr::kable(read.csv("./Sample_Sheet.csv"))

## ---- echo=FALSE, eval=TRUE---------------------------------------------------
knitr::kable(data.frame(MAC = c(1,0,1,0,0,1), MO = c(0,1,0,1,1,0), DonorB = c(0,1,0,0,0,1), DonorC=c(0,0,1,0,1,0)))

## ---- echo=FALSE, eval=TRUE---------------------------------------------------
knitr::kable(data.frame(contrast = "MAC-MO", Hypermethylated=2272, Hypomethylated=21, total=2293))

## ---- echo=FALSE, eval=TRUE---------------------------------------------------
knitr::kable(data.frame(contrast = rep("MAC-MO",3), Hypermethylated=c(54,112,48), Hypomethylated=c(76,95,56),total=c(130,207,104)))

## ----echo=FALSE---------------------------------------------------------------
sessionInfo()

