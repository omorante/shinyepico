library(minfiData)
library(data.table)

baseDir <- system.file("extdata", package = "minfiData")
sample_sheet <- minfi::read.metharray.sheet(baseDir)


#Testing reading iDATs
rgset = read_idats(sample_sheet)
test_that("iDAT reading",{
  expect_is(rgset, "RGChannelSet")
  expect_equal(length(minfi::featureNames(rgset)), 622399)
})


#Testing clean sample sheet calculation
clean_sample_sheet <- generate_clean_samplesheet(sample_sheet, donorvar = "person")

test_that("Clean Sample Sheet",{
  expect_equal(ncol(clean_sample_sheet), 7)
  expect_equal(nrow(clean_sample_sheet), 6)
  expect_equal(unlist(lapply(clean_sample_sheet, class), use.names=FALSE), c("factor","factor","numeric", "factor","factor","factor","factor"))
  
})

#Testing normalization
gset <- normalize_rgset(rgset, "Quantile", 0.01, TRUE, 0, TRUE, TRUE)

test_that("Normalization",{
  expect_error(normalize_rgset(rgset, "Unexistent", 0.01, TRUE, 0, TRUE, TRUE))
  expect_error(normalize_rgset(rgset, "Unexistent", 0.01, "TRUE", 0, TRUE, TRUE))
  expect_error(normalize_rgset(rgset, "Unexistent", "0.01", TRUE, 0, TRUE, TRUE))
  expect_is(gset, "GenomicRatioSet")
  expect_equal(length(minfi::featureNames(gset)), 443192)
})

#Testing Design
design <- generate_design(voi = "status", sample_sheet = clean_sample_sheet, sample_name = "Sample_Name", covariables = c("age","sex"), interactions = NULL)

test_that("Design",{
  expect_error(generate_design(voi = c("status","error"), sample_sheet = clean_sample_sheet, sample_name = "Sample_Name", covariables = c("age","sex"), interactions = NULL))
  expect_error(generate_design(voi = "status", sample_sheet = c("Incorrect_Type"), sample_name = "Sample_Name", covariables = c("age","sex"), interactions = NULL))
  expect_is(design, "matrix")
  expect_equal(ncol(design), 4)
  expect_equal(nrow(design), 6)
})

#Testing limma model
contrasts <- generate_contrasts(as.factor(clean_sample_sheet[["status"]]))
Mvalues <- minfi::getM(gset)
Bvalues <- as.data.frame(minfi::getBeta(gset))
colnames(Bvalues) <- sample_sheet[["Sample_Name"]]
fit <- generate_limma_fit(Mvalues, design, weighting = FALSE)

test_that("Limma Model",{
  expect_is(fit, "MArrayLM")
})

#Testing global difs
global_difs <- calculate_global_difs(Bvalues, as.factor(clean_sample_sheet[["status"]]), contrasts, 1)

test_that("Global Difs",{
  expect_is(global_difs, "data.frame")
  expect_equal(global_difs[order(global_difs$dif_cancer_normal),]$cpg[1],"cg14642259") 
  expect_equal(global_difs[order(global_difs$dif_cancer_normal),]$cpg[200], "cg08756887")
})


#Testing DMPs calculation
ebayes_tables <- find_dif_cpgs(design, fit, contrasts, FALSE, FALSE, 1)

test_that("DMP Calculation",{
  expect_is(ebayes_tables, "list")
  expect_equal(length(ebayes_tables), 1)
  expect_equal(ebayes_tables[[1]][order(ebayes_tables[[1]][["adj.P.Val"]]),]$cpg[1000], "cg03203223")
  expect_equal(ebayes_tables[[1]][order(ebayes_tables[[1]][["adj.P.Val"]]),]$cpg[20000], "cg00363312")
})

#Testing Filtered List
filtered_list <- create_filtered_list(ebayes_tables, global_difs, 0.2, 0.05, 1, 1)

test_that("DMP Calculation",{
  expect_is(filtered_list, "list")
  expect_equal(length(filtered_list), 1)
  expect_equal(nrow(filtered_list[[1]]), 9227)
  expect_equal(nrow(filtered_list[[1]][filtered_list[[1]]$dif_current < 0,]), 3277)
})

#Testing DMR calculation

set.seed(123) #seed for the calculation

mcsea_list <- find_dmrs(ebayes_tables, 5, "450k", as.factor(clean_sample_sheet[["status"]]), "promoters", contrasts, Bvalues, 5000, 1)
mcsea_list <- add_dmrs_globaldifs(mcsea_list, global_difs, "promoters")
mcsea_filtered <- filter_dmrs(mcsea_list, 0.05,1,0.05,"promoters",contrasts)

test_that("DMR Calculation",{
  expect_is(mcsea_list, "list")
  expect_is(mcsea_filtered, "list")
  expect_equal(length(mcsea_list), 1)
  expect_equal(length(mcsea_filtered), 1)
  expect_equal(nrow(mcsea_filtered[[1]][["promoters"]]), 1402)
  expect_equal(row.names(mcsea_filtered[[1]][["promoters"]])[1], "FASTK")
})



