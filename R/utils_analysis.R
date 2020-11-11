globalVariables("cont")
globalVariables("contrast")
globalVariables("grupo")
globalVariables("type")
globalVariables("table")
globalVariables("ranking")
globalVariables("association")

# READING AND NORMALIZATION FUNCTIONS
read_idats <- function(targets, detectP) {
  RGSet <- minfi::read.metharray.exp(
    targets = targets, verbose = TRUE,
    force = TRUE
  )

  RGSet[(rowMeans(as.matrix(minfi::detectionP(RGSet))) < detectP), ]
}

generate_clean_samplesheet <- function(target_samplesheet, donorvar){
  
  suppressWarnings({
    clean_sample_sheet <- as.data.frame(lapply(target_samplesheet, function(x) {
      if (sum(is.na(as.numeric(x))) < length(x) * 0.75 & stats::sd(x, na.rm = TRUE) > 0) {
        return(as.numeric(x))
      } # if NAs produced with as.numeric are less than 75%, and SD is greater than 0 we consider the variable numeric
      else if (length(unique(x)) > 1 &
               length(unique(x)) < length(x)) {
        return(as.factor(make.names(x))) # returning syntactically valid names
      } # if the variable is a character, it should have unique values more than 1 and less than total
      else {
        return(rep(NA, length(x)))
      } # if the requirements are not fulfilled, we discard the variable
    }),
    stringsAsFactors = FALSE
    )
  })
  
  # adding Slide as factor
  clean_sample_sheet[["Slide"]] <- as.factor(clean_sample_sheet[["Slide"]])
  
  # adding donorvar as factor
  clean_sample_sheet[[donorvar]] <- as.factor(clean_sample_sheet[[donorvar]]) 
  
  #Removing xMed and yMed variables
  clean_sample_sheet[["xMed"]] <- NULL
  clean_sample_sheet[["yMed"]] <- NULL
  
  clean_sample_sheet <- clean_sample_sheet[, colSums(is.na(clean_sample_sheet)) < nrow(clean_sample_sheet)] # cleaning all NA variables
  
  clean_sample_sheet
  
}

normalize_rgset <- function(rgset, normalization_mode, dropSNPs,
                            maf, dropCpHs, dropSex) {
  # Normalizing data by the selected method
  if (normalization_mode == "Illumina") {
    gset <- minfi::mapToGenome(minfi::ratioConvert(
      type = "Illumina",
      betaThreshold = 0.001,
      minfi::preprocessIllumina(rgset,
        bg.correct = TRUE,
        normalize = "controls"
      )
    ))
  } else if (normalization_mode == "Raw") {
    gset <- minfi::mapToGenome(minfi::ratioConvert(minfi::preprocessRaw(rgset)))
  } else if (normalization_mode == "Noob") {
    gset <- minfi::mapToGenome(minfi::ratioConvert(minfi::preprocessNoob(rgset)))
  } else if (normalization_mode == "Funnorm") {
    gset <- minfi::preprocessFunnorm(rgset)
  } else if (normalization_mode == "SWAN") {
    gset <- minfi::preprocessSWAN(rgset, mSet = minfi::mapToGenome(minfi::preprocessRaw(rgset))) # MethylSet or GenomicMethylSet?
  } else if (normalization_mode == "Quantile") {
    gset <- minfi::preprocessQuantile(rgset)
  } else if (normalization_mode == "Noob+Quantile") {
    gset <- minfi::preprocessQuantile(minfi::preprocessNoob(rgset))
  }

  # prior CpG probes
  probes <- length(minfi::featureNames(gset))

  # remove SNPs to proceed with the analysis and add sex
  # column
  if (dropSNPs) {
    gset <- minfi::dropLociWithSnps(gset, maf = maf)
  }

  # remove CpHs
  if (dropCpHs) {
    gset <- minfi::dropMethylationLoci(gset,
      dropRS = TRUE,
      dropCH = TRUE
    )
  }

  # Add sex info
  gset <- minfi::addSex(gset)

  # remove chromosomes
  if (dropSex) {
    gset <- gset[rownames(minfi::getAnnotation(gset))[!(minfi::getAnnotation(gset)$chr %in%
      c("chrX", "chrY"))], ]
  }

  # Info CpGs removed
  message(paste(
    probes - length(minfi::featureNames(gset)),
    "probes were filtered out."
  ))

  gset
}

# MODEL GENERATION

generate_contrasts <- function(voi) {
  conts <- utils::combn(levels(voi), 2)

  sprintf("%s-%s", conts[1, ], conts[2, ])
}

generate_design <- function(voi, sample_name, covariables, interactions,
                            sample_sheet) {
  voi_factor <- factor(make.names(sample_sheet[[voi]]))

  formula <- stats::as.formula(paste0("~ 0 + ", paste(c(
    voi,
    covariables, interactions
  ), collapse = "+")))

  # Bulding the design matrix
  design <- stats::model.matrix(formula, data = sample_sheet)


  colnames(design)[seq_len(length(unique(voi_factor)))] <- levels(voi_factor)
  row.names(design) <- sample_sheet[[sample_name]]
  colnames(design) <- make.names(colnames(design), unique = TRUE)

  design
}

generate_limma_fit <- function(Mvalues, design, weighting) {
  if (weighting) {
    weights <- limma::arrayWeights(Mvalues, design = design)
  } else {
    weights <- NULL
  }

  limma::lmFit(Mvalues, design, weights = weights)
}

# DMPs CORE FUNCTIONS
calculate_global_difs <- function(Bvalues_totales, grupos, contrasts,
                                  cores) {
  if (!is.null(Bvalues_totales$cpg)) {
    rownames(Bvalues_totales) <- Bvalues_totales$cpg
    Bvalues_totales$cpg <- NULL
  }


  doParallel::registerDoParallel(cores)

  # calcular medias
  all.means <- foreach::foreach(grupo = levels(grupos), .combine = cbind) %dopar% {
    isolate({
      single.means <- data.frame(cpg = rownames(Bvalues_totales))
      nombre_media <- paste("mean", grupo, sep = "_")
      single.means[nombre_media] <- rowMeans(Bvalues_totales[,
        grupos == grupo,
        drop = FALSE
      ])
      single.means[, -1, drop = FALSE]
    })
  }

  all.means$cpg <- row.names(Bvalues_totales)

  # calcular diferencias
  all.dif <- foreach::foreach(cont = contrasts, .combine = cbind) %dopar% {
    isolate({
      single_dif <- data.frame(cpg = rownames(Bvalues_totales))
      grupo1 <- limma::strsplit2(cont, "-")[1]
      grupo2 <- limma::strsplit2(cont, "-")[2]
      nombre_dif <- paste("dif", grupo1, grupo2, sep = "_")
      single_dif[nombre_dif] <- all.means[paste("mean",
        grupo1,
        sep = "_"
      )] - all.means[paste("mean",
        grupo2,
        sep = "_"
      )]
      single_dif[, -1, drop = FALSE]
    })
  }

  all.dif$cpg <- row.names(Bvalues_totales)




  doParallel::stopImplicitCluster()

  cbind(all.dif, all.means)
}


find_dif_cpgs <- function(design, fit, contrasts, trend = FALSE,
                          robust = FALSE, cores) {
  doParallel::registerDoParallel(cores)

  tabla_global <- foreach::foreach(contrast = contrasts) %do% {
    isolate({
      contraste <- limma::makeContrasts(
        contrasts = contrast,
        levels = design
      )
      fitting <- limma::contrasts.fit(fit, contraste)
      fitting <- limma::eBayes(fitting,
        trend = trend,
        robust = robust
      )

      tt_global <- limma::topTable(fitting,
        coef = 1,
        adjust.method = "fdr", number = Inf, sort.by = "none"
      )
      tt_global <- tt_global[, -1]
      tt_global$cpg <- rownames(tt_global)

      rm(contraste, fitting) # Erase not necessary objects

      tt_global
    })
  }


  doParallel::stopImplicitCluster() # stop cluster

  names(tabla_global) <- contrasts

  tabla_global
}



create_filtered_list <- function(limma_list, global_difs, deltaB,
                                 adjp_max, p.value, cores) {
  force(limma_list)
  force(global_difs)
  force(deltaB)
  force(adjp_max)
  force(p.value)

  doParallel::registerDoParallel(cores)

  filtered_list <- foreach::foreach(cont = names(limma_list)) %dopar% {
    isolate({
      dif_target <- paste("dif", limma::strsplit2(
        cont,
        "-"
      )[1], limma::strsplit2(cont, "-")[2], sep = "_")

      tt_global <- limma_list[[cont]]

      tt_global$dif_current <- global_difs[[dif_target]] # indicamos que contraste se aplica

      tt_global <- tt_global[abs(global_difs[[dif_target]]) >
        deltaB & tt_global$adj.P.Val < adjp_max &
        tt_global$P.Value < p.value, ]

      tt_global <- tt_global[stats::complete.cases(tt_global), ]

      data.table::setDT(tt_global)

      doParallel::stopImplicitCluster()

      tt_global
    })
  }


  names(filtered_list) <- names(limma_list)

  filtered_list
}

# DMRs CORE FUNCTIONS
find_dmrs <- function(find_dif_cpgs, minCpGs, platform, voi,
                      regionsTypes, contrasts, bvalues, permutations, ncores) {
  force(find_dif_cpgs)
  force(minCpGs)
  force(platform)
  force(voi)
  force(regionsTypes)
  force(contrasts)
  force(bvalues)
  force(permutations)

  associationTypes <- paste0(regionsTypes, "_association")

  # Limiting ncores
  if (ncores > 2) {
    ncores <- 2
  }

  # filtering only selected contrasts
  find_dif_cpgs <- find_dif_cpgs[contrasts]

  # Creating rankings to use with mCSEA from limma list
  ranking_list <- lapply(find_dif_cpgs, function(x) {
    ranking <- x[["t"]]
    names(ranking) <- x[["cpg"]]
    ranking
  })

  # Generating design matrix to use with mCSEA
  design_matrix <- data.frame(
    row.names = colnames(bvalues),
    group = voi
  )

  foreach::foreach(ranking = ranking_list, .final = function(x) {
    stats::setNames(
      x,
      names(ranking_list)
    )
  }) %do% {
    isolate({
      # Generating mCSEA result
      result <- mCSEA::mCSEATest(ranking,
        minCpGs = minCpGs,
        methData = bvalues, platform = platform, pheno = design_matrix,
        regionsTypes = regionsTypes, nperm = permutations,
        nproc = ncores
      )

      # Processing associations to match with regions
      for (i in seq_along(associationTypes)) {
        result[[associationTypes[i]]] <- result[[associationTypes[i]]][row.names(result[[regionsTypes[i]]])]
      }

      result
    })
  }
}

add_dmrs_globaldifs <- function(mcsea_result, cpg_globaldifs,
                                regionsTypes) {
  force(mcsea_result)
  force(cpg_globaldifs)
  force(regionsTypes)

  cpg_globaldifs <- data.table::as.data.table(cpg_globaldifs)
  data.table::setkeyv(cpg_globaldifs, "cpg")

  associationTypes <- paste0(regionsTypes, "_association")

  dmrs_globaldifs <- foreach::foreach(
    association = associationTypes
  ) %do% {
    as.data.frame(t(vapply(
      mcsea_result[[1]][[association]],
      function(x) {
        colMeans(cpg_globaldifs[list(x), -c("cpg"),
          nomatch = NULL, mult = "all"
        ],
        na.rm = TRUE
        )
      }, numeric(ncol(cpg_globaldifs[, -c("cpg")]))
    )))
  }

  names(dmrs_globaldifs) <- regionsTypes

  for (contrast in names(mcsea_result)) {
    dif_target <- paste("dif", limma::strsplit2(
      contrast,
      "-"
    )[1], limma::strsplit2(contrast, "-")[2], sep = "_")

    for (region in regionsTypes) {
      mcsea_result[[contrast]][[region]] <- mcsea_result[[contrast]][[region]][row.names(mcsea_result[[1]][[region]]), ]
      mcsea_result[[contrast]][[region]][["dif_beta"]] <- dmrs_globaldifs[[region]][[dif_target]]
    }
  }


  mcsea_result
}

filter_dmrs <- function(mcsea_list, fdr, pval, dif_beta, regionsTypes,
                        contrasts) {
  force(fdr)
  force(pval)
  force(dif_beta)
  force(regionsTypes)
  force(mcsea_list)
  force(contrasts)

  # calculating association names from regionsTypes
  associationTypes <- paste0(regionsTypes, "_association")

  foreach::foreach(result = mcsea_list, .final = function(x) {
    stats::setNames(
      x,
      names(mcsea_list)
    )
  }) %do% {
    isolate({
      # Filtering regions by fdr
      for (region in regionsTypes) {
        result[[region]] <- result[[region]][result[[region]]$padj <
          fdr & result[[region]]$pval < pval & abs(result[[region]]$dif_beta) >=
          dif_beta, ]
      }

      # Filtering associations
      for (i in seq_along(associationTypes)) {
        association <- associationTypes[i]
        result[[association]] <- result[[association]][row.names(result[[regionsTypes[i]]])]
      }

      result
    })
  }
}
