globalVariables("cont")
globalVariables("contrast")
globalVariables("grupo")
globalVariables("type")
globalVariables("table")
globalVariables("ranking")

# NORMALIZATION FUNCTIONS

normalize_rgset = function(rgset, normalization_mode){
  
  if (normalization_mode == "Illumina") {
    gset = minfi::mapToGenome(minfi::ratioConvert(
      betaThreshold = 0.001,
      minfi::preprocessIllumina(
        rgset,
        bg.correct = TRUE,
        normalize = "controls"
      )
    ))
  }
  
  else if (normalization_mode == "Raw") {
    gset = minfi::mapToGenome(minfi::ratioConvert(minfi::preprocessRaw(rgset)))
  }
  
  else if (normalization_mode == "Noob") {
    gset = minfi::mapToGenome(minfi::ratioConvert(minfi::preprocessNoob(rgset)))
  }
  
  else if (normalization_mode == "Funnorm") {
    gset = minfi::preprocessFunnorm(rgset)
  }
  
  else if (normalization_mode == "SWAN") {
    gset = minfi::preprocessSWAN(rgset,
                                 mSet = minfi::mapToGenome(minfi::preprocessRaw(rgset))) #MethylSet or GenomicMethylSet?
  }
  
  else if (normalization_mode == "Quantile") {
    gset = minfi::preprocessQuantile(rgset)
  }
  
  else if (normalization_mode == "Noob+Quantile") {
    gset = minfi::preprocessQuantile(minfi::preprocessNoob(rgset))
  }
  
  gset
}


# MODEL GENERATION

generate_limma_fit = function(Mvalues, design, weighting){
  
  if (weighting) {
      weights = limma::arrayWeights(Mvalues, design = design)
  }
  else {
      weights = NULL
  }
  
  limma::lmFit(Mvalues, design, weights = weights)
}

# DMPs CORE FUNCTIONS

calculate_global_difs = function(Bvalues_totales, grupos, contrasts, cores) {
  if (!is.null(Bvalues_totales$cpg)) {
    rownames(Bvalues_totales) = Bvalues_totales$cpg
    Bvalues_totales$cpg = NULL
  }
  
  
  doParallel::registerDoParallel(cores)
  
  #calcular medias
  all.means = foreach::foreach(grupo = levels(grupos), .combine = cbind) %dopar% {
    isolate({
      single.means = data.frame(cpg = rownames(Bvalues_totales))
      nombre_media = paste("mean", grupo, sep = "_")
      single.means[nombre_media] = rowMeans(Bvalues_totales[, grupos == grupo, drop=FALSE])
      single.means[,-1, drop = FALSE]
    })
  }
  
  all.means$cpg  = row.names(Bvalues_totales)
  
  #calcular diferencias
  all.dif = foreach::foreach(cont = contrasts, .combine = cbind) %dopar% {
    isolate({
      single_dif = data.frame(cpg = rownames(Bvalues_totales))
      grupo1 = limma::strsplit2(cont, "-")[1]
      grupo2 = limma::strsplit2(cont, "-")[2]
      nombre_dif = paste("dif", grupo1, grupo2, sep = "_")
      single_dif[nombre_dif] = all.means[paste("mean", grupo1, sep = "_")] - all.means[paste("mean", grupo2, sep =
                                                                                               "_")]
      single_dif[,-1, drop = FALSE]
    })
  }
  
  all.dif$cpg = row.names(Bvalues_totales)
  
  
  
  
  doParallel::stopImplicitCluster()
  
  cbind(all.dif, all.means)
  
}


find_dif_cpgs = function (groups,
                          design,
                          fit,
                          contrasts,
                          trend = FALSE,
                          robust = FALSE,
                          cores) {
  doParallel::registerDoParallel(cores)
  
  tabla_global = foreach::foreach(contrast = contrasts) %do% {
    isolate({
      contraste = limma::makeContrasts(contrasts = contrast, levels = design)
      fitting = limma::contrasts.fit(fit, contraste)
      fitting = limma::eBayes(fitting, trend = trend, robust = robust)
      
      tt_global = limma::topTable(
        fitting,
        coef = 1,
        adjust.method = "fdr",
        number =  Inf,
        sort.by = "none"
      )
      tt_global = tt_global[,-1]
      tt_global$cpg = rownames(tt_global)
      
      rm(contraste, fitting) #Erase not necessary objects
      
      tt_global
      
      
    })
  }
  
  
  doParallel::stopImplicitCluster() #stop cluster
  
  names(tabla_global) = contrasts
  
  tabla_global
  
}



create_filtered_list = function(limma_list,
                                global_difs,
                                deltaB,
                                adjp_max,
                                p.value,
                                cores) {
  force(limma_list)
  force(global_difs)
  force(deltaB)
  force(adjp_max)
  force(p.value)
  
  doParallel::registerDoParallel(cores)
  
  filtered_list = foreach::foreach(cont = names(limma_list)) %dopar% {
    isolate({
      dif_target = paste("dif",
                         limma::strsplit2(cont, "-")[1],
                         limma::strsplit2(cont, "-")[2],
                         sep = "_")
      
      tt_global = limma_list[[cont]]
      
      tt_global$dif_current = global_difs[[dif_target]] # indicamos que contraste se aplica
      
      tt_global = tt_global[abs(global_difs[[dif_target]]) > deltaB &
                              tt_global$adj.P.Val < adjp_max &
                              tt_global$P.Value < p.value, ]
      
      tt_global = tt_global[stats::complete.cases(tt_global),]
      
      data.table::setDT(tt_global)
      
      doParallel::stopImplicitCluster()
      
      tt_global
    })
  }
  
  
  names(filtered_list) = names(limma_list)
  
  filtered_list
  
}

# DMRs CORE FUNCTIONS


find_dmrs = function(find_dif_cpgs,
                     minCpGs,
                     platform,
                     voi,
                     regionsTypes,
                     contrasts,
                     bvalues,
                     permutations,
                     ncores) {
  force(find_dif_cpgs)
  force(minCpGs)
  force(platform)
  force(voi)
  force(regionsTypes)
  force(contrasts)
  force(bvalues)
  force(permutations)
  
  associationTypes = paste0(regionsTypes, "_association")
  
  #Limiting ncores
  if (ncores > 2) {
    ncores = 2
  }
  
  #filtering only selected contrasts
  find_dif_cpgs = find_dif_cpgs[contrasts]
  
  #Creating rankings to use with mCSEA from limma list
  ranking_list = lapply(find_dif_cpgs, function(x) {
    ranking = x[["t"]]
    names(ranking) = x[["cpg"]]
    ranking
  })
  
  #Generating design matrix to use with mCSEA
  design_matrix = data.frame(row.names = colnames(bvalues), group = voi)
  
  foreach::foreach(
    ranking = ranking_list,
    .final = function(x)
      stats::setNames(x, names(ranking_list))
  ) %do% {
    isolate({
      #Specified seed to obtain always the same results
      set.seed(123)
      
      #Generating mCSEA result
      result = mCSEA::mCSEATest(
        ranking,
        minCpGs = minCpGs,
        methData = bvalues,
        platform = platform,
        pheno = design_matrix,
        regionsTypes = regionsTypes,
        nperm = permutations,
        nproc = ncores
      )
      
      #Processing associations to match with regions
      for(i in seq_along(associationTypes)){
        result[[associationTypes[i]]] = result[[associationTypes[i]]][row.names(result[[regionsTypes[i]]])]
      }
      
      result
    })
  }
}


add_dmrs_globaldifs = function(mcsea_result, cpg_globaldifs, regionsTypes, cores){
  
  force(mcsea_result)
  force(cpg_globaldifs)
  force(regionsTypes)
  force(cores)
  
  
  doParallel::registerDoParallel(cores)
  
  cpg_globaldifs = data.table::as.data.table(cpg_globaldifs) #changed from setDT because of bug with bvalues rownames
  data.table::setkeyv(cpg_globaldifs, "cpg")
  
  associationTypes = paste0(regionsTypes, "_association")
  
  dmrs_globaldifs = foreach::foreach(
    cont = names(mcsea_result),
    .final = function(x)
      stats::setNames(x, names(mcsea_result))
  ) %dopar% {
    isolate({
      dif_target = paste("dif",
                         limma::strsplit2(cont, "-")[1],
                         limma::strsplit2(cont, "-")[2],
                         sep = "_")
      
      result = list()
      
      for (i in seq_along(regionsTypes)) {
        result[[regionsTypes[i]]] = vapply(mcsea_result[[cont]][[associationTypes[i]]], function(x) {
          mean(as.numeric(cpg_globaldifs[list(x), nomatch = NULL, mult = "all"][[dif_target]]), na.rm = TRUE)
        }, numeric(1))
        
      }
      
      result

    })
  }
  
  doParallel::stopImplicitCluster()
  
  isolate({
    for (cont in names(mcsea_result)) {
      for (region in regionsTypes)
      {
        mcsea_result[[cont]][[region]]$dif_beta = dmrs_globaldifs[[cont]][[region]]
      }
    }
  })
  
  mcsea_result
}


filter_dmrs = function(mcsea_list,
                       fdr,
                       pval,
                       dif_beta,
                       regionsTypes,
                       contrasts) {
  force(fdr)
  force(pval)
  force(dif_beta)
  force(regionsTypes)
  force(mcsea_list)
  force(contrasts)
  
  #calculating association names from regionsTypes
  associationTypes = paste0(regionsTypes, "_association")
  
  foreach::foreach(
    result = mcsea_list,
    .final = function(x)
      stats::setNames(x, names(mcsea_list))
  ) %do% {
    isolate({
      #Filtering regions by fdr
      for (region in regionsTypes) {
        result[[region]] = result[[region]][result[[region]]$padj < fdr &
                                              result[[region]]$pval < pval &
                                            abs(result[[region]]$dif_beta) >= dif_beta,] 
        
      }
      
      #Filtering associations
      for (i in seq_along(associationTypes)) {
        association = associationTypes[i]
        result[[association]] = result[[association]][row.names(result[[regionsTypes[i]]])]
      }
      
      result
    })
  }
  
}


create_dmrs_heatdata = function(mcsea_result, bvalues, regions, contrasts) {

  mcsea_result = mcsea_result[contrasts]
  associations = paste0(regions, "_association")
  
  
  bvalues$cpg = row.names(bvalues)
  data.table::setDT(bvalues)
  data.table::setkeyv(bvalues, "cpg")
  
  list_heatdata = list()
  identifiers = c()
  
  for (contrast in contrasts) {
    for (i in seq_along(regions)) {
      list_heatdata[[paste0(contrast, "|", regions[i])]] = data.table::as.data.table(t(vapply(mcsea_result[[contrast]][[associations[i]]], function(x) {
        colMeans(bvalues[list(x), -c("cpg"), nomatch = NULL, mult = "all"])
      }, numeric(ncol(bvalues) - 1))))
      identifiers = c(identifiers, paste(contrast, regions[i], names(mcsea_result[[contrast]][[associations[i]]]), sep="|"))
      
    }
  }
  
  heatdata = data.table::rbindlist(list_heatdata)
  colnames(heatdata) = colnames(bvalues)[-ncol(bvalues)]
  heatdata = as.data.frame(heatdata)
  row.names(heatdata) = identifiers
  
  heatdata
}
  







