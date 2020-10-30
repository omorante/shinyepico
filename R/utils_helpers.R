globalVariables("cont")
globalVariables("contrast")
globalVariables("grupo")
globalVariables("type")
globalVariables("table")

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
      setNames(x, names(ranking_list))
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
      setNames(x, names(mcsea_result))
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
      setNames(x, names(mcsea_list))
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
  
  for (contrast in contrasts) {
    for (i in seq_along(regions)) {
      list_heatdata[[paste0(contrast, "|", regions[i])]] = data.table::as.data.table(t(vapply(mcsea_result[[contrast]][[associations[i]]], function(x) {
        colMeans(bvalues[list(x), -c("cpg"), nomatch = NULL, mult = "all"])
      }, numeric(ncol(bvalues) - 1))))
    }
  }
  
  heatdata = data.table::rbindlist(list_heatdata)
  colnames(heatdata) = colnames(bvalues)[-ncol(bvalues)]
  
  heatdata
}
  







# DOWNLOAD HANDLER FUNCTIONS
fwrite_bed = function(bed_file, file_name, DMR = FALSE) {
  if (!DMR) {
    bed_file = data.table::data.table(
      chr = bed_file$chr,
      start = bed_file$pos - 1,
      end = bed_file$pos,
      name = bed_file$cpg
    )
  }
  
  bed_file$start = trimws(format(as.numeric(bed_file$start), scientific = FALSE))
  bed_file$end = trimws(format(as.numeric(bed_file$end), scientific = FALSE))
  
  data.table::fwrite(
    bed_file,
    file = file_name,
    sep = "\t",
    quote = FALSE,
    col.names = FALSE,
    row.names = FALSE
  )
}


create_filtered_beds_dmrs = function(mcsea_filtered,
                                     regionsTypes,
                                     annotation,
                                     directory) {
  annotation$cpg = row.names(annotation)
  annotation = data.table::as.data.table(annotation)
  data.table::setkeyv(annotation, "cpg")
  
  associationTypes = paste0(regionsTypes, "_association")
  
  for (contrast in names(mcsea_filtered)) {
    for (association in associationTypes) {
      temp = data.table::as.data.table(t(vapply(names(mcsea_filtered[[contrast]][[association]]), function(name) {
        target = annotation[list(mcsea_filtered[[contrast]][[association]][[name]]), nomatch = NULL, mult = "all"]
        chr = unique(as.character(target$chr))
        ini = min(as.numeric(target$pos))
        fin = max(as.numeric(target$pos))
        name = paste(contrast, association, name, sep = "_")
        
        if (length(chr) > 1) {
          chr = chr[1]
        }
        
        res = c(
          chr = chr,
          start = ini,
          end = fin,
          name = name
        )
        
        if (length(res) != 4) {
          message(
            paste(
              name,
              "DMR is not correctly annotated. NAs will be introduced in final bed"
            )
          )
          res = c(
            chr = "NA",
            start = "NA",
            end = "NA",
            name = name
          )
        }
        
        res
        
      }, character(4))))
      
      fwrite_bed(
        temp,
        file_name = paste0(directory, "/", contrast, "_", association, ".bed"),
        DMR = TRUE
      )
      
    }
  }
}

create_filtered_bed_dmrs_clusters = function(dendro_data, annotation, directory){
  
  
  
  
}

create_filtered_beds = function(filtered_data, annotation, directory) {
  
  annotation$cpg = row.names(annotation)
  annotation = data.table::as.data.table(annotation)
  
  #saving hypo and hyper results individually
  lapply(names(filtered_data), function(name) {
    temp = data.table::merge.data.table(filtered_data[[name]][filtered_data[[name]]$dif_current < 0,],
                                        annotation,
                                        by = "cpg",
                                        all.x = TRUE)
    
    fwrite_bed(temp,
               file_name = paste0(directory, "/", name, "_hypermethylated.bed"))
  })
  
  lapply(names(filtered_data), function(name) {
    temp = data.table::merge.data.table(filtered_data[[name]][filtered_data[[name]]$dif_current > 0, ],
                                        annotation,
                                        by = "cpg",
                                        all.x = TRUE)
    
    fwrite_bed(temp,
               file_name = paste0(directory, "/", name, "_hypomethylated.bed"))
  })
  
  
  fwrite_bed(annotation, file_name = paste0(directory, "/", "annotation.bed"))
  
}

create_filtered_bed_clusters = function(dendro_data, annotation, directory) {
  annotation$cpg = row.names(annotation)
  annotation = data.table::setDT(annotation)
  
  #saving results by cluster
  lapply(unique(dendro_data), function(cluster) {
    temp = annotation[annotation$cpg %in% names(dendro_data)[dendro_data == cluster], ]
    
    fwrite_bed(temp,
               file_name = paste0(
                 directory,
                 "/",
                 "Cluster_",
                 which(unique(dendro_data) %in% cluster),
                 ".bed"
               ))
  })
  
  fwrite_bed(annotation, file_name = paste0(directory, "/", "annotation.bed"))
}

#GRAPHIC FUNCTIONS
create_heatmap = function(plot_data,
                          factorgroups,
                          groups2plot,
                          Colv = TRUE,
                          ColSideColors = FALSE,
                          RowSideColors = NULL,
                          clusteralg = "average",
                          distance = "pearson",
                          scale = "row",
                          static = TRUE) {
  
  heatdata = as.matrix(plot_data)
  heatdata = heatdata[stats::complete.cases(heatdata), ]
  class(heatdata) = "numeric"
  
  #subsetting heatdata groups2plot
  heatdata = heatdata[, groups2plot]
  
  #optional colside colors
  colside_scale = grDevices::rainbow(length(levels(factorgroups)))
  names(colside_scale) = levels(factorgroups)
  
  #order heatdata by groups
  sample_order = c()
  color_order = c()
  group_order = c()
  for (name in levels(factorgroups)) {
    sample_order = c(sample_order, colnames(heatdata)[factorgroups %in% name])
    color_order = c(color_order, rep(colside_scale[[name]], length(colnames(heatdata)[factorgroups %in% name])))
    group_order = c(group_order, rep(name, length(colnames(heatdata)[factorgroups %in% name])))
  }

  
  heatdata = heatdata[, sample_order]
  
  buylrd = c(
    "#313695",
    "#4575B4",
    "#74ADD1",
    "#ABD9E9",
    "#E0F3F8",
    "#FFFFBF",
    "#FEE090",
    "#FDAE61",
    "#F46D43",
    "#D73027",
    "#A50026"
  )
  colors.martin = grDevices::colorRampPalette(buylrd)(100)

  
  if(!ColSideColors){
    color_order = rlang::missing_arg()
    group_order = rlang::missing_arg()
  }
  
  if(is.null(RowSideColors))
    row_order = rlang::missing_arg()
  else
    row_order = RowSideColors

  if (static) {
    if (distance == "euclidean")
      distfun = stats::dist
    else
      distfun = function(x)
        stats::as.dist(1 - stats::cor(t(x), method = distance))
    
    gplots::heatmap.2 (
      heatdata,
      col = colors.martin,
      Colv = Colv,
      ColSideColors = color_order,
      RowSideColors = row_order,
      key.xlab = "B values",
      na.rm = TRUE,
      colsep = 0,
      rowsep = 0.001,
      sepcolor = "white",
      sepwidth = c(0.05, 0.05),
      colRow = NULL,
      colCol = NULL,
      cexRow = 1,
      cexCol = 1,
      margins = c(12, 1),
      labCol = NULL,
      srtRow = NULL,
      srtCol = NULL,
      adjRow = c(0, NA),
      adjCol = c(NA, 0),
      offsetRow = 0.5,
      offsetCol = 0.5,
      key = TRUE,
      keysize = 1,
      key.title = "",
      key.ylab = "count",
      density.info = "none",
      trace = "none",
      labRow = "",
      scale = scale,
      dendrogram = "both",
      distfun = distfun,
      hclustfun = function(x)
        stats::hclust(x, method = clusteralg)
    )
  }
  
  else{
    
    if (distance == "euclidean") {
      distance = stats::dist
    }
    heatmaply::heatmaply(
      heatdata,
      col = colors.martin,
      Colv = Colv,
      col_side_colors = group_order,
      row_side_colors = row_order,
      key.title = "",
      na.rm = TRUE,
      dendogram = "both",
      scale = scale,
      distfun = distance,
      hclustfun = function(x)
      stats::hclust(x, method = clusteralg),
      seriate = "mean",
      row_dend_left = TRUE,
      showticklabels = c(TRUE, FALSE),
      branches_lwd = 0.3,
      plot_method = "plotly",
      colorbar_xpos = -0.01,
      colorbar_ypos = 0.3,
      margins = c(25, 25, NA, 0)
    )
    
  }
  
}


create_pca = function(Bvalues,
                      pheno_info,
                      pc_x = "PC1",
                      pc_y = "PC2",
                      group,
                      color = NULL) {
  pca_data = stats::prcomp(t(stats::na.omit(Bvalues)))
  pca_info = as.data.frame(summary(pca_data)[["importance"]])
  
  pca_data = as.data.frame(pca_data$x)
  pca_data = cbind(pca_data, pheno_info)
  
  
  pca_graph = plotly::ggplotly(
    ggplot2::ggplot(pca_data) +
      ggplot2::geom_point(
        ggplot2::aes_string(
          x = pc_x,
          y = pc_y,
          group = group,
          color = color
        ),
        size = 3
      ) +
      ggplot2::theme_bw()
  ) %>%
    plotly_config()
  
  return(list(info = pca_info, graph = pca_graph))
  
  
}


create_corrplot = function(Bvalues, clean_sample_sheet, sample_target_sheet) {
  #bvalues should be a dataframe
  
  pca_data = as.data.frame(stats::prcomp(t(stats::na.omit(Bvalues)))$x)
  cor_data = as.data.frame(cor3(pca_data, clean_sample_sheet))
  cor_data$Var1 = row.names(cor_data)
  
  data.table::setDT(cor_data)
  
  cor_data = data.table::melt(
    cor_data,
    id.vars = "Var1",
    variable.name = "Var2",
    value.name = "cor"
  )
  
  cor_data$Var1 = factor(cor_data$Var1, levels = colnames(pca_data))
  
  corr_graph = plotly::ggplotly(
    ggplot2::ggplot(cor_data, ggplot2::aes_string("Var1", "Var2", fill = "cor")) + ggplot2::geom_tile(color =
                                                                                                        "darkgrey", size = 1) +
      ggplot2::scale_fill_gradient2(
        high = "#6D9EC1",
        low = "#E46726",
        mid = "white",
        midpoint = 0,
        limit = c(-1, 1),
        space = "Lab",
        name = "Correlation"
      ) +
      ggplot2::theme_bw() +
      ggplot2::theme(
        panel.grid = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_text(
          angle = 45,
          vjust = 1,
          size = 12,
          hjust = 1
        )
      ) +
      ggplot2::labs(x = "", y = "")
  ) %>% plotly_config(fixedrange = FALSE)
  
  
  corr_info = data.frame(
    Variable = colnames(clean_sample_sheet),
    Type = unlist(lapply(clean_sample_sheet, class))
  )
  
  corr_info$Correlation = ifelse(corr_info$Type == "numeric", "Pearson", "R-squared")
  
  if (length(colnames(sample_target_sheet)[!(colnames(sample_target_sheet) %in% colnames(clean_sample_sheet))]) > 0)
  {
    corr_disc = data.frame(
      Variable = colnames(sample_target_sheet)[!(colnames(sample_target_sheet) %in% colnames(clean_sample_sheet))],
      Type = "Discarded",
      Correlation = "None"
    )
  }
  else{
    corr_disc = data.frame()
  }
  
  corr_info = rbind(corr_info, corr_disc)
  
  return(list(info = corr_info, graph = corr_graph))
  
  }


create_plotqc = function(rgset, sample_names, badSampleCutoff = 10) {
  plotly::ggplotly(
    as.data.frame(minfi::getQC(minfi::preprocessRaw(rgset))) %>%
      dplyr::mutate(
        Sample = sample_names,
        Quality = ifelse(
          .data$mMed < badSampleCutoff |
            .data$uMed < badSampleCutoff,
          "Suboptimal",
          "OK"
        )
      ) %>%
      ggplot2::ggplot(
        ggplot2::aes_string(
          x = "mMed",
          y = "uMed",
          group = "Sample",
          color = "Quality"
        )
      ) +
      ggplot2::geom_point() +
      ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = 1.20)) +
      ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = 1.20)) +
      ggplot2::scale_color_manual(values = c("black", "darkred")) +
      ggplot2::geom_abline(slope = 1, intercept = 0) +
      ggplot2::geom_abline(
        slope = 1,
        intercept = 0.5,
        linetype = "dotted",
        color = "darkgrey"
      ) +
      ggplot2::geom_abline(
        slope = 1,
        intercept = -0.5,
        linetype = "dotted",
        color = "darkgrey"
      ) +
      ggplot2::theme_bw() +
      ggplot2::theme(
        panel.grid.major = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank()
      )
  ) %>%
    plotly_config()
  
  
}

create_densityplot = function(Bvalues, n = 200000) {
  #Creating density plot using a sample of n CpGs
  plotly::ggplotly(
    Bvalues[sample(seq_len(nrow(Bvalues)), n), ] %>%
      tidyr::pivot_longer(
        cols = seq_len(ncol(Bvalues)),
        names_to = "sample",
        values_to = "Bvalues"
      ) %>%
      ggplot2::ggplot(ggplot2::aes(
        x = .data$Bvalues, color = .data$sample
      )) + ggplot2::stat_density(position = "identity", geom = "line") +
      ggplot2::theme_bw() + ggplot2::theme(
        panel.grid.major = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank()
      )
  ) %>%
    plotly_config() #%>% plotly::toWebGL()
  
}

create_boxplot = function(Bvalues, n = 200000) {
  Bvalues[sample(seq_len(nrow(Bvalues)), n), ] %>%
    tidyr::pivot_longer(cols = seq_len(ncol(Bvalues)),
                        names_to = "sample",
                        values_to = "Bvalues") %>%
    ggplot2::ggplot(ggplot2::aes(x = .data$sample, y = .data$Bvalues)) +
    ggplot2::geom_boxplot(
      outlier.colour = "black",
      outlier.shape = 16,
      outlier.size = 2,
      notch = FALSE,
      fill = "#56B1F7"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank()
    ) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = -45, hjust =
                                                         0))
  
}

create_dendrogram = function(plot_data,
                             factorgroups,
                             groups2plot,
                             clusteralg = "average",
                             distance = "pearson",
                             scale_selection = "row",
                             k_number) {
  heatdata = as.matrix(plot_data)
  heatdata = heatdata[stats::complete.cases(heatdata),]
  class(heatdata) = "numeric"
  
  #subsetting heatdata groups2plot
  heatdata = heatdata[, groups2plot]
  
  #scaling data if option is selected
  if (scale_selection == "row")
    heatdata = t(scale(t(heatdata)))
  
  
  #order heatdata by groups
  sample_order = c()
  for (name in levels(factorgroups)) {
    sample_order = c(sample_order, colnames(heatdata)[factorgroups %in% name])
  }
  
  heatdata = heatdata[, sample_order]
  
  if (distance == "euclidean")
    distfun = stats::dist
  else
    distfun = function(x)
      stats::as.dist(1 - stats::cor(t(x), method = distance))
  
  
  distance = distfun(heatdata)
  clustering = stats::hclust(distance, clusteralg)
  clusters = stats::cutree(clustering, k = k_number)
  
  rowside_scale = grDevices::rainbow(length(unique(clusters)))
  
  for (number in unique(clusters)) {
    clusters[clusters == number] = rowside_scale[number]
  }
  
  clusters
}


create_snpheatmap = function(snps, sample_names, color_groups) {
  buylrd = c(
    "#313695",
    "#4575B4",
    "#74ADD1",
    "#ABD9E9",
    "#E0F3F8",
    "#FFFFBF",
    "#FEE090",
    "#FDAE61",
    "#F46D43",
    "#D73027",
    "#A50026"
  )
  
  
  colors.martin = grDevices::colorRampPalette(buylrd)(100)
  
  colnames(snps) = sample_names
  
  heatmaply::heatmaply(
    snps,
    col = colors.martin,
    Colv = TRUE,
    key.title = "",
    na.rm = TRUE,
    dendogram = "both",
    scale = "row",
    col_side_colors = color_groups,
    distfun = "pearson",
    hclustfun = function(x)
      stats::hclust(x, method = "average"),
    seriate = "mean",
    row_dend_left = TRUE,
    showticklabels = c(TRUE, FALSE),
    branches_lwd = 0.3,
    plot_method = "plotly",
    colorbar_xpos = -0.01,
    colorbar_ypos = 0.3,
    margins = c(25, 25, NA, 0),
  ) %>% plotly_config()
}

create_bisulfiteplot = function(rgset, sample_names, threshold = 1.5) {
  colnames(rgset) = sample_names
  
  ctrlAddress = minfi::getControlAddress(rgset, controlType = "BISULFITE CONVERSION II")
  ctlWide = as.matrix(minfi::getRed(rgset)[ctrlAddress, , drop = FALSE])
  ctlR = reshape2::melt(ctlWide, varnames = c("address", "sample"))
  ctlWide = as.matrix(minfi::getGreen(rgset)[ctrlAddress, , drop = FALSE])
  ctlG = reshape2::melt(ctlWide, varnames = c("address", "sample"))
  ctl = rbind(cbind(channel = "Red", ctlR), cbind(channel = "Green",
                                                  ctlG))
  
  plot_data = tidyr::pivot_wider(
    ctl,
    id_cols = 2:4,
    values_from = "value",
    names_from = "channel"
  )
  plot_data$Ratio = plot_data$Red / plot_data$Green
  plot_data = stats::aggregate(stats::as.formula("Ratio ~ sample"),
                               data = plot_data,
                               FUN = min)
  plot_data$Status = ifelse(plot_data$Ratio > threshold, "OK", "Failed")
  
  plotly::ggplotly(
    ggplot2::ggplot(plot_data) +
      ggplot2::geom_segment(
        ggplot2::aes_string(
          y = 0,
          x = "sample",
          yend = "Ratio",
          xend = "sample"
        ),
        color = "darkgrey"
      ) +
      ggplot2::geom_point(
        ggplot2::aes_string(x = "sample", y = "Ratio", fill = "Status"),
        size = 4
      ) +
      ggplot2::theme_bw() +
      ggplot2::scale_fill_manual(values = c("#2a9d8f", "#e76f51")) +
      ggplot2::geom_hline(yintercept = 1.5, linetype = "dashed") +
      ggplot2::xlab("") +
      ggplot2::ylab("Minimum Ratio Converted/Non-Converted") +
      ggplot2::coord_flip()
  ) %>% plotly_config()
  
  
  
}

create_sexplot = function(gset, sample_names) {
  
  sex_info = as.data.frame(minfi::pData(gset)[,c("xMed","yMed","predictedSex")])
  sex_info$sample = sample_names
  
  plotly::ggplotly(
    ggplot2::ggplot(
      sex_info,
      ggplot2::aes_string(x = "xMed", y = "yMed", color = "predictedSex")
    ) +
      ggplot2::geom_point(size = 3) +
      ggplot2::geom_abline(
        intercept = -2,
        slope = 1,
        color = "darkgrey",
        linetype = "dashed"
      ) +
      ggplot2::theme_bw()
  ) %>% plotly_config()
  
}

create_plotSA = function(fit) {
  plot_data = data.frame(Amean = fit$Amean, sigma = sqrt(fit$sigma))
  
  ggplot2::ggplot(plot_data, ggplot2::aes_string(x = "Amean", y = "sigma")) +
    ggplot2::geom_bin2d(bins = 500) +
    ggplot2::labs(x = "Amean", y = "sqrt(sigma)") +
    ggplot2::theme_bw()
}

plotly_config = function(plotly_object, fixedrange = TRUE) {
  plotly_object %>%
    plotly::config(
      displaylogo = FALSE,
      modeBarButtonsToRemove = c("zoomIn2d", "zoomOut2d"),
      toImageButtonOptions = list(format = "svg")
    ) %>%
    plotly::layout(xaxis = list(fixedrange = fixedrange)) %>%
    plotly::layout(yaxis = list(fixedrange = fixedrange))
}


#PLOT AUX FUNCTIONS

cor3 = function(df1, df2) {
  #function based on cor2 function of https://gist.github.com/talegari
  #This function handle the comparison fo numeric variables (PCs) 
  #with numeric or factor variables.
  
  stopifnot(inherits(df1, "data.frame"))
  stopifnot(
    vapply(df1, class, FUN.VALUE = character(1)) %in%
      c("integer"
        , "numeric"
        , "factor"
        , "character")
  )
  
  stopifnot(inherits(df2, "data.frame"))
  stopifnot(
    vapply(df2, class, FUN.VALUE = character(1)) %in%
      c("integer"
        , "numeric"
        , "factor"
        , "character")
  )
  
  cor_fun <- function(pos_df1, pos_df2) {
    # both are numeric
    if (class(df1[[pos_df1]]) %in% c("integer", "numeric") &&
        class(df2[[pos_df2]]) %in% c("integer", "numeric")) {
      r <- stats::cor(df1[[pos_df1]]
                      , df2[[pos_df2]]
                      , use = "pairwise.complete.obs")
    }
    
    # one is numeric and other is a factor/character
    if (class(df1[[pos_df1]]) %in% c("integer", "numeric") &&
        class(df2[[pos_df2]]) %in% c("factor", "character")) {
      r <- sqrt(summary(stats::lm(df1[[pos_df1]] ~ as.factor(df2[[pos_df2]])))[["r.squared"]])
    }
    
    if (class(df2[[pos_df2]]) %in% c("integer", "numeric") &&
        class(df1[[pos_df1]]) %in% c("factor", "character")) {
      r <- sqrt(summary(stats::lm(df2[[pos_df2]] ~ as.factor(df1[[pos_df1]])))[["r.squared"]])
    }
    
    return(r)
  }
  
  cor_fun <- Vectorize(cor_fun)
  
  # now compute corr matrix
  corrmat <- outer(seq_len(ncol(df1))
                   , seq_len(ncol(df2))
                   , function(x, y)
                     cor_fun(x, y))
  
  rownames(corrmat) <- colnames(df1)
  colnames(corrmat) <- colnames(df2)
  
  return(corrmat)
}