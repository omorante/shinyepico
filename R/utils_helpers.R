globalVariables("cont")
globalVariables("contrast")
globalVariables("grupo")
globalVariables("type")
globalVariables("table")


#calculation of global difs table
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
      #nombre_sd = paste("sd",grupo,sep="_")
      single.means[nombre_media] = rowMeans(Bvalues_totales[, grupos == grupo])
      #single.means[nombre_sd] = apply(Bvalues_totales[grupos == grupo], 1, sd) # We are not using standard deviation for now
      single.means[, -1, drop = FALSE]
    })
  }
  
  #all.means = as.data.frame(all.means)   #desconectado por prueba
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
      single_dif[, -1, drop = FALSE]
    })
  }
  
  all.dif$cpg = row.names(Bvalues_totales)
  
  
  
  
  doParallel::stopImplicitCluster()
  
  cbind(all.dif, all.means) #cambiamos dplyr por cbind
  
  
  #all.dif.means[order(all.dif.means$cpg),] #desactivar orden alfabetico
}


find_dif_cpgs = function (groups,
                          design,
                          fit,
                          contrasts,
                          trend = FALSE,
                          robust = FALSE,
                          cores) {
  doParallel::registerDoParallel(cores)
  
  tabla_global = foreach::foreach(contrast = contrasts) %dopar% {
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
      tt_global = tt_global[, -1]
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
                                sd_cort = Inf,
                                cores) {
  #global_difs = global_difs[order(global_difs$cpg),] #ordenamos global_difs para poder comparar
  
  force(limma_list)
  force(global_difs)
  force(deltaB)
  force(adjp_max)
  force(p.value)
  print("create_filtered_list")
  
  doParallel::registerDoParallel(cores)
  
  filtered_list = foreach::foreach(cont = names(limma_list)) %dopar% {
    isolate({
      dif_target = paste("dif",
                         limma::strsplit2(cont, "-")[1],
                         limma::strsplit2(cont, "-")[2],
                         sep = "_")
      #sd_target1 = paste("sd",limma::strsplit2(cont,"-")[1],sep="_")
      #sd_target2 = paste("sd",limma::strsplit2(cont,"-")[2],sep="_")
      tt_global = limma_list[[cont]]
      
      #tt_global = tt_global[order(tt_global$cpg),] #ordenamos
      
      
      tt_global$dif_current = global_difs[[dif_target]] # indicamos que contraste se aplica
      #tt_global = tt_global[abs(global_difs[[dif_target]]) > deltaB & tt_global$adj.P.Val < adjp_max & tt_global$P.Value < p.value & global_difs[[sd_target1]] < sd_cort & global_difs[[sd_target2]] < sd_cort, ]
      tt_global = tt_global[abs(global_difs[[dif_target]]) > deltaB &
                              tt_global$adj.P.Val < adjp_max &
                              tt_global$P.Value < p.value,] #sin filtro de sd
      tt_global = tt_global[stats::complete.cases(tt_global), ]
      
      doParallel::stopImplicitCluster()
      
      tt_global
      #dplyr::left_join(tt_global, global_difs, by="cpg")   #Por ahora, esta función devuelve solo la tabla filtrada procedente de topTable. En el futuro, si hace falta se le puede añadir anotación, global-difs...
    })
  }
  
  
  names(filtered_list) = names(limma_list)
  
  filtered_list
  
}



create_heatmap = function(filtered_data,
                          global_Bvalues,
                          factorgroups,
                          groups2plot,
                          contrasts2plot,
                          Colv = TRUE,
                          clusteralg = "average",
                          distance = "pearson",
                          scale = "row",
                          static = FALSE) {
  #filtered_data = create_filtered_list( findcpgdata, global_difs, deltaB=deltaB, adjp_max = adjp_max, p.value=p.value, sd_cort=sd_cort)
  
  
  filtered_data = filtered_data[contrasts2plot] # filter contrasts2plot
  dif_cpgs = unique(data.table::rbindlist(filtered_data)$cpg)
  
  if (length(dif_cpgs) == 0 ){
    message("No differences found with this criteria. Heatmap can't be plotted")
    return()
  }

  
  join_table = global_Bvalues[dif_cpgs, ]
  join_table$cpg = NULL
  heatdata = as.matrix(join_table)
  heatdata = heatdata[stats::complete.cases(heatdata),]
  
  class(heatdata) = "numeric"
  
  #subsetting heatdata groups2plot
  heatdata = heatdata[, groups2plot]
  
  #order heatdata by groups
  
  sample_order = c()
  for (name in levels(factorgroups)) {
    sample_order = c(sample_order, colnames(heatdata)[factorgroups %in% name])
  }
  
  heatdata = heatdata[, sample_order]
  message(sample_order)
  
  
  buylrd = c("#313695", "#4575B4", "#74ADD1", "#ABD9E9", "#E0F3F8", "#FFFFBF","#FEE090",
            "#FDAE61", "#F46D43", "#D73027", "#A50026"
  )
  colors.martin = grDevices::colorRampPalette(buylrd)(100)
  
  message("create_heatmap")
  
  
  if (static) {
    message("Static plot selected")
    if (distance == "euclidean")
      distfun = stats::dist
    else
      distfun = function(x)
        stats::as.dist(1 - stats::cor(t(x), method = distance))
    
    gplots::heatmap.2 (
      heatdata,
      col = colors.martin,
      Colv = Colv,
      key.xlab = "B values"
      ,
      na.rm = TRUE,
      colsep = 0,
      rowsep = 0.001,
      sepcolor = "white",
      sepwidth = c(0.05, 0.05)
      ,
      colRow = NULL,
      colCol = NULL,
      cexRow = 1,
      cexCol = 1,
      margins = c(15, 2)
      ,
      labCol = NULL,
      srtRow = NULL,
      srtCol = NULL,
      adjRow = c(0, NA),
      adjCol = c(NA, 0)
      ,
      offsetRow = 0.5,
      offsetCol = 0.5,
      key = TRUE,
      keysize = 1,
      key.title = "",
      key.ylab = "count",
      density.info = "none"
      ,
      trace = "none",
      labRow = "",
      scale = scale,
      dendrogram = "both"
      ,
      distfun = distfun,
      hclustfun = function(x)
        stats::hclust(x, method = clusteralg)
    )
  }
  
  else{
    message("interactive plot selected")
    if (distance == "euclidean") {
      distance = stats::dist
    }
    heatmaply::heatmaply(
      heatdata,
      col = colors.martin,
      Colv = Colv,
      key.title = "",
      na.rm = T,
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


create_filtered_beds = function(filtered_data, annotation, cores) {

  doParallel::registerDoParallel(cores)
  
  annotation$cpg = row.names(annotation)
  annotation = data.table::as.data.table(annotation)
  

  #saving hypo and hyper results individually
  
  result = foreach::foreach (
    table = filtered_data,
    .final = function(x)
      stats::setNames(x, paste0(names(filtered_data),"_hypermethylated"))
  ) %dopar% {
    temp = dplyr::left_join(table[table$dif_current < 0,], annotation, by = "cpg")

    data.table::data.table(
      chr = temp$chr,
      start = temp$pos - 1,
      end = temp$pos,
      name = temp$cpg
    )
  }
  
  result = c(result,
                foreach::foreach (
                  table = filtered_data,
                  .final = function(x)
                    stats::setNames(x, paste0(names(filtered_data),"_hypomethylated"))
                ) %dopar% {
                  temp = dplyr::left_join(table[table$dif_current > 0,], annotation, by = "cpg")
                  data.table::data.table(
                    chr = temp$chr,
                    start = temp$pos - 1,
                    end = temp$pos,
                    name = temp$cpg
                  )
                } )
  

  result[["annotation"]] = data.table::data.table(
    chr = annotation$chr,
    start = annotation$pos - 1,
    end = annotation$pos,
    name = annotation$cpg
  )
  
  result
  
}
