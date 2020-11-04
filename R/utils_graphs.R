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


create_corrplot = function(Bvalues, clean_sample_sheet, sample_target_sheet, p.value) {
  #bvalues should be a dataframe
  
  pca_data = as.data.frame(stats::prcomp(t(stats::na.omit(Bvalues)))$x)
  cor_data = as.data.frame(cor3(pca_data, clean_sample_sheet, p.value = p.value))
  cor_data$Var1 = row.names(cor_data)
  
  data.table::setDT(cor_data)
  
  cor_data = data.table::melt(
    cor_data,
    id.vars = "Var1",
    variable.name = "Var2",
    value.name = "value"
  )
  
  cor_data$Var1 = factor(cor_data$Var1, levels = colnames(pca_data))
  
  if(p.value){
    limit = c(0,1)
    name = "p.value"
    low = "#6D9EC1"
    high = "white"
  }
  else{
    limit = c(-1,1)
    name = "Correlation"
    high = "#6D9EC1"
    low = "#E46726"
    mid = "white"
  }
  
  corr_graph = plotly::ggplotly(
    ggplot2::ggplot(cor_data, ggplot2::aes_string("Var1", "Var2", fill = "value")) + ggplot2::geom_tile(color =
                                                                                                          "darkgrey", size = 1) +
      {if(p.value){
        ggplot2::scale_fill_gradient(
          high = high,
          low = low,
          limit = limit,
          space = "Lab",
          name = name
        )
      }
        else{
          ggplot2::scale_fill_gradient2(
            high = high,
            low = low,
            mid = mid,
            midpoint = 0,
            limit = limit,
            space = "Lab",
            name = name
          )
        }} +
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

cor3 = function(df1, df2, p.value = FALSE) {
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
  
  cor_fun = function(pos_df1, pos_df2) {
    # both are numeric
    if (class(df1[[pos_df1]]) %in% c("integer", "numeric") &&
        class(df2[[pos_df2]]) %in% c("integer", "numeric")) {
      r = stats::cor(df1[[pos_df1]]
                     , df2[[pos_df2]]
                     , use = "pairwise.complete.obs")
    }
    
    # one is numeric and other is a factor/character
    if (class(df1[[pos_df1]]) %in% c("integer", "numeric") &&
        class(df2[[pos_df2]]) %in% c("factor", "character")) {
      r = sqrt(summary(stats::lm(df1[[pos_df1]] ~ as.factor(df2[[pos_df2]])))[["r.squared"]])
    }
    
    if (class(df2[[pos_df2]]) %in% c("integer", "numeric") &&
        class(df1[[pos_df1]]) %in% c("factor", "character")) {
      r = sqrt(summary(stats::lm(df2[[pos_df2]] ~ as.factor(df1[[pos_df1]])))[["r.squared"]])
    }
    
    return(r)
  }
  
  cor_fun_pvalue = function(pos_df1, pos_df2) {
    # both are numeric
    if (class(df1[[pos_df1]]) %in% c("integer", "numeric") &&
        class(df2[[pos_df2]]) %in% c("integer", "numeric")) {
      r = stats::cor.test(df1[[pos_df1]]
                          , df2[[pos_df2]]
                          , use = "pairwise.complete.obs")[["p.value"]]
    }
    
    # one is numeric and other is a factor/character
    if (class(df1[[pos_df1]]) %in% c("integer", "numeric") &&
        class(df2[[pos_df2]]) %in% c("factor", "character")) {
      r = sqrt(summary(stats::lm(df1[[pos_df1]] ~ as.factor(df2[[pos_df2]])))[["coefficients"]][2,4])
    }
    
    if (class(df2[[pos_df2]]) %in% c("integer", "numeric") &&
        class(df1[[pos_df1]]) %in% c("factor", "character")) {
      r = sqrt(summary(stats::lm(df2[[pos_df2]] ~ as.factor(df1[[pos_df1]])))[["coefficients"]][2,4])
    }
    
    return(r)
  }
  
  cor_fun = Vectorize(cor_fun)
  cor_fun_pvalue = Vectorize(cor_fun_pvalue)
  
  # now compute corr matrix
  
  if(p.value){
    corrmat = outer(seq_len(ncol(df1))
                    , seq_len(ncol(df2))
                    , function(x, y)
                      cor_fun_pvalue(x, y))
  }
  else{
    corrmat = outer(seq_len(ncol(df1))
                    , seq_len(ncol(df2))
                    , function(x, y)
                      cor_fun(x, y))
  }
  
  rownames(corrmat) = colnames(df1)
  colnames(corrmat) = colnames(df2)
  
  return(corrmat)
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

create_individual_boxplot = function(Bvalues, cpg, voi) {
  Bvalues[cpg,,drop=FALSE] %>%
    tidyr::pivot_longer(cols = seq_len(ncol(Bvalues)),
                        names_to = "sample",
                        values_to = "Bvalues") %>%
    dplyr::mutate(group = voi) %>%
    ggplot2::ggplot(ggplot2::aes_string(x = "group", y = "Bvalues")) +
    ggplot2::geom_boxplot(
      outlier.colour = "black",
      outlier.shape = 16,
      outlier.size = 2,
      notch = FALSE,
      fill = "#56B1F7"
    ) +
    ggplot2::ggtitle(cpg)+
    ggplot2::theme_bw() +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank()
    ) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = -45, hjust =
                                                         0),
                   plot.title = ggplot2::element_text(hjust = 0.5))
  
}