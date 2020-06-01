#' The application server-side
#'
#' @param input,output,session Internal parameters for {shiny}.
#'     DO NOT REMOVE.
#' @import shiny
#' @importFrom rlang .data
#'
#' @noRd



app_server = function(input, output, session) {

  #PERFORMANCE SETTINGS
  n_cores = golem::get_golem_options("n_cores")
  options(shiny.maxRequestSize = golem::get_golem_options("max_upload_size"))
  
  #INPUT
  
  #Load button only shows if file is uploaded
  output$ui_button_input_load = renderUI({
    if (!is.null(input$fileinput_input$datapath))
      return(actionButton("button_input_load", "Load Data"))
    else
      return()
  })
  
  #When you press button_input_load, the data is unzipped and the metharray sheet is loaded
  rval_sheet_target = eventReactive(input$button_input_load, 
          rval_sheet()[rval_sheet()[, input$select_input_samplenamevar]  %in% input$selected_samples, ])
  
  rval_sheet = eventReactive(input$button_input_load, {
    utils::unzip(input$fileinput_input$datapath,
                 exdir = paste0(tempdir(), "/experiment_data"))
    minfi::read.metharray.sheet(paste0(tempdir(), "/experiment_data"))
  })
  
  #When you press button_input_load, the form options are updated
  observeEvent(input$button_input_load, {
    updateSelectInput(
      session,
      "select_input_samplenamevar",
      label = "Select Sample Names Column:",
      choices = colnames(rval_sheet())
    )
    updateSelectInput(
      session,
      "select_input_groupingvar",
      label = "Select Grouping Variable:",
      choices = colnames(rval_sheet())
    )
    updateSelectInput(
      session,
      "select_input_donorvar",
      label = "Select Donor Variable:",
      choices = colnames(rval_sheet())
    )
  })
  
  #the checkbox of samples to process is updated when samplenamevar changes
  observeEvent(
    input$select_input_samplenamevar,
    updateCheckboxGroupInput(
      session,
      "selected_samples",
      label = "Select Samples to Process:",
      choices = rval_sheet()[, input$select_input_samplenamevar],
      selected = rval_sheet()[, input$select_input_samplenamevar]
    )
  )
  
  #The dataframe is rendered
  output$samples_table =  DT::renderDT(rval_sheet())
  
  #rval_rgset loads RGSet using read.metharray.exp and the sample sheet (rval_sheet())
  rval_rgset = eventReactive(input$button_input_next, {
    targets = rval_sheet()[rval_sheet()[, input$select_input_samplenamevar]  %in% input$selected_samples, ]
    RGSet = minfi::read.metharray.exp(targets = targets, verbose = TRUE)
    RGSet [(rowMeans(as.matrix(minfi::detectionP(RGSet))) < 0.01), ]
  })
  
  #We change the page to the next one
  observeEvent(input$button_input_next,
      { withProgress(message = "Loading data...", value = 2, max = 5,
                              {rval_rgset()
                                updateNavbarPage(session, "navbar_epic", "Minfi Norm.")
                              })
      })
  
  
  #MINFI NORMALIZATION
  
  
  #Calculation of minfi normalized data
  rval_gset = reactive({
    
    if (input$select_minfi_norm == "Illumina") {
      gset = minfi::mapToGenome(minfi::ratioConvert( 
        type="Illumina",
        minfi::preprocessIllumina(
          rval_rgset(),
          bg.correct = TRUE,
          normalize = "controls"
        )
      ))
    }
    
    else if (input$select_minfi_norm == "Raw") {
      gset = minfi::mapToGenome(minfi::ratioConvert(minfi::preprocessRaw(rval_rgset())))
    }
    
    else if (input$select_minfi_norm == "Noob") {
      gset = minfi::mapToGenome(minfi::ratioConvert(minfi::preprocessNoob(rval_rgset())))
    }
    
    else if (input$select_minfi_norm == "Funnorm") {
      gset = minfi::preprocessFunnorm((rval_rgset()))
    }
    else if (input$select_minfi_norm == "Quantile") {
      gset = minfi::preprocessQuantile(rval_rgset())
    }
    
    else if (input$select_minfi_norm == "Noob+Quantile") {
      gset = minfi::preprocessQuantile(minfi::preprocessNoob(rval_rgset()))
    }
    
    #remove SNPs to proceed with the analysis
    minfi::addSex(minfi::dropLociWithSnps(gset, maf = 0))
    
  })
  
  
  #getBeta/getM reactives
  rval_rgset_getBeta = reactive({
    bvalues = as.data.frame(minfi::getBeta(rval_rgset()))
    colnames(bvalues) = input$selected_samples
    bvalues
  })
  
  rval_gset_getBeta = reactive({
    bvalues = as.data.frame(minfi::getBeta(rval_gset()))
    colnames(bvalues) = input$selected_samples
    bvalues
  })
  
  
  rval_gset_getM = reactive({
    minfi::getM(rval_gset())
  })
  ###
  
  
  
  #Plotting of QC graphs
  output$graph_minfi_qcraw = renderCachedPlot(minfi::plotQC(minfi::getQC(minfi::preprocessRaw(rval_rgset(
  )))), paste0("Raw", "QC"))
  output$graph_minfi_densityplotraw = renderCachedPlot(minfi::densityPlot(rval_rgset()),
                                                       paste0("Raw", "densityplot"))
  output$graph_minfi_densitybeanplotraw = renderCachedPlot(minfi::densityBeanPlot(rval_rgset()),
                                                           paste0("Raw", "densitybeanplot"))
  output$graph_minfi_mdsplotraw = renderCachedPlot(
    minfi::mdsPlot(rval_rgset(), 
                   sampGroups = rval_sheet_target()[,input$select_input_groupingvar],
                   sampNames = rval_sheet_target()[,input$select_input_samplenamevar]),
    paste0("Raw", "mdsplot")
   )
  
  output$graph_minfi_boxplotraw = renderCachedPlot(graphics::boxplot(as.matrix(rval_rgset_getBeta())),
                                                   paste0(input$select_minfi_norm, "boxplot"))
  
  
  output$graph_minfi_densityplot = renderCachedPlot(minfi::densityPlot(as.matrix(rval_gset_getBeta())),
                                                    paste0(input$select_minfi_norm, "densityplot"))
  output$graph_minfi_densitybeanplot = renderCachedPlot(
    minfi::densityBeanPlot(as.matrix(rval_gset_getBeta())),
    paste0(input$select_minfi_norm, "densitybeanplot")
  )
  
  output$graph_minfi_mdsplot = renderCachedPlot(
    minfi::mdsPlot(
      as.matrix(rval_gset_getBeta()),
      sampGroups = rval_sheet_target()[,input$select_input_groupingvar],
      sampNames = rval_sheet_target()[, input$select_input_samplenamevar]
    ),
    paste0(input$select_minfi_norm, "mdsplot")
  )
  output$graph_minfi_boxplot = renderCachedPlot(graphics::boxplot(as.matrix(rval_gset_getBeta())),
                                                paste0(input$select_minfi_norm, "boxplot"))
  
  output$graph_minfi_sex = renderCachedPlot(minfi::plotSex(rval_gset(), 
                                            id = rval_sheet_target()[, input$select_input_samplenamevar] ),
                                            paste0(input$select_minfi_norm, "sex"))
  
  #SNPs heatmap
  
  rval_plot_minfi_snps = reactive({
    buylrd = c("#313695","#4575B4", "#74ADD1", "#ABD9E9", "#E0F3F8",
               "#FFFFBF", "#FEE090","#FDAE61","#F46D43","#D73027","#A50026")
    

    colors.martin = grDevices::colorRampPalette(buylrd)(100)
    snps = minfi::getSnpBeta(rval_rgset())
    rval_sheet_target =  rval_sheet()[rval_sheet()[, input$select_input_samplenamevar]  %in% input$selected_samples, ]
    
    colnames(snps) = rval_sheet_target[, input$select_input_samplenamevar]
    heatmaply::heatmaply(
      snps,
      col = colors.martin,
      Colv = T,
      key.title = "",
      na.rm = T,
      dendogram = "both",
      scale = "row",
      col_side_colors = rval_sheet_target[, input$select_input_donorvar],
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
      margins = c(25, 25, NA, 0)
    )
    
  })
  
  
  output$graph_minfi_snps = plotly::renderPlotly(rval_plot_minfi_snps())
  
  
  
  
  
  
  #Update of next form and move to Limma
  observeEvent(input$button_minfi_select,
               {
                 updateSelectInput(
                   session,
                   "select_limma_voi",
                   label = "Select Variable of Interest",
                   choices = colnames(minfi::pData(rval_gset())),
                   selected = input$select_input_groupingvar
                 )
                 updateCheckboxGroupInput(
                   session,
                   "checkbox_limma_covariables",
                   label = "Select Covariables to block",
                   choices = colnames(minfi::pData(rval_gset())),
                   selected = input$select_input_donorvar
                 )
                 updateCheckboxGroupInput(
                   session,
                   "checkbox_limma_groups",
                   label = "Select Groups to compare",
                   choices = unique(minfi::pData(rval_gset())[, input$select_input_groupingvar])
                 )
                 updateNavbarPage(session, "navbar_epic", "Limma")
                 
                 updateSelectInput(session,
                                   "select_limma_trend",
                                   label = "eBayes trend option:",
                                   c("TRUE", "FALSE"),
                                   "FALSE")
                 updateSelectInput(session,
                                   "select_limma_robust",
                                   label = "eBayes robust option:",
                                   c("TRUE", "FALSE"),
                                   "FALSE")
               })
  
  
  
  
  #LIMMA
  
  #Variable of interest
  rval_voi = reactive(factor(sub("-","_",minfi::pData(rval_gset())[, input$select_limma_voi]))) #add substitution of "-" for "_", avoiding conflicts

  #Design calculation
  rval_design = reactive({
    #design = model.matrix(~ 0 + rval_voi())
    pdata = as.data.frame(apply(minfi::pData(rval_gset()), 2, sub, pattern = "-", replacement = "_")) #avoiding any "-" in the data
    formula = stats::as.formula(paste0("~ 0 + ", paste(
      c(
        input$select_limma_voi,
        input$checkbox_limma_covariables
      ),
      collapse = "+"
    )))
    print(paste0("~ 0 + ", paste(
      c(
        input$select_limma_voi,
        input$checkbox_limma_covariables
      ),
      collapse = "+"
    )))
    design = stats::model.matrix(formula, data = pdata)
    print(design)
    colnames(design)[1:length(unique(rval_voi()))] = levels(rval_voi())
    design
  })
  
  #Calculation of contrasts
  rval_contrasts = reactive({
    contrastes = c() #calculamos lista de contrastes de todos contra todos
    valor = 1
    all.groups = levels(rval_voi())
    for (i in 1:(length(all.groups) - 1)) {
      for (z in 1:(length(all.groups) - valor)) {
        contrastes = c(contrastes, paste(all.groups[i], all.groups[i + z], sep = "-"))
      }
      valor <- valor + 1
    }
    contrastes
  })
  
  
  #Calculation of limma model
  rval_fit = eventReactive(input$button_limma_calculatemodel, {
    
    if (as.logical(input$select_limma_weights)){
      weights = limma::arrayWeights(rval_gset_getM())
    }
    else { weights = NULL}
    
    limma::lmFit(rval_gset_getM(), rval_design(), weights = weights)
  })
  
  
  #f rval_fit() has NAs, we remove the option to trend or robust in eBayes to prevent failure
  observeEvent(input$button_limma_calculatemodel, {
    if (any(is.na(rval_fit()$coefficients))) {
      updateSelectInput(session, "select_limma_trend", label = "Norm. method is not compatible with trend", "FALSE", "FALSE")
      updateSelectInput(session, "select_limma_robust", label = "Norm. method is not compatible with robust", "FALSE", "FALSE")
    }
  })
  
  
  #Calculation of global difs
  rval_globaldifs = eventReactive(input$button_limma_calculatedifs, {
    calculate_global_difs(rval_gset_getBeta(), rval_voi(), rval_contrasts(), cores=n_cores)
  })
  
  #Calculation of differences (eBayes)
  rval_finddifcpgs = eventReactive(input$button_limma_calculatedifs, {
    find_dif_cpgs(
      rval_voi(),
      rval_design(),
      rval_fit(),
      rval_contrasts(),
      trend = as.logical(input$select_limma_trend),
      robust = as.logical(input$select_limma_robust),
      cores = n_cores
    )
  })
  
  #Update of heatmap controls
  observeEvent(input$button_limma_calculatedifs, {
    updateTabsetPanel(session, "tabset_limma", "differential_cpgs")
    updateSelectInput(
      session,
      "select_limma_groups2plot",
      label = "Groups to plot",
      choices = levels(rval_voi()),
      selected = levels(rval_voi())
    )
    updateSelectInput(
      session,
      "select_limma_contrasts2plot",
      label = "Contrasts to plot",
      choices = rval_contrasts(),
      selected = rval_contrasts()
    )

    #force rval_filteredlist
    rval_filteredlist()
  })
  
  #Calculation of filtered list
  rval_filteredlist = reactive({
    withProgress(message = "Performing contrasts calculations...",
                 value = 1,
                 max = 6,
                 {
                   setProgress(message = "Calculating global difs...", value = 1)
                   rval_globaldifs()
                   setProgress(message = "Calculating eBayes...", value = 4)
                   rval_finddifcpgs()
                   setProgress(message = "Calculating filtered list...", value = 5)
                   create_filtered_list(
                     rval_finddifcpgs(),
                     rval_globaldifs(),
                     deltaB = input$slider_limma_deltab,
                     adjp_max = input$slider_limma_adjpvalue,
                     p.value = input$slider_limma_pvalue,
                     cores = n_cores
                   )
                 })
  })
  
  
  #render of plots and tables
  output$graph_limma_plotSA = renderPlot(limma::plotSA(rval_fit()))
  output$graph_limma_plotMA = renderPlot(limma::plotMA(rval_gset_getM()))
  
  
  plot_heatmap = eventReactive(
    input$button_limma_heatmapcalc,
    create_heatmap(
      rval_filteredlist(),
      rval_gset_getBeta(),
      factorgroups =  factor(rval_voi()[rval_voi() %in% input$select_limma_groups2plot],
                             levels =
                               input$select_limma_groups2plot),
      groups2plot = rval_voi() %in% input$select_limma_groups2plot,
      contrasts2plot = rval_contrasts() %in% input$select_limma_contrasts2plot,
      Colv = as.logical(input$select_limma_colv),
      clusteralg = input$select_limma_clusteralg,
      distance = input$select_limma_clusterdist,
      scale = input$select_limma_scale,
      static = as.logical(input$select_limma_graphstatic)
    )
  )
  make_table = eventReactive(
    input$button_limma_heatmapcalc,
    data.table::rbindlist(rval_filteredlist()[rval_contrasts() %in% input$select_limma_contrasts2plot],
                          idcol = "contrast") %>% 
                          dplyr::mutate(type = ifelse(.data$dif_current < 0, "Hypermethylated", "Hypomethylated")) %>% 
                            dplyr::group_by(.data$contrast, .data$type)  %>%
                            dplyr::summarise(CpGs = dplyr::n()) %>% tidyr::pivot_wider(names_from =
                            .data$type, values_from = .data$CpGs) %>% dplyr::mutate(total = .data$Hypermethylated + .data$Hypomethylated)
  )
  
  observeEvent(input$button_limma_heatmapcalc, {
    #Render the correct plot depending on the selected
    output$graph_limma_heatmapcontainer = renderUI({
      if (!as.logical(input$select_limma_graphstatic))
        return(
          plotly::plotlyOutput(
            "graph_limma_heatmap_interactive",
            width = "600px",
            height = "500px"
          )  %>% shinycssloaders::withSpinner()
        )
      else
        return(
          plotOutput(
            "graph_limma_heatmap_static",
            width = "600px",
            height = "630px"
          ) %>% shinycssloaders::withSpinner()
        )
    })
    
    
  
    if (!as.logical(input$select_limma_graphstatic))
      output$graph_limma_heatmap_interactive = plotly::renderPlotly(plot_heatmap())
    else
      output$graph_limma_heatmap_static = renderPlot(plot_heatmap())
  })
  
  output$table_limma_difcpgs = renderTable(make_table())
  
  

  
  
  #EXPORT
  
  rval_annotation = reactive(minfi::getAnnotation(rval_gset()))
  
  #R Objects
  output$download_export_robjects = downloadHandler(
    "R_Objects.zip",
    content = function(file) {
      
      rgset = rval_rgset()
      gset = rval_gset()
      fit = rval_fit()
      design = rval_design()
      ebayes_tables = rval_finddifcpgs()
      bvalues = rval_gset_getBeta()
      mvalues = rval_gset_getM()
      #annotation = rval_annotation()
      global_difs = rval_globaldifs()
      
      save(bvalues,file=paste(tempdir(),"Bvalues.RData", sep="/"))
      save(mvalues,file=paste(tempdir(),"Mvalues.RData", sep="/"))
      save(rgset , file = paste(tempdir(),"RGSet.RData", sep="/"))
      save(gset, file=paste(tempdir(),"Normalized_genomicratioset.RData", sep="/"))
      save(fit, file=paste(tempdir(),"fit.RData", sep="/"))
      save(design, file=paste(tempdir(),"design.RData", sep="/"))
      save(global_difs, file=paste(tempdir(),"global_difs.RData", sep="/"))
      #save(annotation,file=paste(tempdir(),"annotation.RData", sep="/"))
      save(ebayes_tables,file=paste(tempdir(),"ebayestables.RData", sep="/"))
      
      objects = list.files(tempdir(), pattern = "*.RData", full.names = TRUE)
      
      utils::zip(file, objects, flags = "-j9X")
      
    }
  )
  
  
  #Filtered BEDs
  output$download_export_filteredbeds  = downloadHandler("filtered_beds.zip",
    content = function(file) {
      
      filtered_beds = create_filtered_beds(rval_filteredlist(), rval_annotation(), cores=n_cores)

      
      lapply(names(filtered_beds), function(x) {
        data.table::fwrite(
          filtered_beds[[x]],
          file = paste0(dirname(file), "/", x, ".bed"),
          sep = "\t",
          quote = FALSE,
          col.names = FALSE,
          row.names = FALSE
        )
      })
      
      objects = list.files(path = dirname(file),
                           full.names = TRUE,
                           pattern = "*.bed")
      
      utils::zip(file, objects, flags = "-j9X")
    }
  )
  
  
  #Markdown Report
  
  output$download_export_markdown = downloadHandler(
    filename= "Markdown.html",
    content = function(file) {
       params = list(
        name_var = input$select_input_samplenamevar,
        grouping_var = input$select_input_groupingvar,
        donor_var = input$select_input_donorvar,
        normalization_mode = input$select_minfi_norm,
        rval_rgset = rval_rgset(),
        rval_gset = rval_gset(),
        rval_sheet_target = rval_sheet_target(),
        rval_fit = rval_fit(),
        rval_design = rval_design(),
        rval_filteredlist = rval_filteredlist(),
        rval_contrasts = rval_contrasts(),
        limma_ebayes_trend = input$select_limma_trend,
        limma_ebayes_robust = input$select_limma_robust,
        limma_voi = input$select_limma_voi,
        limma_covar = input$checkbox_limma_covariables,
        groups2plot = input$select_limma_groups2plot,
        contrasts2plot = input$select_limma_contrasts2plot,
        Colv = input$select_limma_colv,
        clusteralg = input$select_limma_clusteralg,
        distance = input$select_limma_clusterdist,
        scale = input$select_limma_scale,
        max_fdr = input$slider_limma_adjpvalue,
        min_deltabeta = input$slider_limma_deltab,
        max_pvalue = input$slider_limma_pvalue
      )
      
      print(file)
      newenv = new.env(parent = globalenv())
      newenv$create_heatmap = create_heatmap
      
      rmarkdown::render(
        input = system.file("report.Rmd", package = "shinyepico"),
        output_file = file,
        params = params,
        envir = newenv
      )
      
    }
  )
  
  
  #custom heatmap
  
  output$download_export_heatmaps = downloadHandler(
    "Custom_heatmap.pdf",
    content = function(file) {
      grDevices::pdf(file = file,
                     height = 11.71,
                     width = 8.6)
      create_heatmap(
        rval_filteredlist(),
        rval_gset_getBeta(),
        factorgroups =  factor(rval_voi()[rval_voi() %in% input$select_limma_groups2plot], levels =
                                 input$select_limma_groups2plot),
        groups2plot = rval_voi() %in% input$select_limma_groups2plot,
        contrasts2plot = rval_contrasts() %in% input$select_limma_contrasts2plot,
        Colv = as.logical(input$select_limma_colv),
        clusteralg = input$select_limma_clusteralg,
        distance = input$select_limma_clusterdist,
        scale = input$select_limma_scale,
        static = TRUE
      )
      grDevices::dev.off()
      
    }
    
  )
  
  #Graphs
  
  
  
  
}
