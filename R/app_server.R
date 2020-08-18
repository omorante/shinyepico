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
  
  #INITIALIZE REACTIVE VARIABLES
  rval_generated_limma_model = reactiveVal(value = FALSE)
  
  #Load button only shows if file is uploaded
  output$ui_button_input_load = renderUI({
    if (!is.null(input$fileinput_input$datapath))
      return(actionButton("button_input_load", "Load Data"))
    else
      return()
  })
  
  #When you press button_input_load, the data is unzipped and the metharray sheet is loaded
  rval_sheet = eventReactive(input$button_input_load, {

    utils::unzip(input$fileinput_input$datapath,
                 exdir = paste0(tempdir(), "/experiment_data"))
    
    sheet = minfi::read.metharray.sheet(paste0(tempdir(), "/experiment_data"))

    #We check if sheet is correct
    #This is to prevent app crashes when zip updated is not correct.
    validate(need(is.data.frame(sheet) & any(colnames(sheet) %in% "Slide") & any(colnames(sheet) %in% "Array"), 
                  "SampleSheet is not correct. Please, check your samplesheet and your zip file." ))
    
    sheet
  })
  
  
  rval_sheet_target = eventReactive(input$button_input_load,
        rval_sheet()[rval_sheet()[, input$select_input_samplenamevar]  %in% input$selected_samples, ])
  
  
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
  

  #The checkbox of samples to process is updated when samplenamevar changes
  observeEvent(
    {input$select_input_samplenamevar
      input$select_input_groupingvar},
    
    updatePickerInput(
      session,
      "selected_samples",
      label = "Select Samples to Process:",
      selected = rval_sheet()[, input$select_input_samplenamevar],
      choices = rval_sheet()[, input$select_input_samplenamevar],
      choicesOpt = list(subtext = paste("Group: ", rval_sheet()[, input$select_input_groupingvar]))
    )
  )
  

  #The dataframe is rendered
  output$samples_table =  DT::renderDT(rval_sheet())
  
  #rval_rgset loads RGSet using read.metharray.exp and the sample sheet (rval_sheet())
  
  rval_rgset = eventReactive(input$button_input_next, {
    targets = rval_sheet()[rval_sheet()[, input$select_input_samplenamevar]  %in% input$selected_samples, ]
    
    #Check prior conditions to read data
    validate(need(anyDuplicated(rval_sheet_target()[,input$select_input_samplenamevar]) == 0, 
                  "Sample Name Variable should not have duplicated values"))
    validate(need(anyDuplicated(rval_sheet_target()[,input$select_input_groupingvar]) > 0, 
                  "Grouping variable should have groups greater than 1"))
    
    #We need to check if this step works
    try({RGSet = minfi::read.metharray.exp(targets = targets, verbose = TRUE, force = TRUE)})
    
    if (!exists("RGSet", inherits = FALSE)){
      showModal(modalDialog(
        title = "reading error",
        "Minfi can't read arrays specified in your samplesheet. Please, check your zipfile and your sampleSheet",
        easyClose = TRUE,
        footer = NULL))
    }
    
    validate(need(exists("RGSet", inherits = FALSE), "Minfi can't read arrays specified in your samplesheet. Please, check your zipfile and your sampleSheet"))  
    
    #We return RGSet filter by the standard detection threshold of Pvalue, 0.01
    RGSet [(rowMeans(as.matrix(minfi::detectionP(RGSet))) < 0.01), ]
  })
  
  #We change the page to the next one
  observeEvent(input$button_input_next,
      { 
        
        #check if variables selected are correct
        if(anyDuplicated(rval_sheet_target()[,input$select_input_samplenamevar]) > 0 | 
           anyDuplicated(rval_sheet_target()[,input$select_input_groupingvar]) == 0 ){
          
          showModal(modalDialog(
            title = "Variable error",
            "Check if selected variables are correct. Sample Name Variable should not have duplicated values 
          and grouping variable should have groups greater than 1.",
            easyClose = TRUE,
            footer = NULL))
        }
        
        updateSelectInput(session, "select_minfi_pcaplot_color", choices = colnames(rval_sheet()), selected = input$select_input_donorvar)
        
        withProgress(message = "Loading data...", value = 2, max = 5,
                              {rval_rgset()
                                updateNavbarPage(session, "navbar_epic", "Normalization")
                              })
      })
  
  
  #MINFI NORMALIZATION
  
  
  #Calculation of minfi normalized data
  rval_gset = reactive({
    
    withProgress(message = "Normalization in progress..",
                 value = 1,
                 max = 4,
    {
      
    if (is.null(rval_rgset())){
      stop() # this function doesn't continue if rval_rgset doesn't exist.
    }
    
   try({
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
      gset = minfi::preprocessFunnorm(rval_rgset())
    }
    
    else if (input$select_minfi_norm == "SWAN") {
      gset = minfi::preprocessSWAN(rval_rgset(), mSet = minfi::mapToGenome(minfi::preprocessRaw(rval_rgset()))) #MethylSet or GenomicMethylSet?
    }
    
    else if (input$select_minfi_norm == "Quantile") {
      gset = minfi::preprocessQuantile(rval_rgset())
    }
    
    else if (input$select_minfi_norm == "Noob+Quantile") {
      gset = minfi::preprocessQuantile(minfi::preprocessNoob(rval_rgset()))
    }
     
    
   })
    
    #check if normalization has worked
    validate(need(exists("gset", inherits = FALSE), "An unexpected error has occurred during minfi normalization. Please, notify the error to the package maintainer."))  
    
    #remove SNPs to proceed with the analysis and add sex column
    minfi::addSex(minfi::dropLociWithSnps(gset, maf = 0))
    
    })
    
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
  ##############
  
  # Minfi Graphics
  
  #Density plots
  output$graph_minfi_densityplotraw = plotly::renderPlotly(create_densityplot(rval_rgset_getBeta(), 200000))
  output$graph_minfi_densityplot =  plotly::renderPlotly(create_densityplot(rval_gset_getBeta(), 200000))
  
  #PCA
  rval_plot_pca = reactive(create_pca(Bvalues=rval_gset_getBeta(), pheno_info=rval_sheet_target(), group=input$select_input_samplenamevar, 
                                      pc_x=input$select_minfi_pcaplot_pcx, pc_y=input$select_minfi_pcaplot_pcy, color=input$select_minfi_pcaplot_color))
  
  output$graph_minfi_pcaplot = plotly::renderPlotly(rval_plot_pca()[["graph"]])
  output$table_minfi_pcaplot = renderTable(rval_plot_pca()[["info"]], rownames = TRUE)
  
  
  #Correlations
  
  output$graph_minfi_corrplot = plotly::renderPlotly(create_corrplot(rval_gset_getBeta(), rval_sheet_target()))
  
  #Boxplots
  output$graph_minfi_boxplotraw = renderCachedPlot(graphics::boxplot(as.matrix(rval_rgset_getBeta())),
                                                  paste0(input$select_minfi_norm, "boxplot", input$selected_samples))
  
  output$graph_minfi_boxplot = renderCachedPlot(graphics::boxplot(as.matrix(rval_gset_getBeta())),
                                               paste0(input$select_minfi_norm, "boxplot", input$selected_samples))
  

  
  
  #QC plots
  output$graph_minfi_qcraw = plotly::renderPlotly(create_plotqc(rval_rgset(), rval_sheet_target()[,input$select_input_samplenamevar])) 
  output$graph_minfi_bisulfiterawII = plotly::renderPlotly(create_bisulfiteplot(rval_rgset(), rval_sheet_target()[,input$select_input_samplenamevar]))

  #Sex prediction
  #output$graph_minfi_sex = renderCachedPlot(minfi::plotSex(rval_gset(), 
  #                                          id = rval_sheet_target()[, input$select_input_samplenamevar] ),
  #                                          paste0(input$select_minfi_norm, "sex", input$selected_samples))
  
  output$graph_minfi_sex = plotly::renderPlotly(create_sexplot(rval_gset(), rval_sheet_target()[,input$select_input_samplenamevar]))
  
  #SNPs heatmap
  output$graph_minfi_snps = plotly::renderPlotly(create_snpheatmap(minfi::getSnpBeta(rval_rgset()), 
                                                 rval_sheet_target()[, input$select_input_samplenamevar],
                                                 rval_sheet_target()[, input$select_input_donorvar]))
  
  
  
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
                 
                 updatePickerInput(
                   session,
                   "checkbox_limma_covariables",
                   label = "Select linear model covariables",
                   choices = colnames(minfi::pData(rval_gset())),
                   selected = input$select_input_donorvar,
                 )
                 

                 updateNavbarPage(session, "navbar_epic", "DMPs")
                 
                 #preparing next form
                 updateSwitchInput(session, "select_limma_trend", 
                                   disabled=FALSE)
                 
                 updateSwitchInput(session, "select_limma_trend", 
                                   disabled=FALSE)
                 
               })
  
  
  
  
  #LIMMA
  
  #Variable of interest
  rval_voi = reactive(factor(sub("-","_",minfi::pData(rval_gset())[, input$select_limma_voi]))) #add substitution of "-" for "_", avoiding conflicts

  #Design calculation
  rval_design = reactive({

    pdata = data.table::as.data.table(apply(minfi::pData(rval_gset()), 2, sub, pattern = "-", replacement = "_")) #avoiding any "-" in the data
    
    #We convert to numeric variable any column of pdata with all numbers
    suppressWarnings({
    pdata = data.table::as.data.table(lapply(pdata, function(x){
      if (!any(is.na(as.numeric(x)))) {
        as.numeric(x)
      }
      else x
      }))
    })
    
    formula = stats::as.formula(paste0("~ 0 + ", paste(
      c(
        input$select_limma_voi,
        input$checkbox_limma_covariables
      ),
      collapse = "+"
    )))
    
    validate(need(anyDuplicated(pdata[,input$select_limma_voi]) > 0, "Variable selected is not valid. Please, change the variable and try again"))
    
    if (length(input$checkbox_limma_covariables > 0)){
    validate(need(all(apply(pdata[,input$checkbox_limma_covariables, drop=FALSE], 2, anyDuplicated)) > 0, 
                  "Covariable(s) selected are not valid. Please, change the covariable and try again" ))
    }
    
    #Bulding the design matrix
    design = stats::model.matrix(formula, data = pdata)
    print(design)
    colnames(design)[seq_len(length(unique(rval_voi())))] = levels(rval_voi())
    
    design
  })
  
  #Calculation of contrasts
  rval_contrasts = reactive({
    contrastes = c() #calculamos lista de contrastes de todos contra todos
    valor = 1
    all.groups = levels(rval_voi())
    for (i in seq_len(length(all.groups) - 1)) {
      for (z in seq_len((length(all.groups)) - valor)) {
        contrastes = c(contrastes, paste(all.groups[i], all.groups[i + z], sep = "-"))
      }
      valor <- valor + 1
    }
    contrastes
  })
  
  
  #Calculation of limma model
  rval_fit = eventReactive(input$button_limma_calculatemodel, {
    
    withProgress(message = "Generating linear model...",
                 value = 3,
                 max = 6,
      {if (as.logical(input$select_limma_weights)){
        try({weights = limma::arrayWeights(rval_gset_getM(), design = rval_design())})
      }
        else { weights = NULL}
    
      try({fit = limma::lmFit(rval_gset_getM(), rval_design(), weights = weights)})
    
      if (!exists("fit", inherits = FALSE)){
        rval_generated_limma_model(FALSE) #disable contrast button
        showModal(modalDialog(
          title = "lmFit error",
          "lmFit model has failed. Please, check your options and try again.",
          easyClose = TRUE,
          footer = NULL))
          }
    })
    
    validate(need(exists("fit", inherits = FALSE), "lmFit model has failed. Please, check your options and try again.")) 
    
    fit #returning linear model
  })
  
  
  # rval_fit() has NAs, we remove the option to trend or robust in eBayes to prevent failure
  observeEvent(input$button_limma_calculatemodel, {
    
    if (any(vapply(rval_fit(), function(x){any(is.na(unlist(x))|unlist(x)==Inf|unlist(x)==-Inf)}, logical(1)))) {
      
      message("NAs or Inf values detected, trend and robust options are disabled.")
      
      updateSwitchInput(session, "select_limma_trend", 
                        value = FALSE, 
                        disabled=TRUE)
      
      updateSwitchInput(session, "select_limma_robust", 
                        value = FALSE, 
                        disabled=TRUE)
    }
    
    else {
      
      message("NAs or Inf values not detected, trend and robust options are enabled")
      
      updateSwitchInput(session, "select_limma_trend", 
                        value = FALSE, 
                        disabled = FALSE)
      
      updateSwitchInput(session, "select_limma_robust", 
                        value = FALSE, 
                        disabled=FALSE)
      
    }
    
    # Creating calculate differences button
    rval_generated_limma_model(TRUE)
  })
  
  
  output$button_limma_calculatedifs_container = renderUI({

    if (rval_generated_limma_model())
      return(tagList(
        
        br(),
        
        switchInput(
          inputId = "select_limma_trend",
          label = "eBayes Trend", 
          labelWidth = "80px",
          value = FALSE),
        
        switchInput(
          inputId = "select_limma_robust",
          label = "eBayes Robust", 
          labelWidth = "80px",
          value = FALSE),
        
        actionButton("button_limma_calculatedifs", "Calculate Contrasts")))
    
    else 
      return()
      })
  
  
  #render of plots and tables
  #output$graph_limma_plotMA = renderPlot(limma::plotMA(rval_gset_getM(), main = ""))
  output$graph_limma_plotSA = renderPlot(limma::plotSA(rval_fit()))
  
  
  
  #Calculation of global difs
  rval_globaldifs = eventReactive(input$button_limma_calculatedifs, {
    calculate_global_difs(rval_gset_getBeta(), rval_voi(), rval_contrasts(), cores=n_cores)
  })
  
  #Calculation of differences (eBayes)
  rval_finddifcpgs = eventReactive(input$button_limma_calculatedifs, {
    
    try({dif_cpgs = find_dif_cpgs(
      rval_voi(),
      rval_design(),
      rval_fit(),
      rval_contrasts(),
      trend = as.logical(input$select_limma_trend),
      robust = as.logical(input$select_limma_robust),
      cores = n_cores)
      })
    
    if (!exists("dif_cpgs", inherits = FALSE)){
      
      showModal(modalDialog(
        title = "Contrasts Calculation Error",
        "An unexpected error has ocurred during contrasts calculation. Please, generate the model again. 
        If the problem persists, report the error to the mantainer.",
        easyClose = TRUE,
        footer = NULL))
  
      rval_generated_limma_model(FALSE)
    }
    
    validate(need(exists("dif_cpgs", inherits = FALSE), 
                  "An unexpected error has ocurred during contrasts calculation. Please, generate the model again. 
                  If the problem persists, report the error to the mantainer."))
    
    dif_cpgs
    
    
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
  
  
  rval_filteredlist2heatmap = reactive({
    
    filtered_data = rval_filteredlist()[rval_contrasts() %in% input$select_limma_contrasts2plot] # filter contrasts2plot
    dif_cpgs = unique(data.table::rbindlist(filtered_data)$cpg)
    join_table = rval_gset_getBeta()[dif_cpgs, ]
    join_table$cpg = NULL
    
    validate(need(!(is.null(join_table)) & nrow(join_table) > 0,  "No differences were found with these criteria" ))
    validate(need(nrow(join_table) <= 10000, "Too many differences with these criteria (>10000), they can not be plotted" )) 
    
    join_table
  })
  
  
  plot_heatmap = eventReactive(input$button_limma_heatmapcalc,
    
     {

      
      create_heatmap(
      rval_filteredlist2heatmap(),
      factorgroups =  factor(rval_voi()[rval_voi() %in% input$select_limma_groups2plot],
                             levels =
                               input$select_limma_groups2plot),
      groups2plot = rval_voi() %in% input$select_limma_groups2plot,
      Colv = as.logical(input$select_limma_colv),
      clusteralg = input$select_limma_clusteralg,
      distance = input$select_limma_clusterdist,
      scale = input$select_limma_scale,
      static = as.logical(input$select_limma_graphstatic)
    )
   }
  )
  make_table = eventReactive(
    input$button_limma_heatmapcalc,
    data.table::rbindlist(rval_filteredlist()[rval_contrasts() %in% input$select_limma_contrasts2plot],
                          idcol = "contrast") %>% 
                          dplyr::mutate(type = factor(ifelse(.data$dif_current < 0, "Hypermethylated", "Hypomethylated"), levels = c("Hypermethylated","Hypomethylated"))) %>% 
                            dplyr::group_by(.data$contrast, .data$type)  %>%
                            dplyr::summarise(CpGs = dplyr::n()) %>% 
                            tidyr::complete(.data$contrast, .data$type, fill=list(CpGs=0)) %>%
                            tidyr::pivot_wider(names_from = .data$type, values_from = .data$CpGs) %>% 
                            dplyr::mutate(total = .data$Hypermethylated + .data$Hypomethylated) %>%
                            dplyr::mutate(dplyr::across(c("Hypermethylated","Hypomethylated", "total"), ~ format(., digits=0)))
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
      
      withProgress(message = "Generating R Objects...",
                   value = 1,
                   max = 4,
      {
        
      rgset = rval_rgset()
      gset = rval_gset()
      fit = rval_fit()
      design = rval_design()
      ebayes_tables = rval_finddifcpgs()
      bvalues = rval_gset_getBeta()
      mvalues = rval_gset_getM()
      #annotation = rval_annotation()
      global_difs = rval_globaldifs()
      
      
      setProgress(value=2, message = "Saving RObjects...")
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
      
      setProgress(value=3, message = "Compressing RObjects...")
      utils::zip(file, objects, flags = "-j9X")
        })
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
      
      withProgress(message = "Generating Report...",
                   value = 1,
                   max = 3,
       {
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
        rval_filteredlist2heatmap = rval_filteredlist2heatmap(),
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
      
      })
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
        rval_filteredlist2heatmap(),
        factorgroups =  factor(rval_voi()[rval_voi() %in% input$select_limma_groups2plot], levels =
                                 input$select_limma_groups2plot),
        groups2plot = rval_voi() %in% input$select_limma_groups2plot,
        Colv = as.logical(input$select_limma_colv),
        clusteralg = input$select_limma_clusteralg,
        distance = input$select_limma_clusterdist,
        scale = input$select_limma_scale,
        static = TRUE
      )
      grDevices::dev.off()
    }
  )
}
