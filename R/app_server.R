#' The application server-side
#'
#' @param input,output,session Internal parameters for {shiny}.
#'     DO NOT REMOVE.
#' @import shiny
#' @importFrom rlang .data
#' @import data.table
#'
#' @noRd

app_server = function(input, output, session) {
  #PERFORMANCE SETTINGS
  n_cores = getShinyOption("n_cores")
  options(shiny.maxRequestSize = getShinyOption("shiny.maxRequestSize"))
  
  #INITIALIZE REACTIVE VARIABLES
  rval_generated_limma_model = reactiveVal(value = FALSE)
  rval_analysis_finished = reactiveVal(value = FALSE)
  rval_filteredlist2heatmap_valid = reactiveVal(value = FALSE)
  rval_filteredmcsea2heatmap_valid = reactiveVal(value = FALSE)
  rval_dmrs_finished = reactiveVal(value = FALSE)
  rval_dmrs_ready2heatmap = reactiveVal(value = FALSE)
  rval_dmrs_ready2mcsea = reactiveVal(value = FALSE)
  
  #Load button only shows if file is uploaded
  output$ui_button_input_load = renderUI({
    if (!is.null(input$fileinput_input$datapath))
      return(actionButton("button_input_load", "Load Data"))
    else
      return()
  })
  
  # Enable load button every time file input is updated
  observeEvent(input$fileinput_input, shinyjs::enable("button_input_load"))
  
  #When you press button_input_load, the data is unzipped and the metharray sheet is loaded
  rval_sheet = eventReactive(input$button_input_load, {
    
    #Check if updated file is .zip
    validate(need(tools::file_ext(input$fileinput_input$datapath) == "zip", "File extension should be .zip"))
    
    shinyjs::disable("button_input_load") #disable the load button to avoid multiple clicks
    
    if (dir.exists(paste0(tempdir(), "/experiment_data"))) {
      unlink(paste0(tempdir(), "/experiment_data"), recursive = TRUE) #remove current files in target directory
    }
        
    zip::unzip(input$fileinput_input$datapath,
                 exdir = paste0(tempdir(), "/experiment_data")) #extracting zip
    
    sheet = minfi::read.metharray.sheet(paste0(tempdir(), "/experiment_data"))
    
    #We check if sheet is correct
    #This is to prevent app crashes when zip updated is not correct.
    validate(
      need(
        is.data.frame(sheet) &
          any(colnames(sheet) %in% "Slide") &
          any(colnames(sheet) %in% "Array"),
        "SampleSheet is not correct. Please, check your samplesheet and your zip file."
      )
    )
    
    validate(
      need(
        anyDuplicated(colnames(sheet)) == 0,
        "Repeated variable names are not allowed. Please, modify your sample sheet."
      )
    )
    
    colnames(sheet) = make.names(colnames(sheet)) # fix possible not-valid colnames
    
    sheet
  })
  
  
  rval_sheet_target = eventReactive(input$button_input_next,
                                    rval_sheet()[rval_sheet()[, input$select_input_samplenamevar]  %in% input$selected_samples,])
  
  
  rval_clean_sheet_target = eventReactive(rval_gset(),{
    
    suppressWarnings({
      clean_sample_sheet = as.data.frame(lapply(minfi::pData(rval_gset()), function(x) {
        if (sum(is.na(as.numeric(x))) < length(x) * 0.75 & stats::sd(x, na.rm=TRUE) > 0) {
          return(as.numeric(x))
        } #if NAs produced with as.numeric are less than 75%, and SD is greater than 0 we consider the variable numeric
        else if (length(unique(x)) > 1 &
                 length(unique(x)) < length(x)) {
              return(as.factor(make.names(x))) #returning syntactically valid names
        } #if the variable is a character, it should have unique values more than 1 and less than total
        else{
          return(rep(NA, length(x)))
        } #if the requeriments are not fullfilled, we discard the variable
      }),
      stringsAsFactors = FALSE)
    })
    
    clean_sample_sheet[["Slide"]] = as.factor(clean_sample_sheet[["Slide"]]) #adding Slide as factor
    clean_sample_sheet[[input$select_input_donorvar]] = as.factor(clean_sample_sheet[[input$select_input_donorvar]]) #adding donorvar as factor
    clean_sample_sheet = clean_sample_sheet[, colSums(is.na(clean_sample_sheet)) < nrow(clean_sample_sheet)] #cleaning all NA variables

    clean_sample_sheet
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
      label = "Select Variable of Interest:",
      choices = colnames(rval_sheet())
    )
    updateSelectInput(
      session,
      "select_input_donorvar",
      label = "Select Donor Variable:",
      choices = colnames(rval_sheet())
    )
    
    shinyjs::enable("button_input_next") #Enable button continue
  })
  
  
  #The checkbox of samples to process is updated when samplenamevar changes
  observeEvent({
    input$select_input_samplenamevar
    input$select_input_groupingvar
  },
  
  updatePickerInput(
    session,
    "selected_samples",
    label = "Select Samples to Process:",
    selected = rval_sheet()[, input$select_input_samplenamevar],
    choices = rval_sheet()[, input$select_input_samplenamevar],
    choicesOpt = list(subtext = paste("Group: ", rval_sheet()[, input$select_input_groupingvar]))
  ))
  
  #when samples selected are changed, continue button is enabled again
  observeEvent(input$selected_samples, shinyjs::enable("button_input_next"))
  
  #The dataframe is rendered
  output$samples_table =  DT::renderDT(
    rval_sheet(),
    rownames = FALSE,
    selection = "single",
    style = "bootstrap",
    options = list(
      pageLength = 10,
      autoWidth = TRUE,
      scrollX = TRUE,
      columnDefs = list(list(
       targets = match("Basename",colnames(rval_sheet())) - 1, visible = FALSE
      ))
    )
  )
  
  
  
  
  #rval_rgset loads RGSet using read.metharray.exp and the sample sheet (rval_sheet())
  rval_rgset = eventReactive(input$button_input_next, ignoreNULL = FALSE, {
    
    validate(need(input$fileinput_input != "", "Data has not been uploaded yet"))
    
    #Prior check to test variable selection
    if (anyDuplicated(rval_sheet_target()[, input$select_input_samplenamevar]) > 0 |
        anyDuplicated(rval_sheet_target()[, input$select_input_groupingvar]) == 0) {
      showModal(
        modalDialog(
          title = "Variable error",
          "Check if selected variables are correct. Sample Name Variable should not have duplicated values
          and the variable of interest should have groups greater than 1.",
          easyClose = TRUE,
          footer = NULL
        )
      )
    }
    
    #Check prior conditions to read data
    validate(need(
      anyDuplicated(rval_sheet_target()[, input$select_input_samplenamevar]) == 0,
      "Sample Name Variable should not have duplicated values"
    ))
    validate(need(
      anyDuplicated(rval_sheet_target()[, input$select_input_groupingvar]) > 0,
      "Grouping variable should have groups greater than 1"
    ))
    
    #disable button to avoid multiple clicks
    shinyjs::disable("button_input_next")
    
    
    #We need to check if this step works
    withProgress(message = "Reading array data...",
                 value = 2,
                 max = 5,
                 {
                   try({
                     RGSet = minfi::read.metharray.exp(targets = rval_sheet_target(),
                                                       verbose = TRUE,
                                                       force = TRUE)
                   })
                   
                   
                   
                   
                   if (!exists("RGSet", inherits = FALSE)) {
                     showModal(
                       modalDialog(
                         title = "reading error",
                         "Minfi can't read arrays specified in your samplesheet. Please, check your zipfile and your sampleSheet",
                         easyClose = TRUE,
                         footer = NULL
                       )
                     )
                     shinyjs::disable("button_minfi_select")
                   }
                   
                   validate(
                     need(
                       exists("RGSet", inherits = FALSE),
                       "Minfi can't read arrays specified in your samplesheet. Please, check your zipfile and your sampleSheet"
                     )
                   )
                   
                   #we return RGSet filtered by the standard detection threshold of Pvalue, 0.01
                   RGSet[(rowMeans(as.matrix(minfi::detectionP(RGSet))) < 0.01), ]
                 })
    
  })
  
  #We change the page to the next one
  observeEvent(input$button_input_next,
               {
                 #check if rgset is loaded
                 req(rval_rgset())
                 
                 #update PCA parameters
                 updateSelectInput(
                   session,
                   "select_minfi_pcaplot_pcx",
                   choices = paste0("PC", seq_len(nrow(rval_sheet_target()))),
                   selected = "PC1"
                 )
                 
                 updateSelectInput(
                   session,
                   "select_minfi_pcaplot_pcy",
                   choices = paste0("PC", seq_len(nrow(rval_sheet_target()))),
                   selected = "PC2"
                 )
                 
                 updateSelectInput(
                   session,
                   "select_minfi_pcaplot_color",
                   choices = c(colnames(rval_sheet_target()),
                               "xMed",
                               "yMed",
                               "predictedSex"),
                   
                   selected = input$select_input_groupingvar
                 )
                 
                 shinyjs::enable("button_minfi_select")
                 updateNavbarPage(session, "navbar_epic", "Normalization")
               })
  
  
  #MINFI NORMALIZATION
  
  
  #Calculation of minfi normalized data
  rval_gset = eventReactive(input$button_minfi_select, {
    
    validate(need(!is.null(rval_rgset()),
                  "Raw data has not been loaded yet."))
    
    shinyjs::disable("button_minfi_select") #disable button to avoid repeat clicking
    
    withProgress(message = "Normalization in progress...",
                 value = 1,
                 max = 4,
                 {
                   try({
                     gset = normalize_rgset(rgset = rval_rgset(), normalization_mode = input$select_minfi_norm)
                   })
                   
                   #check if normalization has worked
                   
                   if(!exists("gset", inherits = FALSE)) {
                     showModal(
                       modalDialog(
                         title = "Normalization failure",
                         "An unexpected error has occurred during minfi normalization. Please, notify the error to the package maintainer.",
                         easyClose = TRUE,
                         footer = NULL
                       )
                     )
                     shinyjs::enable("button_minfi_select") 
                   }
                   
                   validate(
                     need(
                       exists("gset", inherits = FALSE),
                       "An unexpected error has occurred during minfi normalization. Please, notify the error to the package maintainer."
                     )
                   )
                   
                   #prior CpG probes
                   probes = length(minfi::featureNames(gset))
                   
                   #remove SNPs to proceed with the analysis and add sex column
                   if (input$select_minfi_dropsnps)
                    gset = minfi::dropLociWithSnps(gset, maf = input$slider_minfi_maf)
                   
                   #remove CpHs
                   if (input$select_minfi_dropcphs)
                    gset = minfi::dropMethylationLoci(gset, dropRS = TRUE, dropCH = TRUE)
                   
                   #Add sex info
                    gset = minfi::addSex(gset)
                   
                   #remove chromosomes
                   if(input$select_minfi_chromosomes){
                     gset = gset[rownames(minfi::getAnnotation(gset))[!(minfi::getAnnotation(gset)$chr %in% c("chrX","chrY"))],]
                   }
                   
                   #Info CpGs removed
                   message(paste(probes - length(minfi::featureNames(gset)), "probes were filtered out."))
                   
                   #enable button
                   shinyjs::enable("button_minfi_select") 

                   #return gset
                   gset
                   
                 })
  })
  
  
  #filtered probes info
  
  rval_gsetprobes = eventReactive(input$button_minfi_select,{
    req(rval_gset())
    length(minfi::featureNames(rval_gset()))
  })
  
  output$text_minfi_probes = renderText(paste(rval_gsetprobes(),"positions after normalization"))

  #getBeta/getM reactives
  rval_rgset_getBeta = eventReactive(rval_rgset(), {
    bvalues = as.data.frame(minfi::getBeta(rval_rgset()))
    colnames(bvalues) = input$selected_samples
    bvalues
  })
  
  rval_gset_getBeta = eventReactive(rval_gset(), {
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
  rval_plot_densityplotraw = reactive(create_densityplot(rval_rgset_getBeta(), 200000))
  rval_plot_densityplot = reactive(create_densityplot(rval_gset_getBeta(), 200000))
  
  output$graph_minfi_densityplotraw = plotly::renderPlotly(rval_plot_densityplotraw())
  output$graph_minfi_densityplot = plotly::renderPlotly(rval_plot_densityplot())
  
  #PCA
  rval_plot_pca = eventReactive(
    list(input$button_pca_update, input$button_minfi_select),
    create_pca(
      Bvalues = rval_gset_getBeta(),
      pheno_info = as.data.frame(minfi::pData(rval_gset())),
      group = input$select_input_samplenamevar,
      pc_x = input$select_minfi_pcaplot_pcx,
      pc_y = input$select_minfi_pcaplot_pcy,
      color = input$select_minfi_pcaplot_color
    )
  )
  
  output$graph_minfi_pcaplot = plotly::renderPlotly(rval_plot_pca()[["graph"]])
  
  output$table_minfi_pcaplot = DT::renderDT(
    format(rval_plot_pca()[["info"]], digits = 2),
    rownames = TRUE,
    selection = "single",
    style = "bootstrap",
    options = list(
      autoWidth = TRUE,
      paging = FALSE,
      scrollX = TRUE,
      lengthChange = FALSE,
      searching = FALSE,
      info = FALSE
    )
  )
  
  #Correlations
  
  rval_plot_corrplot = reactive(
    create_corrplot(
      rval_gset_getBeta(),
      rval_clean_sheet_target(),
      rval_sheet_target(),
      p.value = input$select_minfi_typecorrplot == "p.value"
    )
  )
  
  output$graph_minfi_corrplot = plotly::renderPlotly(rval_plot_corrplot()[["graph"]])
  output$table_minfi_corrplot = DT::renderDT(rval_plot_corrplot()[["info"]],
                                             rownames = FALSE,
                                             selection = "single",
                                             style = "bootstrap",
                                             caption = "Autodetected variable types:",
                                             options = list(pageLength = 10, autoWidth = TRUE))
  
  #Boxplots
  
  rval_plot_boxplotraw = reactive(create_boxplot(rval_rgset_getBeta()))
  rval_plot_boxplot = reactive(create_boxplot(rval_gset_getBeta()))
  
  output$graph_minfi_boxplotraw = renderPlot(rval_plot_boxplotraw())
  output$graph_minfi_boxplot = renderPlot(rval_plot_boxplot())
  
  
  #QC plots
  
  rval_plot_qcraw = reactive(create_plotqc(rval_rgset(), rval_sheet_target()[, input$select_input_samplenamevar]))
  rval_plot_bisulfiterawII = reactive(create_bisulfiteplot(rval_rgset(), rval_sheet_target()[, input$select_input_samplenamevar]))
  
  output$graph_minfi_qcraw = plotly::renderPlotly(rval_plot_qcraw())
  output$graph_minfi_bisulfiterawII = plotly::renderPlotly(rval_plot_bisulfiterawII())
  
  #Sex prediction
  
  rval_plot_sexprediction = reactive({
    req(rval_gset())
    create_sexplot(rval_gset(), rval_sheet_target()[, input$select_input_samplenamevar])
  })
  
  rval_plot_sextable = reactive({
    req(rval_gset())
    data.frame(name = rval_sheet_target()[[input$select_input_samplenamevar]], sex = as.data.frame(minfi::pData(rval_gset()))[["predictedSex"]])
  })
  
  output$graph_minfi_sex = plotly::renderPlotly(rval_plot_sexprediction())
  
  output$table_minfi_sex = DT::renderDT(
    rval_plot_sextable(),
    rownames = FALSE,
    selection = "single",
    style = "bootstrap",
    caption = "Predicted sex:",
    options = list(
      pageLength = 10,
      scrollX = TRUE,
      autoWidth = TRUE
    )
  )
  
  
  #SNPs heatmap
  
  rval_plot_snpheatmap = reactive(
    create_snpheatmap(
      minfi::getSnpBeta(rval_rgset()),
      rval_sheet_target()[, input$select_input_samplenamevar],
      rval_sheet_target()[, input$select_input_donorvar]
    )
  )
  
  output$graph_minfi_snps = plotly::renderPlotly(rval_plot_snpheatmap())
  
  
  
  #Update of next form and move to Limma
  observeEvent(input$button_minfi_select,
               {
                 #check if normalization has been performed
                 req(rval_gset())
                 
                 updateSelectInput(
                   session,
                   "select_limma_voi",
                   label = "Select Variable of Interest",
                   choices = colnames(rval_clean_sheet_target())[vapply(rval_clean_sheet_target(), is.factor, logical(1))], 
                   selected = input$select_input_groupingvar
                 )
                 
                 updatePickerInput(
                   session,
                   "checkbox_limma_covariables",
                   label = "Select linear model covariables",
                   choices = colnames(rval_clean_sheet_target()),
                   selected = input$select_input_donorvar
                 )
                 
                 #possible interactions of variables:
                 if (length(colnames(rval_clean_sheet_target()) > 2)) {
                   interactions = utils::combn(colnames(rval_clean_sheet_target()), 2)
                   interactions = sprintf('%s:%s', interactions[1, ], interactions[2, ])
                 }
                 else {
                   interactions = c()
                 }
                 
                 
                 updatePickerInput(
                   session,
                   "checkbox_limma_interactions",
                   label = "Select linear model interactions",
                   choices = interactions,
                   selected = input$select_input_donorvar
                 )
                 
                 shinyjs::enable("button_limma_calculatemodel")
               })
  
  
  
  
  #LIMMA
  
  #Variable of interest
  rval_voi = reactive(factor(make.names(minfi::pData(rval_gset(
  ))[, input$select_limma_voi])))
  
  #Design calculation
  rval_design = eventReactive(input$button_limma_calculatemodel, {
    
    req(input$select_limma_voi) # a variable of interest is required
  
    #Bulding the design matrix
    try({
      design = generate_design(voi=input$select_limma_voi, sample_name=input$select_input_samplenamevar,
                               covariables=input$checkbox_limma_covariables, interactions=input$checkbox_limma_interactions, 
                               sample_sheet=rval_clean_sheet_target())
    })
    
    validate(
      need(
        exists("design", inherits = FALSE),
        "Design matrix has failed. Please, check your options and try again."
      ),
      need(
        nrow(design) == length(rval_sheet_target()[[input$select_input_samplenamevar]]),
        "The design matrix contains missing values. Please, check the selected variable and covariables."
      )
    )
    
    #checking colinearity
    if (qr(design)$rank < ncol(design)) {
      showModal(
        modalDialog(
          title = "Colinearity warning",
          "The design matrix presents colinear columns. Even though it is possible to try to continue the analysis, checking the selected covariables is recommended.",
          easyClose = TRUE,
          footer = NULL
        )
      )
      
    }
    
    design
  })
  
  #Calculation of contrasts
  rval_contrasts = reactive({
    conts = utils::combn(levels(rval_voi()), 2)
    
    sprintf('%s-%s', conts[1,], conts[2,])
  })
  
  
  #Calculation of limma model
  rval_fit = eventReactive(input$button_limma_calculatemodel, ignoreNULL = FALSE, {
    
    validate(
      need(input$fileinput_input != "", "DMP calculation has not been performed or data has not been uploaded."))
    
    req(rval_design())
    
    shinyjs::disable("button_limma_calculatemodel") #disable button to avoid repeat clicking
    
    withProgress(message = "Generating linear model...",
                 value = 3,
                 max = 6,
                 {
                   try({
                     fit = generate_limma_fit(Mvalues = rval_gset_getM(), design = rval_design(),
                                        weighting = as.logical(input$select_limma_weights) )
                   })
                   
                   if (!exists("fit", inherits = FALSE)) {
                     rval_generated_limma_model(FALSE) #disable contrast button
                     rval_analysis_finished(FALSE) #disable download buttons
                     
                     showModal(
                       modalDialog(
                         title = "lmFit error",
                         "lmFit model has failed. Please, check your options and try again.",
                         easyClose = TRUE,
                         footer = NULL
                       )
                     )
                   }
                 })
    
    shinyjs::enable("button_limma_calculatemodel") #enable button again to allow repeting calculation
    
    validate(need(
      exists("fit", inherits = FALSE),
      "lmFit model has failed. Please, check your options and try again."
    ))
    
    fit #returning linear model
  })
  
  
  # rval_fit() has NAs, we remove the option to trend or robust in eBayes to prevent failure
  observeEvent(input$button_limma_calculatemodel, {
    if (any(vapply(rval_fit(), function(x) {
      any(is.na(unlist(x)) |
          unlist(x) == Inf | unlist(x) == -Inf)
    }, logical(1)))) {
      message("NAs or Inf values detected, trend and robust options are disabled.")
      
      updateSwitchInput(session,
                        "select_limma_trend",
                        value = FALSE,
                        disabled = TRUE)
      
      updateSwitchInput(session,
                        "select_limma_robust",
                        value = FALSE,
                        disabled = TRUE)
    }
    
    else {
      message("NAs or Inf values not detected, trend and robust options are enabled")
      
      updateSwitchInput(session,
                        "select_limma_trend",
                        value = FALSE,
                        disabled = FALSE)
      
      updateSwitchInput(session,
                        "select_limma_robust",
                        value = FALSE,
                        disabled = FALSE)
      
    }
    
    # Creating calculate differences button
    rval_generated_limma_model(TRUE)
  })
  
  
  output$button_limma_calculatedifs_container = renderUI({
    if (rval_generated_limma_model())
      return(tagList(
        br(),
        
        h4("Contrasts options"),
        
        switchInput(
          inputId = "select_limma_trend",
          label = "eBayes Trend",
          labelWidth = "80px",
          value = FALSE
        ),
        
        switchInput(
          inputId = "select_limma_robust",
          label = "eBayes Robust",
          labelWidth = "80px",
          value = FALSE
        ),
        
        actionButton("button_limma_calculatedifs", "Calc. Contrasts")
      ))
    
    else
      return()
  })
  
  
  #render of plots and tables
  
  rval_plot_plotSA = reactive(create_plotSA(rval_fit()))
  output$graph_limma_plotSA = renderPlot(rval_plot_plotSA())
  output$table_limma_design =  DT::renderDT(
    rval_design(),
    rownames = TRUE,
    selection = "single",
    style = "bootstrap",
    options = list(
      autoWidth = TRUE,
      scrollX = TRUE,
      lengthChange = FALSE,
      searching=FALSE)
    )
  

  #Calculation of global difs
  rval_globaldifs = eventReactive(input$button_limma_calculatedifs, {
    try({
      globaldifs = calculate_global_difs(rval_gset_getBeta(), rval_voi(), rval_contrasts(), cores =
                                           n_cores)
    })
    
    if (!exists("globaldifs", inherits = FALSE)) {
      showModal(
        modalDialog(
          title = "Global difs calculation error",
          "An unexpected error has ocurred during global diffs and means calculation. Please, check your selected samples.",
          easyClose = TRUE,
          footer = NULL
        )
      )
    }
    
    validate(
      need(
        exists("globaldifs", inherits = FALSE),
        "An unexpected error has ocurred during global diffs and means calculation."
      )
    )
    
    globaldifs
    
    
  })
  
  #Calculation of differences (eBayes)
  rval_finddifcpgs = eventReactive(input$button_limma_calculatedifs, {
    try({
      dif_cpgs = find_dif_cpgs(
        design = rval_design(),
        fit = rval_fit(),
        contrasts = rval_contrasts(),
        trend = as.logical(input$select_limma_trend),
        robust = as.logical(input$select_limma_robust),
        cores = n_cores
      )
    })
    
    if (!exists("dif_cpgs", inherits = FALSE)) {
      showModal(
        modalDialog(
          title = "Contrasts Calculation Error",
          "An unexpected error has ocurred during contrasts calculation. Please, generate the model again and check if it is correct.
        If the problem persists, report the error to the maintainer",
          easyClose = TRUE,
          footer = NULL
        )
      )
      
      rval_generated_limma_model(FALSE)
      rval_analysis_finished(FALSE)
      shinyjs::disable("button_limma_heatmapcalc")
      
    }
    
    validate(
      need(
        exists("dif_cpgs", inherits = FALSE),
        "An unexpected error has ocurred during contrasts calculation. Please, generate the model again and check if it is correct.
        If the problem persists, report the error to the maintainer"
      )
    )
    
    rval_analysis_finished(TRUE)
    
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
    
    updateSelectInput(
      session,
      "select_limma_anncontrast",
      label = "Contrast",
      choices = rval_contrasts()
    )
    
    # disable button to avoid repeat clicking
    shinyjs::disable("button_limma_calculatedifs")
    
    #force rval_filteredlist
    rval_filteredlist()
    
    #enable or disable removebatch option
    covariables_design = as.matrix(rval_design()[, -seq_len(length(unique(rval_voi())))])
    
    if (ncol(covariables_design) > 0) {
      updateSwitchInput(session,
                        "select_limma_removebatch",
                        value = FALSE,
                        disabled = FALSE)
      updateSwitchInput(session,
                        "select_dmrs_removebatch",
                        value = FALSE,
                        disabled = FALSE)
    }
    else{
      updateSwitchInput(session,
                        "select_limma_removebatch",
                        value = FALSE,
                        disabled = TRUE)
      updateSwitchInput(session,
                        "select_dmrs_removebatch",
                        value = FALSE,
                        disabled = FALSE)
    }
    
    #enable and click heatmap button to get default graph
    shinyjs::enable("button_limma_heatmapcalc")
    shinyjs::click("button_limma_heatmapcalc")
    
    # enable again the button to allow repeat calculation
    shinyjs::enable("button_limma_calculatedifs") 
    
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
  
  rval_filteredlist2heatmap = reactive(
    {
      filtered_data = rval_filteredlist()[rval_contrasts() %in% input$select_limma_contrasts2plot] # filter contrasts2plot
      dif_cpgs = unique(data.table::rbindlist(filtered_data)$cpg)
      
      #If remove batch option is enabled, limma:removebatcheffects is applied using the design info data.
      if (!input$select_limma_removebatch)
      {
        join_table = rval_gset_getBeta()[dif_cpgs,]
      }
      else
      {
        voi_design = as.matrix(rval_design()[, seq_len(length(unique(rval_voi())))])
        covariables_design = as.matrix(rval_design()[, -seq_len(length(unique(rval_voi())))])
        
        join_table = as.data.frame(
          limma::removeBatchEffect(
            rval_gset_getBeta(),
            design = voi_design,
            covariates = covariables_design
          )[dif_cpgs,]
        )
      }
        
      join_table$cpg = NULL
      
      #If the number of CpGs is not in the plotting range, return NULL to avoid errors in plot_heatmap and disable download
      if (is.null(join_table) |
          nrow(join_table) < 2 | nrow(join_table) > 12000)
      {
        rval_filteredlist2heatmap_valid(FALSE)
        NULL
      }
      else
      {
        rval_filteredlist2heatmap_valid(TRUE)
        join_table
      }
      
  })
  
  rval_cpgcount_heatmap = eventReactive(input$button_limma_heatmapcalc, nrow(rval_filteredlist2heatmap()))
  
  rval_dendrogram = eventReactive(input$button_limma_heatmapcalc,
                                  {
                                    if (input$select_limma_rowsidecolors)
                                      
                                      #check if dendrogram cutting works (k should be minor than heatmap rows)
                                      try({
                                        dendrogram = create_dendrogram(
                                          rval_filteredlist2heatmap(),
                                          factorgroups =  factor(rval_voi()[rval_voi() %in% input$select_limma_groups2plot],
                                                                 levels = input$select_limma_groups2plot),
                                          groups2plot = rval_voi() %in% input$select_limma_groups2plot,
                                          clusteralg = input$select_limma_clusteralg,
                                          distance = input$select_limma_clusterdist,
                                          scale_selection = input$select_limma_scale,
                                          k_number = input$select_limma_knumber
                                        )
                                      })
                                    
                                    else
                                      dendrogram = NULL
                                    
                                    
                                    if (!exists("dendrogram", inherits = FALSE)) {
                                      dendrogram = NULL
                                      showModal(
                                        modalDialog(
                                          title = "Row clustering error",
                                          "An error has ocurred during cluster cutting from row dendrogram. Maybe the number of clusters selected is too high.",
                                          easyClose = TRUE,
                                          footer = NULL
                                        )
                                      )
                                    }
                                    
                                    dendrogram #returning the dendrogram classification
                                    
                                  })
  
  plot_heatmap = eventReactive(input$button_limma_heatmapcalc,
                               {
                                 validate(need(
                                   !is.null(rval_filteredlist2heatmap()),
                                   "Differences are not in the plotting range (<12000, >1)"
                                 ),
                                 need(input$select_limma_groups2plot != "", 
                                      "Select at least one group to plot."))
                                 
                                 create_heatmap(
                                   rval_filteredlist2heatmap(),
                                   factorgroups =  factor(rval_voi()[rval_voi() %in% input$select_limma_groups2plot],
                                                          levels = input$select_limma_groups2plot),
                                   groups2plot = rval_voi() %in% input$select_limma_groups2plot,
                                   Colv = as.logical(input$select_limma_colv),
                                   ColSideColors = input$select_limma_colsidecolors,
                                   RowSideColors = rval_dendrogram(),
                                   clusteralg = input$select_limma_clusteralg,
                                   distance = input$select_limma_clusterdist,
                                   scale = input$select_limma_scale,
                                   static = as.logical(input$select_limma_graphstatic)
                                 )
                               })
  
  make_table = eventReactive(input$button_limma_heatmapcalc,
                             {
                               default_df = data.frame(
                                 contrast = rval_contrasts(),
                                 Hypermethylated = 0,
                                 Hypomethylated = 0,
                                 total = 0
                               )
                               
                               count_df = data.table::rbindlist(rval_filteredlist(), idcol = "contrast")
                               
                               if (nrow(count_df) > 0 & !is.null(count_df)) {
                                 count_df = count_df %>%
                                   dplyr::mutate(type = factor(
                                     ifelse(
                                       .data$dif_current < 0,
                                       "Hypermethylated",
                                       "Hypomethylated"
                                     ),
                                     levels = c("Hypermethylated", "Hypomethylated")
                                   )) %>%
                                   dplyr::group_by(.data$contrast, .data$type)  %>%
                                   dplyr::summarise(CpGs = dplyr::n()) %>%
                                   tidyr::complete(.data$contrast, .data$type, fill = list(CpGs =
                                                                                             0)) %>%
                                   tidyr::pivot_wider(names_from = .data$type,
                                                      values_from = .data$CpGs) %>%
                                   dplyr::mutate(total = .data$Hypermethylated + .data$Hypomethylated) %>%
                                   dplyr::mutate(dplyr::across(
                                     c("Hypermethylated", "Hypomethylated", "total"),
                                     ~ format(., scientific = FALSE, digits = 0)
                                   ))
                               }
                               
                               rbind(as.data.frame(count_df), default_df[!(default_df[["contrast"]] %in% count_df[["contrast"]]), ])
                               
                             })
  
  
  
  observeEvent(input$button_limma_heatmapcalc,
               {
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
                       div(
                        min_width=400,
                         plotOutput(
                         "graph_limma_heatmap_static",
                         width = "600px",
                         height = "600px"
                       ) %>% shinycssloaders::withSpinner()
                     ))
                 })
                 
                 if (!as.logical(input$select_limma_graphstatic))
                   output$graph_limma_heatmap_interactive = plotly::renderPlotly(plot_heatmap())
                 else
                   output$graph_limma_heatmap_static = renderPlot(plot_heatmap())
               })
  
  output$text_limma_heatmapcount = renderText(paste("DMPs in heatmap:", rval_cpgcount_heatmap()))
  
  output$table_limma_difcpgs = renderTable(make_table(), digits = 0)
  
  table_annotation = eventReactive(input$button_limma_anncalc,{
    
    req(rval_filteredlist())
    
    dif_target = paste("dif",
                       limma::strsplit2(input$select_limma_anncontrast, "-")[1],
                       limma::strsplit2(input$select_limma_anncontrast, "-")[2],
                       sep = "_")
    
    temp = as.data.frame(minfi::getAnnotation(rval_gset()))
    
    int_cols = c(
      "Name",
      "chr",
      "pos",
      "strand",
      "Relation_to_Island",
      "UCSC_RefGene_Name",
      "UCSC_RefGene_Group"
    )
    
    temp = temp[row.names(temp) %in% rval_filteredlist()[[input$select_limma_anncontrast]]$cpg, int_cols]
    
    temp$dif_beta = rval_globaldifs()[[dif_target]][rval_globaldifs()[["cpg"]] %in% row.names(temp)]
    temp$fdr = rval_filteredlist()[[input$select_limma_anncontrast]][["adj.P.Val"]][rval_filteredlist()[[input$select_limma_anncontrast]][["cpg"]] %in% row.names(temp)]

    temp
  })
  
  output$table_limma_ann = DT::renderDT(
    table_annotation(),
    extensions = 'Buttons',
    rownames = FALSE,
    selection = "single",
    style = "bootstrap",
    caption = "DMPs Annotation:",
    options = list(
      dom = 'Blfrtip',
      lengthMenu = list(c(10,25,50,100,25000), c(10,25,50,100,25000)),
      pageLength = 10,
      autoWidth = TRUE,
      scrollX = TRUE,
      buttons = c('csv', 'excel', 'print')
    ),
  )
  
  
  #Disable or enable buttons depending on software state
  observeEvent(rval_analysis_finished(),
               if (rval_analysis_finished())
               {
                 shinyjs::enable("download_export_robjects")
                 shinyjs::enable("download_export_filteredbeds")
                 shinyjs::enable("download_export_markdown")
                 shinyjs::enable("button_dmrs_calculate")
                 
                 updatePickerInput(
                   session,
                   "select_dmrs_contrasts",
                   selected = rval_contrasts(),
                   choices = rval_contrasts(),
                 )
                 updatePickerInput(
                   session,
                   "select_dmrs_platform",
                   selected = if(nrow(rval_finddifcpgs()[[1]])>450001) {"EPIC"} else {"450k"},
                   choices = c("450k", "EPIC")
                 )
                 
               }
               
               else
               {
                 shinyjs::disable("download_export_robjects")
                 shinyjs::disable("download_export_filteredbeds")
                 shinyjs::disable("download_export_markdown")
                 shinyjs::disable("download_export_heatmaps")
                 shinyjs::disable("button_dmrs_calculate")
               })
  
  
  #DMRs
  
  rval_mcsea = eventReactive(input$button_dmrs_calculate, {
    
    validate(
      need(
        rval_analysis_finished(),
        "To calculate DMRs, you have to finish first the DMP calculation."
      )
    )
    
    validate(
      need(
        requireNamespace("mCSEA", quietly = T),
        "mCSEA is not installed. You should install the package to calculate DMRs"
      ),
      need(
        input$select_dmrs_contrasts != "",
        "You should select at least one contrast."
      ),
      need(
        input$select_dmrs_regions != "",
        "You should select at least one DMR type."
      )
    )
    print("calculating rval_mcsea")
    
    updateSelectInput(session,
                      "select_dmrs_selcont",
                      choices = input$select_dmrs_contrasts)
    updateSelectInput(session,
                      "select_dmrs_selreg",
                      choices = input$select_dmrs_regions)
    
    updateSelectInput(
      session,
      "select_dmrs_groups2plot",
      label = "Groups to plot",
      choices = levels(rval_voi()),
      selected = levels(rval_voi())
    )
    updateSelectInput(
      session,
      "select_dmrs_contrasts2plot",
      label = "Contrasts to plot",
      choices = rval_contrasts(),
      selected = rval_contrasts()
    )
    updateSelectInput(
      session,
      "select_dmrs_regions2plot",
      label = "Regions to plot",
      choices = input$select_dmrs_regions,
      selected = input$select_dmrs_regions
    )
    
    shinyjs::disable("button_dmrs_calculate")
    
    withProgress(message = "Calculating DMRs...",
                 max = 3,
                 value = 1,
                 
                 {
                   
                   try({
                     dmrs_result = find_dmrs(
                       rval_finddifcpgs(),
                       minCpGs = input$slider_dmrs_cpgs,
                       platform = "EPIC",
                       voi = rval_voi(),
                       regionsTypes = input$select_dmrs_regions,
                       contrasts = input$select_dmrs_contrasts,
                       bvalues = rval_gset_getBeta(),
                       permutations = input$slider_dmrs_permutations,
                       ncores = n_cores
                     )
                     
                     setProgress(value = 2, message = "Calculating differential of betas...")
                     
                     dmrs_result = add_dmrs_globaldifs(
                       mcsea_result = dmrs_result,
                       cpg_globaldifs = rval_globaldifs(),
                       regionsTypes = input$select_dmrs_regions,
                       cores = n_cores
                     )
                   })
                 })
    
    shinyjs::enable("button_dmrs_calculate")
    
    if (!exists("dmrs_result", inherits = FALSE)) {
      rval_dmrs_finished(FALSE)
      showModal(
        modalDialog(
          title = "DMR calculation error",
          "An unexpected error has occurred during DMRs calculation.",
          easyClose = TRUE,
          footer = NULL
        )
      )
      shinyjs::disable("button_dmrs_heatmapcalc")
    }
    
    validate(need(
      exists("dmrs_result", inherits = FALSE),
      "An unexpected error has occurred during DMRs calculation."
    ))
    
    #enable heatmap button
    shinyjs::enable("button_dmrs_heatmapcalc")
    rval_dmrs_finished(TRUE)
    rval_dmrs_ready2heatmap(TRUE)

    dmrs_result
  })
  
  rval_filteredmcsea = eventReactive(list(input$button_dmrs_calculate, input$button_dmrs_heatmapcalc),{
    
    req(rval_mcsea())
    print("calculating rval_filteredmcsea")
    
    filter_dmrs(
      mcsea_list = rval_mcsea(),
      fdr = input$slider_dmrs_adjpvalue,
      pval = input$slider_dmrs_pvalue,
      dif_beta = input$slider_dmrs_deltab,
      regionsTypes = input$select_dmrs_regions,
      contrasts = input$select_dmrs_contrasts
    )
    
  })
  
  rval_filteredmcsea2heatmap = reactive({
    
    req(input$select_dmrs_regions2plot)
    req(input$select_dmrs_contrasts2plot)
    
    print(input$select_dmrs_regions2plot)
    print(input$select_dmrs_contrasts2plot)
    print(str(rval_filteredmcsea()))
    
    voi_design = as.matrix(rval_design()[, seq_len(length(unique(rval_voi())))])
    covariables_design = as.matrix(rval_design()[,-seq_len(length(unique(rval_voi())))])

    #filtering contrasts
    join_table = create_dmrs_heatdata(
      mcsea_result = rval_filteredmcsea()[input$select_dmrs_contrasts2plot],
      bvalues = if (!input$select_dmrs_removebatch) {
        rval_gset_getBeta()
      } else {
        as.data.frame(
          limma::removeBatchEffect(
            rval_gset_getBeta(),
            design = voi_design,
            covariates = covariables_design
          )
        )
      },
      regions = input$select_dmrs_regions2plot,
      contrasts = input$select_dmrs_contrasts
    )
    
    #If the number of CpGs is not in the plotting range, return NULL to avoid errors in plot_dmrsheatmap
    if (is.null(join_table) |
        nrow(join_table) < 2 | nrow(join_table) > 12000)
    {
      rval_filteredmcsea2heatmap_valid(FALSE)
      NULL
    }
    else
    {
      rval_filteredmcsea2heatmap_valid(TRUE)
      join_table
    }
    
  })
  
  observeEvent(rval_dmrs_ready2heatmap(),{
    if(rval_dmrs_ready2heatmap()){
      rval_dmrs_ready2heatmap(FALSE)
      print("clicked")
      
      shinyjs::click("button_dmrs_heatmapcalc")
    }
  })
  
  #Heatmap DMRs
  
  plot_dmrsheatmap = eventReactive(input$button_dmrs_heatmapcalc, ignoreNULL=FALSE, {

    validate(need(input$fileinput_input != "", "DMR calculation has not been performed or data has not been uploaded."))
    
    validate(need(
      !is.null(rval_filteredmcsea2heatmap()),
      "Differences are not in the plotting range (<12000, >1)"
    ))
    
    create_heatmap(
      rval_filteredmcsea2heatmap(),
      factorgroups =  factor(rval_voi()[rval_voi() %in% input$select_dmrs_groups2plot],
                             levels = input$select_dmrs_groups2plot),
      groups2plot = rval_voi() %in% input$select_dmrs_groups2plot,
      Colv = as.logical(input$select_dmrs_colv),
      ColSideColors = input$select_dmrs_colsidecolors,
      RowSideColors = rval_dendrogram_dmrs(),
      clusteralg = input$select_dmrs_clusteralg,
      distance = input$select_dmrs_clusterdist,
      scale = input$select_dmrs_scale,
      static = as.logical(input$select_dmrs_graphstatic)
    )
  })
  
  rval_dendrogram_dmrs = eventReactive(input$button_dmrs_heatmapcalc,
                                  {
                                    if (input$select_dmrs_rowsidecolors)
                                      
                                      #check if dendrogram cutting works (k should be minor than heatmap rows)
                                      try({
                                        dendrogram = create_dendrogram(
                                          rval_filteredmcsea2heatmap(),
                                          factorgroups =  factor(rval_voi()[rval_voi() %in% input$select_dmrs_groups2plot],
                                                                 levels = input$select_dmrs_groups2plot),
                                          groups2plot = rval_voi() %in% input$select_dmrs_groups2plot,
                                          clusteralg = input$select_dmrs_clusteralg,
                                          distance = input$select_dmrs_clusterdist,
                                          scale_selection = input$select_dmrs_scale,
                                          k_number = input$select_dmrs_knumber
                                        )
                                      })
                                    
                                    else
                                      dendrogram = NULL
                                    
                                    
                                    if (!exists("dendrogram", inherits = FALSE)) {
                                      dendrogram = NULL
                                      showModal(
                                        modalDialog(
                                          title = "Row clustering error",
                                          "An error has ocurred during cluster cutting from row dendrogram. Maybe the number of clusters selected is too high.",
                                          easyClose = TRUE,
                                          footer = NULL
                                        )
                                      )
                                    }
                                    
                                    dendrogram #returning the dendrogram classification
                                    
                                  })
  rval_cpgcount_dmrs_heatmap = eventReactive(input$button_dmrs_heatmapcalc, nrow(rval_filteredmcsea2heatmap()))
  
  output$graph_dmrs_heatmap = renderPlot(plot_dmrsheatmap())
  
  output$text_dmrs_heatmapcount = renderText(paste("DMRs in heatmap:", rval_cpgcount_dmrs_heatmap()))
  
  
  #DMRs count table
  
  make_table_dmrscount = eventReactive(rval_filteredmcsea(),{
    print("calculating table")
    result = data.frame(
      contrast =  rep(
        names(rval_filteredmcsea()),
        each = length(input$select_dmrs_regions)
      ),
      region_type = input$select_dmrs_regions
    )
    
    
    result$Hypermethylated = apply(result, 1, function(x) {
      nrow(rval_filteredmcsea()[[x[1]]][[x[2]]][rval_filteredmcsea()[[x[1]]][[x[2]]]$NES < 0, ])
    })
    
    result$Hypomethylated = apply(result, 1, function(x) {
      nrow(rval_filteredmcsea()[[x[1]]][[x[2]]][rval_filteredmcsea()[[x[1]]][[x[2]]]$NES > 0, ])
    })
    
    result$total = result$Hypermethylated + result$Hypomethylated
    
    result
  })
  
  output$table_dmrs_count = renderTable(make_table_dmrscount(), digits = 0)
  
  #DMRs single plot
  

  plot_singledmr = eventReactive(input$button_dmrs_graphsingle,
                                 {
                                   validate(need(!is.null(input$table_dmrs_table_rows_selected), "A DMR should be selected."))
                                   
                                   selected_dmr = row.names(rval_filteredmcsea()[[input$select_dmrs_selcont]][[input$select_dmrs_selreg]])[input$table_dmrs_table_rows_selected]
                                   
                                   mCSEA::mCSEAPlot(
                                     rval_mcsea()[[input$select_dmrs_selcont]],
                                     regionType = input$select_dmrs_selreg,
                                     dmrName = selected_dmr,
                                     makePDF = F,
                                     transcriptAnnotation = "symbol"
                                   )
                                 })
  

  
  output$graph_dmrs_singledmr = renderPlot(plot_singledmr())
  
  
  rval_table_sigdmrs = reactive({
    rval_filteredmcsea()[[input$select_dmrs_selcont]][[input$select_dmrs_selreg]]
  })
  
  output$table_dmrs_table = DT::renderDT(
    rval_table_sigdmrs(),
    rownames = TRUE,
    extensions = 'Buttons',
    selection = "single",
    style = "bootstrap",
    caption = "Select DMR to plot:",
    options = list(
      pageLength = 10,
      dom = 'Blfrtip',
      lengthMenu = list(c(10,25,50,100,25000), c(10,25,50,100,25000)),
      autoWidth = TRUE,
      scrollX = TRUE,
      buttons = c('csv', 'excel', 'print'),
      columnDefs = list(
        list(
          targets = match("leadingEdge", colnames(rval_filteredmcsea()[[input$select_dmrs_selcont]][[input$select_dmrs_selreg]])),
          visible = FALSE
        ))
    )
  )
  
  
  
  
  #EXPORT
  
  #R Objects
  
  output$download_export_robjects = downloadHandler(
    "R_Objects.zip",
    
    
    content = function(file) {
      shinyjs::disable("download_export_robjects")
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
                     global_difs = rval_globaldifs()
                     
                     objetos = list(RGSet="rgset",GenomicRatioSet="gset",fit="fit",design="design",
                                    ebayestables="ebayes_tables",Bvalues="bvalues",Mvalues="mvalues",global_difs="global_difs")
                     
                     #if(rval_dmrs_finished()) mcsea_result = rval_mcsea()
                     
                     setProgress(value = 2, message = "Saving RObjects...")
                     
                     for(id in input$select_export_objects2download){
                       do.call(function(object, name) save(object,file=paste0(tempdir(),"/",name,".RData")), list(object=as.name(objetos[[id]]), name=objetos[[id]]))
                     }
                     
                     objects2save = list.files(tempdir(), pattern = "*.RData", full.names = TRUE)
                     
                     setProgress(value = 3, message = "Compressing RObjects...")
                     zip::zip(
                       file,
                       objects2save,
                       recurse = FALSE,
                       include_directories = FALSE,
                       mode = "cherry-pick"
                     )
                     shinyjs::enable("download_export_robjects")
                   })
    }
  )
  
  rval_annotation = reactive({
    
      if (input$select_export_genometype == "hg38" &
          (!requireNamespace("rtracklayer", quietly = T)) |
          !requireNamespace("GenomicRanges", quietly = T)){
        
        showModal(
          modalDialog(
            title = "rtracklayer::liftOver not available",
            "Rtracklayer is not installed and it is needed to liftOver annotation from hg19 to hg38 genome. Please, install the package and restart the R session.",
            easyClose = TRUE,
            footer = NULL
          )
        )
        updateSelectInput(session, "select_export_genometype", choices = "hg19", selected = "hg19")
      }
        
    annotation = as.data.frame(minfi::getAnnotation(rval_gset()))
  
    if(input$select_export_genometype == "hg19"){
      data.frame(row.names = row.names(annotation), chr = annotation$chr, pos = annotation$pos)
    }
    else{
      chain = rtracklayer::import.chain(system.file("hg19ToHg38.over.chain", package = "shinyepico"))
      ann_granges = data.frame(chr=annotation$chr, start=annotation$pos - 1, end = annotation$pos, name = row.names(annotation))
      ann_granges = GenomicRanges::makeGRangesFromDataFrame(ann_granges, starts.in.df.are.0based=TRUE, keep.extra.columns = T)
      ann_granges = unlist(rtracklayer::liftOver(ann_granges, chain = chain))
      
      data.frame(row.names = GenomicRanges::mcols(ann_granges)[[1]], chr=GenomicRanges::seqnames(ann_granges), pos = GenomicRanges::start(ann_granges))
    }
      
    
  }
  )

  #Bed downloading enabling/disabling control
  observe({
    if (input$select_export_analysistype == "DMPs") {
      if (rval_analysis_finished() &
          input$select_export_bedtype == "by heatmap cluster" &
          (!rval_filteredlist2heatmap_valid() |
           is.null(rval_dendrogram()))) {
        shinyjs::disable("download_export_filteredbeds")
      } else if (rval_analysis_finished() &
                 rval_filteredlist2heatmap_valid() &
                 !is.null(rval_dendrogram())) {
        shinyjs::enable("download_export_filteredbeds")
      } else if (rval_analysis_finished() &
                 input$select_export_bedtype == "by contrasts") {
        shinyjs::enable("download_export_filteredbeds")
      }
    }
    else{
      if (rval_dmrs_finished() &
          input$select_export_bedtype == "by heatmap cluster" &
          (!rval_filteredmcsea2heatmap_valid() |
           is.null(rval_dendrogram_dmrs()))) {
        shinyjs::disable("download_export_filteredbeds")
      }
      else if (rval_dmrs_finished() &
               rval_filteredmcsea2heatmap_valid() &
               !is.null(rval_dendrogram_dmrs())) {
        shinyjs::enable("download_export_filteredbeds")
      } else if (rval_dmrs_finished() &
                 input$select_export_bedtype == "by contrasts") {
        shinyjs::enable("download_export_filteredbeds")
      }
    }
  })
  
  #heatmap downloading enable/disable
  observe({
    if (input$select_export_heatmaptype == "DMPs") {
      if (rval_analysis_finished() & rval_filteredlist2heatmap_valid())
        shinyjs::enable("download_export_heatmaps")
      else
        shinyjs::disable("download_export_heatmaps")
    }
    else{
      if (rval_dmrs_finished() & rval_filteredmcsea2heatmap_valid())
        shinyjs::enable("download_export_heatmaps")
      else
        shinyjs::disable("download_export_heatmaps")
    }
  })

  
  #Filtered BEDs
  output$download_export_filteredbeds  = downloadHandler(
    "filtered_beds.zip",
    
    content = function(file) {
      shinyjs::disable("download_export_filteredbeds")
      withProgress(message = "Generating BED files...",
                   value = 1,
                   max = 4,
                   {
                     
                     if (input$select_export_bedtype == "by heatmap cluster" & input$select_export_analysistype == "DMPs") {
                       create_filtered_bed_clusters( dendro_data = rval_dendrogram(),
                                                     annotation = rval_annotation(),
                                                     directory = dirname(file)
                       )
                       
                     }
                     else if(input$select_export_bedtype == "by contrasts" & input$select_export_analysistype == "DMPs"){
                       create_filtered_beds(rval_filteredlist(),
                                            rval_annotation(),
                                            directory = dirname(file))
                     }
                     else if(input$select_export_bedtype == "by heatmap cluster" & input$select_export_analysistype == "DMRs"){
                       create_filtered_bed_dmrs_clusters(
                         dendro_data = rval_dendrogram_dmrs(),
                         mcsea_filtered = rval_filteredmcsea(),
                         annotation = rval_annotation(),
                         directory = dirname(file)
                       )
                       create_dmrs_bed_background(
                         mcsea_result = rval_mcsea(),
                         collapse = TRUE,
                         regionsTypes = input$select_dmrs_regions2plot,
                         annotation = rval_annotation(),
                         directory = dirname(file)
                       )
                     }
                     else if(input$select_export_bedtype == "by contrasts" & input$select_export_analysistype == "DMRs"){
                       create_filtered_beds_dmrs(
                         mcsea_filtered = rval_filteredmcsea(),
                         regionsTypes = input$select_dmrs_regions,
                         annotation = rval_annotation(),
                         directory = dirname(file)
                       )
                       create_dmrs_bed_background(
                         mcsea_result = rval_mcsea(),
                         collapse = FALSE,
                         regionsTypes = input$select_dmrs_regions,
                         annotation = rval_annotation(),
                         directory = dirname(file)
                       )

                     }
                     
                     
                     objects = list.files(path = dirname(file),
                                          full.names = TRUE,
                                          pattern = "*.bed")
                     
                     zip::zip(
                       file,
                       objects,
                       recurse = FALSE,
                       include_directories = FALSE,
                       mode = "cherry-pick"
                     )
                     
                     file.remove(objects) #removing objects after exporting
                     
                     shinyjs::enable("download_export_filteredbeds")
                     
                   })
    }
  )
  
  
  #Markdown Report
  output$download_export_markdown = downloadHandler(
    filename = "Report.html",
    content = function(file) {
      shinyjs::disable("download_export_markdown")
      withProgress(message = "Generating Report...",
                   value = 1,
                   max = 3,
                   {
                     params = list(
                       rval_sheet = rval_sheet(),
                       rval_sheet_target = rval_sheet_target(),
                       name_var = input$select_input_samplenamevar,
                       grouping_var = input$select_input_groupingvar,
                       donor_var = input$select_input_donorvar,
                       normalization_mode = input$select_minfi_norm,
                       dropsnps = input$select_minfi_dropsnps,
                       dropcphs = input$select_minfi_dropcphs,
                       dropsex = input$select_minfi_chromosomes,
                       maf = input$slider_minfi_maf,
                       probes = rval_gsetprobes(),
                       limma_voi = input$select_limma_voi,
                       limma_covar = input$checkbox_limma_covariables,
                       limma_inter = input$checkbox_limma_interactions,
                       limma_arrayweights = input$select_limma_weights,
                       limma_ebayes_trend = input$select_limma_trend,
                       limma_ebayes_robust = input$select_limma_robust,
                       rval_design = rval_design(),
                       rval_contrasts = rval_contrasts(),
                       rval_voi = rval_voi(),
                       rval_dendrogram = rval_dendrogram(),
                       min_deltabeta = input$slider_limma_deltab,
                       max_fdr = input$slider_limma_adjpvalue,
                       max_pvalue = input$slider_limma_pvalue,
                       clusteralg = input$select_limma_clusteralg,
                       groups2plot = input$select_limma_groups2plot,
                       contrasts2plot = input$select_limma_contrasts2plot,
                       Colv = input$select_limma_colv,
                       distance = input$select_limma_clusterdist,
                       scale = input$select_limma_scale,
                       removebatch = input$select_limma_removebatch,
                       plot_densityplotraw = rval_plot_densityplotraw(),
                       plot_densityplot = rval_plot_densityplot(),
                       plot_pcaplot = rval_plot_pca()[["graph"]],
                       plot_corrplot = rval_plot_corrplot()[["graph"]],
                       plot_boxplotraw = rval_plot_boxplotraw(),
                       plot_boxplot = rval_plot_boxplot(),
                       plot_qcraw = rval_plot_qcraw(),
                       plot_bisulfiterawII = rval_plot_bisulfiterawII(),
                       plot_sexprediction = rval_plot_sexprediction(),
                       plot_snpheatmap = rval_plot_snpheatmap(),
                       plot_plotSA = rval_plot_plotSA(),
                       table_pcaplot = rval_plot_pca()[["info"]],
                       table_corrplot = rval_plot_corrplot()[["info"]],
                       data_sexprediction = as.data.frame(minfi::pData(rval_gset()))[["predictedSex"]],
                       table_dmps = make_table(),
                       filteredlist2heatmap = rval_filteredlist2heatmap()
                     )
                     
                     #if DMR analysis has been done, we add specific parameters
                     if (rval_dmrs_finished()) {
                       params[["dmrs_contrasts"]] = input$select_dmrs_contrasts
                       params[["dmrs_rval_dendrogram"]] = rval_dendrogram_dmrs()
                       params[["dmrs_min_deltabeta"]] =  input$slider_dmrs_deltab
                       params[["dmrs_max_fdr"]] = input$slider_dmrs_adjpvalue
                       params[["dmrs_max_pvalue"]] = input$slider_dmrs_pvalue
                       params[["dmrs_clusteralg"]] = input$select_dmrs_clusteralg
                       params[["dmrs_groups2plot"]] = input$select_dmrs_groups2plot
                       params[["dmrs_contrasts2plot"]] = input$select_dmrs_contrasts2plot
                       params[["dmrs_regions2plot"]] = input$select_dmrs_regions2plot
                       params[["dmrs_Colv"]] = input$select_dmrs_colv
                       params[["dmrs_distance"]] = input$select_dmrs_clusterdist
                       params[["dmrs_scale"]] = input$select_dmrs_scale
                       params[["dmrs_removebatch"]] = input$select_dmrs_removebatch
                       params[["table_dmrs"]] = make_table_dmrscount()
                       params[["filteredmcsea2heatmap"]] = rval_filteredmcsea2heatmap()
                       
                     }
                     
                     newenv = new.env(parent = globalenv())
                     newenv$create_heatmap = create_heatmap
                     
                     render_file = rmarkdown::render(
                       input = system.file("report.Rmd", package = "shinyepico"),
                       output_file = file,
                       run_pandoc = TRUE,
                       params = params,
                       envir = newenv
                     )
                     
                     shinyjs::enable("download_export_markdown")
                     
                   })
    }
  )
  
  #custom heatmap
  output$download_export_heatmaps = downloadHandler(
    "Custom_heatmap.pdf",
    content = function(file) {
      shinyjs::disable("download_export_heatmaps")
      withProgress(message = "Generating plot...",
                   value = 1,
                   max = 4,
                   {
                     grDevices::pdf(file = file,
                                    height = 11.71,
                                    width = 8.6)
                     if(input$select_export_heatmaptype == "DMPs"){
                       create_heatmap(
                         rval_filteredlist2heatmap(),
                         factorgroups =  factor(rval_voi()[rval_voi() %in% input$select_limma_groups2plot],
                                                levels = input$select_limma_groups2plot),
                         groups2plot = rval_voi() %in% input$select_limma_groups2plot,
                         Colv = as.logical(input$select_limma_colv),
                         ColSideColors = as.logical(input$select_limma_colsidecolors),
                         RowSideColors = rval_dendrogram(),
                         clusteralg = input$select_limma_clusteralg,
                         distance = input$select_limma_clusterdist,
                         scale = input$select_limma_scale,
                         static = TRUE
                       )
                     }
                     else{
                       create_heatmap(
                         rval_filteredmcsea2heatmap(),
                         factorgroups = factor(rval_voi()[rval_voi() %in% input$select_dmrs_groups2plot],
                                                levels = input$select_dmrs_groups2plot),
                         groups2plot = rval_voi() %in% input$select_dmrs_groups2plot,
                         Colv = as.logical(input$select_dmrs_colv),
                         ColSideColors = as.logical(input$select_dmrs_colsidecolors),
                         RowSideColors = rval_dendrogram_dmrs(),
                         clusteralg = input$select_dmrs_clusteralg,
                         distance = input$select_dmrs_clusterdist,
                         scale = input$select_dmrs_scale,
                         static = TRUE
                       )
                     }

                     grDevices::dev.off()
                     shinyjs::enable("download_export_heatmaps")
                   })
    }
  )
}
