# project ÉPICo! 
#Octavio Morante-Palacios

#### packages

library(minfi)
library(IlluminaHumanMethylationEPICmanifest)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(limma)
library(shiny)
library(heatmaply)
library(doParallel)
library(shinythemes)
library(data.table)
library(shinycssloaders)
library(dplyr)
library(tidyr)
library(statmod)
library(ggplot2)
library(gplots)

# ENVIRONMENT VARIABLES
norm_options = c("Raw","Illumina", "Funnorm", "Noob", "Quantile", "Noob+Quantile")
hclust_methods = c("single", "complete", "average", "mcquitty", "median","centroid")
options(shiny.maxRequestSize=500*1024^2)

#PERFORMANCE SETTINGS
if (detectCores() > 1) {cores = round(detectCores()/2, digits=0) } else {cores = 1}


#FUNCTIONS
source("utils.R", local=TRUE)
####

####### UI ################ 

ui <- navbarPage("ÉPICo!", id="navbar_epic", theme=shinytheme("sandstone") ,
    tabPanel("Input",
          titlePanel("Input sheet file"),
          sidebarLayout(
              sidebarPanel(
                  fileInput(inputId="fileinput_input", "Upload Compress Experiment Directory (zip):",multiple = FALSE,accept = "application/zip"),
                  #directoryInput("sample_directory", "Select Sample Directory:", value="/home/octavio/Documentos/Análisis DCDEX/Pruebas limma metilación/Arrays de metilación/"),
                  uiOutput("ui_button_input_load"),
                  h3(),
                  
                  conditionalPanel("typeof output.samples_table != 'undefined'", 
                                   selectInput("select_input_samplenamevar", "", choices=c()),
                                   h3(),
                                   checkboxGroupInput("selected_samples", "", c(), selected=TRUE ),
                                   h3(),
                                   selectInput("select_input_groupingvar", "", c()),
                                   h3(),
                                   selectInput("select_input_donorvar", "", c()),
                                   h3(),
                                   actionButton("button_input_next", "Proceed to the next step"))
                          ),
              mainPanel(DT::DTOutput("samples_table") %>% withSpinner())
          )
                                    
     ), 
    
    
    tabPanel("Minfi Norm.",
          titlePanel("Minfi Normalization"),
          sidebarLayout(
              sidebarPanel(
                  selectInput("select_minfi_norm", "Select Normalization", norm_options),
                  #actionButton("button_minfi_check", "Check"),
                  actionButton("button_minfi_select", "Select")),
              mainPanel(tabsetPanel(
                
                tabPanel("Density plot",
                    h4("Raw"),
                    plotOutput("graph_minfi_densityplotraw")%>% withSpinner(),
                    h4("Processed"),
                    plotOutput("graph_minfi_densityplot") %>% withSpinner()),
                
                tabPanel("Density Bean Plot", 
                    h4("Raw"),
                    plotOutput("graph_minfi_densitybeanplotraw")%>% withSpinner(),
                    h4("Processed"),
                    plotOutput("graph_minfi_densitybeanplot") %>% withSpinner()),
                
                tabPanel("MDS plot", 
                    h4("Raw"),
                    plotOutput("graph_minfi_mdsplotraw")%>% withSpinner(),
                    h4("Processed"),
                    plotOutput("graph_minfi_mdsplot") %>% withSpinner()),
                
                tabPanel("Boxplot", 
                         h4("Raw"),
                         plotOutput("graph_minfi_boxplotraw")%>% withSpinner(),
                         h4("Processed"),
                         plotOutput("graph_minfi_boxplot") %>% withSpinner()),
                
                tabPanel("SNPs Heatmap", h4("SNPs beta-values (Raw)"),plotlyOutput("graph_minfi_snps") %>% withSpinner()),
                
                tabPanel("Sex prediction", h4("plotSex"), plotOutput("graph_minfi_sex") %>% withSpinner()),
                
                tabPanel("QC plot", 
                         h4("Raw"),
                         plotOutput("graph_minfi_qcraw")%>% withSpinner(),
                         h4("Processed"),
                         plotOutput("graph_minfi_qc") %>% withSpinner())
              
              
                ))
           )),
    
    tabPanel("Limma",
             titlePanel("Limma parameters"),
             sidebarLayout(
               sidebarPanel(width=3,
                 selectInput("select_limma_voi", "Select Variable of Interest", c()),
                 checkboxGroupInput("checkbox_limma_covariables", "Select covariables to block", c()),
                 #checkboxGroupInput("checkbox_limma_groups", "Select groups to compare", c()),
                 selectInput("select_limma_trend","eBayes trend option:", c("TRUE","FALSE"), c("FALSE")),
                 selectInput("select_limma_robust","eBayes robust option:", c("TRUE","FALSE"), c("FALSE")),
                 actionButton("button_limma_calculatemodel", "Generate Model"),
                 actionButton("button_limma_calculatedifs", "Calculate Contrasts")
                 
               ),
               mainPanel(width=9,
                         tabsetPanel(id="tabset_limma",
                          tabPanel("Model diagnosis", value="model_diagnosis",
                                   plotOutput("graph_limma_plotSA") %>% withSpinner(), 
                                   plotOutput("graph_limma_plotMA")), 
                          tabPanel("Differential CpGs", value="differential_cpgs",
                                  div(style = 'max-width:800px;margin:auto;',fluidPage(
                                    h4("Heatmap"),
                                    uiOutput("graph_limma_heatmapcontainer"),
                                    h4("CpGs selected with these filters:"),
                                    tableOutput("table_limma_difcpgs") %>% withSpinner(),
                                    column(6,
                                           h4("Group options"),
                                           selectizeInput("select_limma_groups2plot", "Groups to plot", c(), multiple = TRUE), 
                                           selectizeInput("select_limma_contrasts2plot", "Contrasts to plot", c(), multiple=TRUE)
                                    ),
                                    column(6, 
                                           h4("Filter options"),
                                           sliderInput("slider_limma_deltab", "Min. DeltaBeta", 0, 1, 0.2), 
                                           sliderInput("slider_limma_adjpvalue", "Max. FDR", 0,1,0.05)
                                    ),
                                    h4("Clustering options", align="center" ),
                                    column(6,
                                           selectInput("select_limma_clusteralg", "Clustering algorithm:", c("single", "complete", "average", "mcquitty", "median","centroid"), "average"),
                                           selectInput("select_limma_clusterdist", "Distance Function:", c("pearson", "spearman", "kendall", "euclidean"), "pearson"),
                                           selectInput("select_limma_graphstatic", label = "Plot Static Graph:", c(TRUE,FALSE),c(FALSE))
                                           
                                           ),
                                    
                                    column(6,
                                           selectInput("select_limma_colv", "Create Column Dendogram:", c(TRUE, FALSE), c(TRUE)),
                                           selectInput("select_limma_scale", "Scale:", c("row","none"), "row"),
                                           tags$br(),
                                           actionButton("button_limma_heatmapcalc", "Update")
                                           ),
                                    tags$br(),tags$br()
                                    
                                   
                                  )))
                         )
             ))
          ),
    
    tabPanel("Export", 
             
             h3("Download RObjects:"),
             downloadButton("download_export_robjects"),
             p("Press to download the R objects used for the analysis (RGSet, GenomicRatioSet, Bvalues, Mvalues, etc."),
             h3("Download filtered bed files:"),
             downloadButton("download_export_filteredbeds"),
             p("Press to download the created filtered lists of contrasts, with the chosen criteria, in bed format. Useful to use with HOMER, GREAT... Be aware that EPIC is annotated
               with hg19 genome."),
             h3("Download Markdown Report:"),
             downloadButton("download_export_markdown"),
             p("Press to download the RMarkdown report of all the steps follow and selected in the pipeline."),
             h3("Download Heatmap:"),
             downloadButton("download_export_heatmaps"),
             p("Press to download the custom heatmap in the gplots::heatmap.2 version.")
             
             ),
             
      tabPanel("About",
            tags$head(tags$style(HTML("a {color: #0275d8}"))),
            h1("ÉPICo!") ,
            h4(tags$a(href="https://www.gnu.org/licenses/gpl-3.0.html", "GNU GPLv3 License")),
            h4("© 2020 Octavio Morante-Palacios"),
            h4(tags$a(href="mailto:octaviompa@gmail.com?ÉPICo!","octaviompa@gmail.com")),
            p("ÉPICo! is a graphical interface based on Shiny. It is intended for importing and normalizing methylation EPIC array data, and also for following a statistical analysis to
              detect differentially methylated CpGs, and plotting them in useful and customable heatmaps. For this purpose, several packages, such as minfi, limma, gplots, heatmaply, shinycssloaders, 
              doParallel and data.table, are used."),
            p("For suggestions or bug reports, please send me an email or use the", tags$a(href="https://github.com/omorante","GitHub"), "forum.")

            ) 
      
    )






server <- function(input, output, session) {
  
  
  onStop(function() {
    print("Cleaning session...") 
    
    print(ls())
    rm(list=ls())
    gc()
  })
  
  #INPUT
  
  #Load button only shows if file is uploaded
  output$ui_button_input_load = renderUI({
    if(!is.null(input$fileinput_input$datapath)) return(actionButton("button_input_load", "Load Data"))
    else return()
  })
  
  #When you press button_input_load, the data is unzipped and the metharray sheet is loaded
  
  rval_sheet_target = eventReactive(input$button_input_load, rval_sheet()[rval_sheet()[,input$select_input_samplenamevar]  %in% input$selected_samples,])
  
  rval_sheet = eventReactive(input$button_input_load, {
    unzip(input$fileinput_input$datapath, exdir= paste0(tempdir(), "/experiment_data"))
    read.metharray.sheet( paste0(tempdir(), "/experiment_data"))
    #read.metharray.sheet(readDirectoryInput(session, 'sample_directory'))
  } )
  
  #When you press button_input_load, the form options are updated
  observeEvent(input$button_input_load, {
    updateSelectInput(session, "select_input_samplenamevar", label="Select Sample Names Column:", choices = colnames(rval_sheet()))
    updateSelectInput(session, "select_input_groupingvar", label="Select Grouping Variable:", choices = colnames(rval_sheet()))
    updateSelectInput(session, "select_input_donorvar", label="Select Donor Variable:", choices = colnames(rval_sheet()))
    })
  
  #the checkbox of samples to process is updated when samplenamevar changes
  observeEvent(input$select_input_samplenamevar,updateCheckboxGroupInput(session, "selected_samples", label="Select Samples to Process:", choices= rval_sheet()[, input$select_input_samplenamevar], selected = rval_sheet()[, input$select_input_samplenamevar] ))
  
  #The dataframe is rendered
  output$samples_table =  DT::renderDT(rval_sheet())
  
  #rval_rgset loads RGSet using read.metharray.exp and the sample sheet (rval_sheet())
  rval_rgset = eventReactive(input$button_input_next,{
    targets = rval_sheet()[rval_sheet()[,input$select_input_samplenamevar]  %in% input$selected_samples,]
    print(targets)
    RGSet = read.metharray.exp(targets = targets,verbose=TRUE)
    RGSet [(rowMeans(as.matrix(detectionP(RGSet))) <0.01),]
  })
  
  #We change the page to the next one
  observeEvent(input$button_input_next, 
    { withProgress(message="Loading data...",value = 2, max=5, 
      {rval_rgset()
      updateNavbarPage(session, "navbar_epic","Minfi Norm.")})
    })
  
  
  #MINFI NORMALIZATION
  
  
  #Calculation of minfi normalized data
  rval_gset = reactive({
    
    if(input$select_minfi_norm == "Illumina"){
      gset = mapToGenome(ratioConvert(preprocessIllumina(rval_rgset(), bg.correct=TRUE, normalize="controls")))
    }
    
    else if(input$select_minfi_norm == "Raw"){
      gset = mapToGenome(ratioConvert(preprocessRaw(rval_rgset())))
    }  
    
    else if(input$select_minfi_norm == "Noob"){
      gset = mapToGenome(ratioConvert(preprocessNoob(rval_rgset()))) 
    }  
    
    else if(input$select_minfi_norm == "Funnorm"){
      gset = preprocessFunnorm((rval_rgset()))
    }
    else if(input$select_minfi_norm == "Quantile"){
      gset = preprocessQuantile(rval_rgset())
    }
    
    else if (input$select_minfi_norm == "Noob+Quantile"){
      gset = preprocessQuantile(preprocessNoob(rval_rgset()))
    }
    
    #remove SNPs to proceed with the analysis
    addSex(dropLociWithSnps(gset, maf = 0))
    
  })
  
  
  #getBeta/getM reactives
  rval_rgset_getBeta = reactive({
    bvalues= as.data.frame(minfi::getBeta(rval_rgset()))
    colnames(bvalues) = input$selected_samples
    bvalues
    })
  
  rval_gset_getBeta = reactive({
    bvalues = as.data.frame(minfi::getBeta(rval_gset()))
    colnames(bvalues) = input$selected_samples
    bvalues
  })
  

  rval_gset_getM = reactive({
       minfi::getM(rval_gset())})
  ###
  

  
  #Plotting of QC graphs
      output$graph_minfi_qcraw = renderCachedPlot(plotQC(getQC(preprocessRaw(rval_rgset()))), paste0("Raw", "QC"))
      output$graph_minfi_densityplotraw = renderCachedPlot(densityPlot(rval_rgset()),paste0("Raw", "densityplot"))
      output$graph_minfi_densitybeanplotraw = renderCachedPlot(densityBeanPlot(rval_rgset()),paste0("Raw", "densitybeanplot"))
      output$graph_minfi_mdsplotraw = renderCachedPlot(mdsPlot(rval_rgset(),sampGroups = input$select_input_groupingvar), paste0("Raw", "mdsplot"))
      output$graph_minfi_boxplotraw = renderCachedPlot(boxplot(as.matrix(rval_rgset_getBeta())), paste0(input$select_minfi_norm, "boxplot"))
      
    
      output$graph_minfi_qc = renderCachedPlot(plotQC(getQC(rval_gset())), paste0(input$select_minfi_norm, "QC"))
      output$graph_minfi_densityplot = renderCachedPlot(densityPlot(as.matrix(rval_gset_getBeta())),paste0(input$select_minfi_norm, "densityplot"))
      output$graph_minfi_densitybeanplot = renderCachedPlot(densityBeanPlot(as.matrix(rval_gset_getBeta())),paste0(input$select_minfi_norm, "densitybeanplot"))
      output$graph_minfi_mdsplot = renderCachedPlot(mdsPlot(as.matrix(rval_gset_getBeta()),sampGroups = input$select_input_groupingvar), paste0(input$select_minfi_norm, "mdsplot"))
      output$graph_minfi_boxplot = renderCachedPlot(boxplot(as.matrix(rval_gset_getBeta())), paste0(input$select_minfi_norm, "boxplot"))
      
      output$graph_minfi_sex = renderCachedPlot(plotSex(rval_gset()), paste0(input$select_minfi_norm, "sex"))
      
      #SNPs heatmap
      
      rval_plot_minfi_snps = reactive({
        
        buylrd = c("#313695", "#4575B4", "#74ADD1", "#ABD9E9", "#E0F3F8", "#FFFFBF",
                    "#FEE090", "#FDAE61", "#F46D43", "#D73027", "#A50026")
        colors.martin = colorRampPalette(buylrd)(100)
        snps = getSnpBeta(rval_rgset())
        rval_sheet_target =  rval_sheet()[rval_sheet()[,input$select_input_samplenamevar]  %in% input$selected_samples,]
        
        colnames(snps) = rval_sheet_target[, input$select_input_samplenamevar]
        heatmaply::heatmaply(snps,
                             col = colors.martin, Colv = T, key.title = "", na.rm = T, dendogram = "both", scale = "row", col_side_colors = rval_sheet_target[, input$select_input_donorvar],
                             distfun = "pearson", hclustfun = function(x) hclust(x, method = "average"),
                             seriate = "mean", row_dend_left = TRUE, showticklabels = c(TRUE, FALSE), branches_lwd = 0.3, plot_method = "plotly", colorbar_xpos = -0.01, colorbar_ypos=0.3, margins=c(25,25,NA,0) )
        
      })
      
      
      output$graph_minfi_snps = renderPlotly(rval_plot_minfi_snps())

      
     

  
  
  #Update of next form and move to Limma
  observeEvent(input$button_minfi_select,
               {
                updateSelectInput(session, "select_limma_voi", label="Select Variable of Interest", choices = colnames(pData(rval_gset())),selected = input$select_input_groupingvar)
                updateCheckboxGroupInput(session, "checkbox_limma_covariables", label="Select Covariables to block", choices = colnames(pData(rval_gset())),selected = input$select_input_donorvar)
                updateCheckboxGroupInput(session, "checkbox_limma_groups", label="Select Groups to compare", choices = unique(pData(rval_gset())[,input$select_input_groupingvar]))
                updateNavbarPage(session, "navbar_epic","Limma")
      
                updateSelectInput(session,"select_limma_trend", label="eBayes trend option:", c("TRUE","FALSE"), "FALSE")
                updateSelectInput(session,"select_limma_robust", label="eBayes robust option:", c("TRUE","FALSE"),"FALSE")
                })
  
  
  

  #LIMMA
  
  #Variable of interest
  rval_voi = reactive(factor(pData(rval_gset())[,input$select_limma_voi]))
  
  #Design calculation
  rval_design = reactive({
    #design = model.matrix(~ 0 + rval_voi())
    pdata = pData(rval_gset())
    formula = as.formula(paste0("~ 0 + ", paste(c(input$select_limma_voi, input$checkbox_limma_covariables), collapse = "+" )))
    print(paste0("~ 0 + ", paste(c(input$select_limma_voi, input$checkbox_limma_covariables), collapse = "+" )))
    design = model.matrix(formula,data=pdata) 
    print(design)
    colnames(design)[1:length(unique(rval_voi()))] = levels(rval_voi())
    design
   })
  
  #Calculation of contrasts
  rval_contrasts = reactive({
    contrastes = c() #calculamos lista de contrastes de todos contra todos
    valor = 1
    all.groups = levels(rval_voi())
    for (i in 1:(length(all.groups)-1)){
      for (z in 1:(length(all.groups) - valor)){
        contrastes = c(contrastes, paste(all.groups[i], all.groups[i+z],sep = "-"))
      }
      valor <- valor + 1}
      contrastes
  })
  
  
  #Calculation of limma model
  rval_fit = eventReactive(input$button_limma_calculatemodel, {
    limma::lmFit(rval_gset_getM(), rval_design())
  })
  
  
  #f rval_fit() has NAs, we remove the option to trend or robust in eBayes to prevent failure
  observeEvent(input$button_limma_calculatemodel,{
    if (any(is.na(rval_fit()$coefficients))){
      updateSelectInput(session,"select_limma_trend", label="Norm. method is not compatible with trend", "FALSE", "FALSE")
      updateSelectInput(session,"select_limma_robust", label="Norm. method is not compatible with robust", "FALSE","FALSE")
      }
  })
  
  
  #Calculation of global difs
  rval_globaldifs = eventReactive(input$button_limma_calculatedifs,{
    calculate_global_difs(rval_gset_getBeta(), rval_voi(), rval_contrasts())}
  )
  
  #Calculation of differences (eBayes)
  rval_finddifcpgs = eventReactive(input$button_limma_calculatedifs,{
    find_dif_cpgs(rval_voi(), rval_design(), rval_fit(), rval_contrasts(), trend = as.logical(input$select_limma_trend), robust = as.logical(input$select_limma_robust))
  })
  
  #Update of heatmap controls
  observeEvent(input$button_limma_calculatedifs, {
    updateTabsetPanel(session,"tabset_limma","differential_cpgs")
    updateSelectInput(session, "select_limma_groups2plot", label="Groups to plot", choices = levels(rval_voi()), selected = levels(rval_voi()))
    updateSelectInput(session, "select_limma_contrasts2plot", label="Contrasts to plot", choices = rval_contrasts(), selected = rval_contrasts())
    updateSliderInput(session, "slider_limma_deltab", "min DeltaBeta", min=0, max=1, value=0.2)
    #force rval_filteredlist
    rval_filteredlist()
  })
  
  #Calculation of filtered list
  rval_filteredlist = reactive({
    withProgress(message="Performing contrasts calculations...", value=1, max=6, {
      setProgress(message ="Calculating global difs...", value=1)
      rval_globaldifs()
      setProgress(message = "Calculating eBayes...", value=4)
      rval_finddifcpgs()
      setProgress(message = "Calculating filtered list...", value=5)
      create_filtered_list(rval_finddifcpgs(), rval_globaldifs(), deltaB = input$slider_limma_deltab, adjp_max = input$slider_limma_adjpvalue, p.value = 0.05)})
    })


  #render of plots and tables
  output$graph_limma_plotSA = renderPlot(limma::plotSA(rval_fit()))
  output$graph_limma_plotMA = renderPlot(limma::plotMA(rval_gset_getM()))
 
  
  plot_heatmap = eventReactive(input$button_limma_heatmapcalc,
                               create_heatmap(rval_filteredlist(), rval_gset_getBeta(), factorgroups =  factor(rval_voi()[rval_voi() %in% input$select_limma_groups2plot], 
                              levels=input$select_limma_groups2plot),groups2plot = rval_voi() %in% input$select_limma_groups2plot, 
                              contrasts2plot = rval_contrasts() %in% input$select_limma_contrasts2plot, 
                              Colv = as.logical(input$select_limma_colv), clusteralg = input$select_limma_clusteralg, 
                              distance = input$select_limma_clusterdist, scale = input$select_limma_scale, static=as.logical(input$select_limma_graphstatic))
                              )
  make_table = eventReactive(input$button_limma_heatmapcalc, 
                             rbindlist(rval_filteredlist()[rval_contrasts() %in% input$select_limma_contrasts2plot], 
                             idcol="contrast") %>% mutate(type=ifelse(dif_current > 0, "Hypermethylated", "Hypomethylated")) %>% group_by(contrast, type)  %>% 
                              summarise(CpGs=dplyr::n()) %>% pivot_wider(names_from=type, values_from=CpGs) %>% mutate(total = Hypermethylated + Hypomethylated)
                             )
  
  observeEvent(input$button_limma_heatmapcalc,{
    #Render the correct plot depending on the selected
    output$graph_limma_heatmapcontainer = renderUI({
      if (!as.logical(input$select_limma_graphstatic)) return( plotlyOutput("graph_limma_heatmap_interactive",width = "600px", height="500px" )  %>% withSpinner() )
      else return( plotOutput("graph_limma_heatmap_static",width = "600px", height="630px" )%>% withSpinner())
    })
    
    if (!as.logical(input$select_limma_graphstatic)) 
      output$graph_limma_heatmap_interactive = renderPlotly(plot_heatmap())
    else 
      output$graph_limma_heatmap_static = renderPlot(plot_heatmap())
  })
  

  
  

  output$table_limma_difcpgs = renderTable(make_table())
  
  
  #EXPORT
  
  rval_annotation = reactive( getAnnotation(rval_gset()) )
  
  #R Objects
  output$download_export_robjects = downloadHandler("R_Objects.zip", content= function(file) {
    
    oldwd = getwd()
    setwd(tempdir())
    
    rgset = rval_rgset()
    gset = rval_gset()
    fit = rval_fit()
    design = rval_design()
    ebayes_tables = rval_finddifcpgs()
    bvalues = rval_gset_getBeta()
    mvalues = rval_gset_getM()
    annotation = rval_annotation()
      
    #save(bvalues,file="Bvalues.RData")
    #save(mvalues,file="Mvalues.RData")
    save(rgset , file="RGSet.RData")
    save(gset, file="Normalized_GenomicRatioSet.RData")
    save(fit, file="fit.RData")
    save(design, file="design.RData")
    #save(annotation,file="annotation.RData")
    save(ebayes_tables, file="ebayestables.RData")
    
    objects = list.files(pattern="*.RData", full.names=T)
    
    setwd(oldwd)
    
    zip(file, objects, flags="-j9X")
  })
  
  
  #Filtered BEDs
  output$download_export_filteredbeds  = downloadHandler("filtered_beds.zip", content = function(file){
    
    filtered_beds = create_filtered_beds(rval_filteredlist(), rval_annotation())
    
    lapply(names(filtered_beds), function(x){fwrite(filtered_beds[[x]], file = paste0(dirname(file),"/", x,".bed"), sep="\t", quote=F, col.names=F, row.names=F)})
    
    objects = list.files(path=dirname(file), full.names=T, pattern="*.bed")
    zip(file, objects, flags="-j9X")
  })
  
  
  #Markdown Report
  
  output$download_export_markdown = downloadHandler("Markdown.html", content = function(file){
    
    params = list(name_var = input$select_input_samplenamevar, grouping_var = input$select_input_groupingvar, donor_var = input$select_input_donorvar, normalization_mode = input$select_minfi_norm, 
                  rval_rgset = rval_rgset(), rval_gset = rval_gset(), rval_sheet_target = rval_sheet_target(), rval_fit = rval_fit(), rval_design = rval_design(), rval_filteredlist = rval_filteredlist(), rval_contrasts=rval_contrasts(),
                  limma_ebayes_trend = input$select_limma_trend, limma_ebayes_robust = input$select_limma_robust, limma_voi = input$select_limma_voi, limma_covar = input$checkbox_limma_covariables, groups2plot = input$select_limma_groups2plot, 
                  contrasts2plot = input$select_limma_contrasts2plot, Colv = input$select_limma_colv, clusteralg = input$select_limma_clusteralg, distance = input$select_limma_clusterdist, scale = input$select_limma_scale, max_fdr=input$slider_limma_adjpvalue, min_deltabeta = input$slider_limma_deltab  )
   
    rmarkdown::render("report.Rmd",output_file=file, params=params, envir= new.env(parent = globalenv()))
    
  })
  
  
  #custom heatmap
  
  output$download_export_heatmaps = downloadHandler("Custom_heatmap.pdf", content=function(file){
    
    pdf(file = file, height = 11.71, width=8.6)
    create_heatmap(rval_filteredlist(), rval_gset_getBeta(), factorgroups =  factor(rval_voi()[rval_voi() %in% input$select_limma_groups2plot], levels=input$select_limma_groups2plot),groups2plot = rval_voi() %in% input$select_limma_groups2plot, contrasts2plot = rval_contrasts() %in% input$select_limma_contrasts2plot, Colv = as.logical(input$select_limma_colv), clusteralg = input$select_limma_clusteralg, distance = input$select_limma_clusterdist, scale = input$select_limma_scale, static=TRUE)
    dev.off()
    
  }
    
  )
  
  #Graphs
  
  

  
}


shinyApp(ui, server)

