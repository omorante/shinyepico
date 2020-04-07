

# ENVIRONMENT VARIABLES
norm_options = c("Raw",
                 "Illumina",
                 "Funnorm",
                 "Noob",
                 "Quantile",
                 "Noob+Quantile")

hclust_methods = c("single",
                   "complete",
                   "average",
                   "mcquitty",
                   "median",
                   "centroid")


`%dopar%` = foreach::`%dopar%`
`%>%` = magrittr::`%>%`

#' The application User-Interface
#'
#' @param request Internal parameter for `{shiny}`.
#'     DO NOT REMOVE.
#' @import shiny
#' @noRd
app_ui <- function(request) {
  navbarPage(
    "\u00C9PICo!",
    id = "navbar_epic",
    theme = shinythemes::shinytheme("sandstone") ,
    tabPanel(
      "Input",
      titlePanel("Input sheet file"),
      sidebarLayout(
        sidebarPanel(
          fileInput(
            inputId = "fileinput_input",
            "Upload Compress Experiment Directory (zip):",
            multiple = FALSE,
            accept = "application/zip"
          ),
          #directoryInput("sample_directory", "Select Sample Directory:", value="/home/octavio/Documentos/Análisis DCDEX/Pruebas limma metilación/Arrays de metilación/"),
          uiOutput("ui_button_input_load"),
          h3(),
          
          conditionalPanel(
            "typeof output.samples_table != 'undefined'",
            selectInput("select_input_samplenamevar", "", choices =
                          c()),
            h3(),
            checkboxGroupInput("selected_samples", "", c(), selected =
                                 TRUE),
            h3(),
            selectInput("select_input_groupingvar", "", c()),
            h3(),
            selectInput("select_input_donorvar", "", c()),
            h3(),
            actionButton("button_input_next", "Proceed to the next step")
          )
        ),
        mainPanel(
          DT::DTOutput("samples_table") %>% shinycssloaders::withSpinner()
        )
      )
      
    ),
    
    
    tabPanel(
      "Minfi Norm.",
      titlePanel("Minfi Normalization"),
      sidebarLayout(
        sidebarPanel(
          selectInput("select_minfi_norm", "Select Normalization", norm_options),
          #actionButton("button_minfi_check", "Check"),
          actionButton("button_minfi_select", "Select")
        ),
        mainPanel(
          tabsetPanel(
            tabPanel(
              "Density plot",
              h4("Raw"),
              plotOutput("graph_minfi_densityplotraw") %>% shinycssloaders::withSpinner(),
              h4("Processed"),
              plotOutput("graph_minfi_densityplot") %>% shinycssloaders::withSpinner()
            ),
            
            tabPanel(
              "Density Bean Plot",
              h4("Raw"),
              plotOutput("graph_minfi_densitybeanplotraw") %>% shinycssloaders::withSpinner(),
              h4("Processed"),
              plotOutput("graph_minfi_densitybeanplot") %>% shinycssloaders::withSpinner()
            ),
            
            tabPanel(
              "MDS plot",
              h4("Raw"),
              plotOutput("graph_minfi_mdsplotraw") %>% shinycssloaders::withSpinner(),
              h4("Processed"),
              plotOutput("graph_minfi_mdsplot") %>% shinycssloaders::withSpinner()
            ),
            
            tabPanel(
              "Boxplot",
              h4("Raw"),
              plotOutput("graph_minfi_boxplotraw") %>% shinycssloaders::withSpinner(),
              h4("Processed"),
              plotOutput("graph_minfi_boxplot") %>% shinycssloaders::withSpinner()
            ),
            
            tabPanel(
              "SNPs Heatmap",
              h4("SNPs beta-values (Raw)"),
              plotly::plotlyOutput("graph_minfi_snps") %>% shinycssloaders::withSpinner()
            ),
            
            tabPanel(
              "Sex prediction",
              h4("plotSex"),
              plotOutput("graph_minfi_sex") %>% shinycssloaders::withSpinner()
            ),
            
            tabPanel(
              "QC plot",
              h4("Raw"),
              plotOutput("graph_minfi_qcraw") %>% shinycssloaders::withSpinner()
            )
            
            
          )
        )
      )
    ),
    
    tabPanel(
      "Limma",
      titlePanel("Limma parameters"),
      sidebarLayout(
        sidebarPanel(
          width = 3,
          selectInput("select_limma_voi", "Select Variable of Interest", c()),
          checkboxGroupInput(
            "checkbox_limma_covariables",
            "Select covariables to block",
            c()
          ),
          #checkboxGroupInput("checkbox_limma_groups", "Select groups to compare", c()),
          selectInput(
            "select_limma_trend",
            "eBayes trend option:",
            c("TRUE", "FALSE"),
            c("FALSE")
          ),
          selectInput(
            "select_limma_robust",
            "eBayes robust option:",
            c("TRUE", "FALSE"),
            c("FALSE")
          ),
          actionButton("button_limma_calculatemodel", "Generate Model"),
          actionButton("button_limma_calculatedifs", "Calculate Contrasts")
          
        ),
        mainPanel(
          width = 9,
          tabsetPanel(
            id = "tabset_limma",
            tabPanel(
              "Model diagnosis",
              value = "model_diagnosis",
              plotOutput("graph_limma_plotSA") %>% shinycssloaders::withSpinner(),
              plotOutput("graph_limma_plotMA")
            ),
            tabPanel(
              "Differential CpGs",
              value = "differential_cpgs",
              div(
                style = 'max-width:800px;margin:auto;',
                fluidPage(
                  h4("Heatmap"),
                  uiOutput("graph_limma_heatmapcontainer"),
                  h4("CpGs selected with these filters:"),
                  tableOutput("table_limma_difcpgs") %>% shinycssloaders::withSpinner(),
                  column(
                    6,
                    h4("Group options"),
                    selectizeInput("select_limma_groups2plot", "Groups to plot", c(), multiple = TRUE),
                    selectizeInput(
                      "select_limma_contrasts2plot",
                      "Contrasts to plot",
                      c(),
                      multiple = TRUE
                    )
                  ),
                  column(
                    6,
                    h4("Filtering options"),
                    sliderInput("slider_limma_deltab", "Min. DeltaBeta", 0, 1, 0.2),
                    sliderInput("slider_limma_adjpvalue", "Max. FDR", 0, 1, 0.05),
                    sliderInput("slider_limma_pvalue", "Max. p-value", 0, 1, 1)

                  ),
                  h4("Clustering options", align =
                       "center"),
                  column(
                    6,
                    selectInput(
                      "select_limma_clusteralg",
                      "Clustering algorithm:",
                      c(
                        "single",
                        "complete",
                        "average",
                        "mcquitty",
                        "median",
                        "centroid"
                      ),
                      "average"
                    ),
                    selectInput(
                      "select_limma_clusterdist",
                      "Distance Function:",
                      c("pearson", "spearman", "kendall", "euclidean"),
                      "pearson"
                    ),
                    selectInput("select_limma_graphstatic", label = "Plot Static Graph:", c(TRUE, FALSE), c(TRUE))
                    
                  ),
                  
                  column(
                    6,
                    selectInput(
                      "select_limma_colv",
                      "Create Column Dendogram:",
                      c(TRUE, FALSE),
                      c(TRUE)
                    ),
                    selectInput("select_limma_scale", "Scale:", c("row", "none"), "row"),
                    tags$br(),
                    actionButton("button_limma_heatmapcalc", "Update")
                  ),
                  tags$br(),
                  tags$br()
                  
                  
                )
              )
            )
          )
        )
      )
    ),
    
    tabPanel(
      "Export",
      
      h3("Download RObjects:"),
      downloadButton("download_export_robjects"),
      p(
        "Press to download the R objects used for the analysis (RGSet, GenomicRatioSet, Bvalues, Mvalues, etc."
      ),
      h3("Download filtered bed files:"),
      downloadButton("download_export_filteredbeds"),
      p(
        "Press to download the created filtered lists of contrasts, with the chosen criteria, in bed format. Useful to use with HOMER, GREAT... Be aware that EPIC is annotated
        with hg19 genome."
      ),
      h3("Download Markdown Report:"),
      downloadButton("download_export_markdown"),
      p(
        "Press to download the RMarkdown report of all the steps follow and selected in the pipeline."
      ),
      h3("Download Heatmap:"),
      downloadButton("download_export_heatmaps"),
      p(
        "Press to download the custom heatmap in the gplots::heatmap.2 version."
      )
      
    ),
    
    tabPanel(
      "About",
      tags$head(tags$style(HTML(
        "a {color: #0275d8}"
      ))),
      h1("\u00C9PICo!") ,
      h4(
        tags$a(href = "https://www.gnu.org/licenses/gpl-3.0.html", "GNU GPLv3 License")
      ),
      h4("\u00A9 2020 Octavio Morante-Palacios"),
      h4(
        tags$a(href = "mailto:octaviompa@gmail.com?\u00C9PICo!", "octaviompa@gmail.com")
      ),
      p(
        "\u00C9PICo! is a graphical interface based on Shiny. It is intended for importing and normalizing methylation EPIC array data, and also for following a statistical analysis to
        detect differentially methylated CpGs, and plotting them in useful and customable heatmaps. For this purpose, several packages, such as minfi, limma, gplots, heatmaply, shinycssloaders,
        doParallel and data.table, are used."
      ),
      p(
        "For suggestions or bug reports, please send me an email or use the",
        tags$a(href = "https://github.com/omorante", "GitHub"),
        "forum."
      )
      
    )
    
  )
  
}

#' Add external Resources to the Application
#'
#' This function is internally used to add external
#' resources inside the Shiny application.
#'
#' @import shiny
#' @importFrom golem add_resource_path activate_js favicon bundle_resources
#' @noRd
#' 
golem_add_external_resources <- function() {
  add_resource_path('www', app_sys('app/www'))
  
  tags$head(favicon(),
            bundle_resources(path = app_sys('app/www'),
                             app_title = 'shiny.epico.golem'))
  # Add here other external resources
  # for example, you can add shinyalert::useShinyalert() )
}

