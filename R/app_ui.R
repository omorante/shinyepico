

# ENVIRONMENT VARIABLES
norm_options = c("Raw",
                 "Illumina",
                 "Funnorm",
                 "Noob",
                 "SWAN",
                 "Quantile",
                 "Noob+Quantile")

hclust_methods = c("single",
                   "complete",
                   "average",
                   "mcquitty",
                   "median",
                   "centroid")


`%dopar%` = foreach::`%dopar%`
`%do%` = foreach::`%do%`
`%>%` = magrittr::`%>%`

#' The application User-Interface
#'
#' @param request Internal parameter for `{shiny}`.
#'     DO NOT REMOVE.
#' @import shiny
#' @import shinyWidgets
#' @noRd
app_ui <- function(request) {
  
  navbarPage(
    "\u00C9PICo!",
    id = "navbar_epic",
    theme = shinythemes::shinytheme("sandstone") ,
    tabPanel(
      "Input",
      titlePanel("Load iDATs and sample sheet"),
      sidebarLayout(
        sidebarPanel(width = 3,
          fileInput(
            inputId = "fileinput_input",
            "Upload Compress Experiment Directory (zip):",
            multiple = FALSE,
            accept = c("application/zip", "application/octet-stream", "application/x-zip-compressed")
          ),
          
          uiOutput("ui_button_input_load"),
          h3(),
          
          conditionalPanel(
            "typeof output.samples_table != 'undefined'",
            selectInput("select_input_samplenamevar", "", choices =
                          c()),

            h3(),
            selectInput("select_input_groupingvar", "", c()),
            h3(),
            selectInput("select_input_donorvar", "", c()),
            h3(),
            pickerInput(
              inputId = "selected_samples",
              label = "",
              choices = c(),
              options = list(
                `actions-box` = TRUE,
                size = 10,
                `selected-text-format` = "count > 3"
              ),
              multiple = TRUE
            ),
            h3(),
            actionButton("button_input_next", "Continue")
          )
        ),
        mainPanel(width = 9,
          DT::DTOutput("samples_table") %>% shinycssloaders::withSpinner()
        )
      )
      
    ),
    
    
    tabPanel(
      "Normalization",
      titlePanel("Normalization"),
      sidebarLayout(
        sidebarPanel(width = 3,
          selectInput("select_minfi_norm", "Select Normalization", norm_options),
          shinyjs::disabled(actionButton("button_minfi_select", "Select"))
        ),
        mainPanel(width=9,
          tabsetPanel(
            
            tabPanel(
              "Quality Control",
              h4("Overall Signal"),
              plotly::plotlyOutput("graph_minfi_qcraw") %>% shinycssloaders::withSpinner(),
              h4("Bisulfite Conversion II"),
              plotly::plotlyOutput("graph_minfi_bisulfiterawII")  %>% shinycssloaders::withSpinner()
            ),
            
            tabPanel(
              "Density plot",
              h4("Raw"),
              plotly::plotlyOutput("graph_minfi_densityplotraw") %>% shinycssloaders::withSpinner(),
              h4("Processed"),
              plotly::plotlyOutput("graph_minfi_densityplot") %>% shinycssloaders::withSpinner()
            ),
            
            tabPanel(
              "Boxplot",
              h4("Raw"),
              plotOutput("graph_minfi_boxplotraw") %>% shinycssloaders::withSpinner(),
              h4("Processed"),
              plotOutput("graph_minfi_boxplot") %>% shinycssloaders::withSpinner()
            ),
            
            tabPanel(
              "Principal Component Analysis",
              h4("Processed"),
              plotly::plotlyOutput("graph_minfi_pcaplot") %>% shinycssloaders::withSpinner(),
              tableOutput("table_minfi_pcaplot") %>% shinycssloaders::withSpinner(),
              column(6,
                     selectInput(inputId = "select_minfi_pcaplot_pcx", 
                                 choices = c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10"), 
                                 selected = "PC1", label = "Select x variable"),
                     
                     selectInput(inputId = "select_minfi_pcaplot_color", choices = c(), 
                                 label = "Select color variable:")
                     ),
              column(6,
                     selectInput(inputId = "select_minfi_pcaplot_pcy", 
                                 choices = c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10"), 
                                 selected = "PC2", 
                                 label = "Select y variable")),
                     actionButton("button_pca_update", "Update")
                     ),
                    
                     
            
            tabPanel(
              "Correlations",
              h4("Processed"),
              plotly::plotlyOutput("graph_minfi_corrplot") %>% shinycssloaders::withSpinner(),
            ),
            
            
            tabPanel(
              "SNPs Heatmap",
              h4("SNPs beta-values (Raw)"),
              plotly::plotlyOutput("graph_minfi_snps") %>% shinycssloaders::withSpinner()
            ),
            
            tabPanel(
              "Sex prediction",
              h4("plotSex"),
              plotly::plotlyOutput("graph_minfi_sex") %>% shinycssloaders::withSpinner()
            )
          )
        )
      )
    ),
    
    tabPanel(
      "DMPs",
      titlePanel("Differentially Methylated Positions"),
      sidebarLayout(
        sidebarPanel(width = 3,
          selectInput("select_limma_voi", "Select Variable of Interest", c()),
          
          pickerInput(
            inputId = "checkbox_limma_covariables",
            label = "Select linear model covariables",
            choices = c(),
            multiple = TRUE,
            options = list(
              `actions-box` = TRUE,
              size = 10,
              `selected-text-format` = "count > 3")
            ),

          #checkboxGroupInput("checkbox_limma_groups", "Select groups to compare", c()),
          

          switchInput(
              inputId = "select_limma_weights",
            label = "Array Weights", 
            labelWidth = "80px",
            value = FALSE,
          ),
          
          
          shinyjs::disabled(actionButton("button_limma_calculatemodel", "Generate Model")),
          tags$br(),
          uiOutput("button_limma_calculatedifs_container")
        ),
        mainPanel(width = 9,
          tabsetPanel(
            id = "tabset_limma",
            tabPanel(
              "Model diagnosis",
              value = "model_diagnosis",
              h4("Sigma vs A plot"),
              plotOutput("graph_limma_plotSA") %>% shinycssloaders::withSpinner(),
              #h4("Log-intensity ratio vs Average Plot"),
              #plotOutput("graph_limma_plotMA") %>% shinycssloaders::withSpinner()
            ),
            tabPanel(
              "Differential CpGs",
              value = "differential_cpgs",
              div(
                style = 'max-width:800px;margin:auto;',
                fluidPage(
                  h4("Heatmap"),
                  uiOutput("graph_limma_heatmapcontainer"),
                  h4("DMP counts in each contrast:"),
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
                    
                    selectInput("select_limma_scale", "Scale:", c("row", "none"), "row"),
                    tags$br(),
                  ),
                  
                  column(
                    6,
                    
                    tags$br(),
                    
                    switchInput(
                      inputId = "select_limma_colv",
                      label = "Column Dendogram", 
                      labelWidth = "80px",
                      value = TRUE,
                    ),
                    
                    tags$br(),
                    
                    switchInput(
                      inputId = "select_limma_graphstatic",
                      label = "Static Graph", 
                      labelWidth = "80px",
                      value = TRUE,
                    ),
                    
                  tags$br(),
            
                  shinyjs::disabled(actionButton("button_limma_heatmapcalc", "Update"))
                  ),
                  
                  
                  
                )
              )
            )
          )
        )
      )
    ),
    
    tabPanel(
      "Export",
      shinyjs::useShinyjs(),
      h3("Download RObjects:"),
      downloadButton("download_export_robjects"),
      p(
        "Press to download the R objects used for the analysis (RGSet, GenomicRatioSet, Bvalues, Mvalues, etc."
      ),
      h3("Download filtered bed files:"),
      downloadButton("download_export_filteredbeds"),
      p(
        "Press to download the created filtered lists of contrasts, 
        with the chosen criteria, in bed format. Useful to use with HOMER, GREAT... 
        Be aware that EPIC is annotated with hg19 genome."
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
        detect differentially methylated CpGs, and plotting them in useful and customable heatmaps."
      ),
      p(
        "For suggestions or bug reports, please send me an email or use the",
        tags$a(href = "https://github.com/omorante", "GitHub"), "forum."
      )
      
    )
    
  )
}

