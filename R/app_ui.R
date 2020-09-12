


# AVAILABLE METHODS
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
    "shiny\u00C9PICo!",
    id = "navbar_epic",
    theme = shinythemes::shinytheme("sandstone") ,
    tabPanel(
      "Input",
      titlePanel("Load iDATs and sample sheet"),
      sidebarLayout(
        sidebarPanel(
          width = 3,
          fileInput(
            inputId = "fileinput_input",
            "Upload Compress Experiment Directory (zip):",
            multiple = FALSE,
            accept = c(
              "application/zip",
              "application/x-zip-compressed"
            )
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
        mainPanel(
          width = 9,
          DT::DTOutput("samples_table") %>% shinycssloaders::withSpinner()
        )
      )
      
    ),
    
    
    tabPanel(
      "Normalization",
      titlePanel("Normalization"),
      sidebarLayout(
        sidebarPanel(
          width = 3,
          selectInput("select_minfi_norm", "Select Normalization", norm_options),
          shinyjs::disabled(actionButton("button_minfi_select", "Select"))
        ),
        mainPanel(
          width = 9,
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
              "SNPs Heatmap",
              h4("SNPs beta-values (Raw)"),
              plotly::plotlyOutput("graph_minfi_snps") %>% shinycssloaders::withSpinner()
             ),
            
            tabPanel(
              "Sex prediction",
              h4("X vs Y chromosomes signal plot"),
              plotly::plotlyOutput("graph_minfi_sex") %>% shinycssloaders::withSpinner(),
              DT::DTOutput("table_minfi_sex") %>% shinycssloaders::withSpinner()
            ),
            
            tabPanel(
              "Principal Component Analysis",
              h4("Processed"),
              plotly::plotlyOutput("graph_minfi_pcaplot") %>% shinycssloaders::withSpinner(),
              DT::DTOutput("table_minfi_pcaplot") %>% shinycssloaders::withSpinner(),
              column(
                6,
                selectInput(
                  inputId = "select_minfi_pcaplot_pcx",
                  choices = c(),
                  label = "Select x variable"
                ),
                
                selectInput(
                  inputId = "select_minfi_pcaplot_color",
                  choices = c(),
                  label = "Select color variable:"
                )
              ),
              column(
                6,
                selectInput(
                  inputId = "select_minfi_pcaplot_pcy",
                  choices = c(),
                  label = "Select y variable"
                )
              ),
              actionButton("button_pca_update", "Update")
            ),
            
            tabPanel(
              "Correlations",
              h4("Processed"),
              plotly::plotlyOutput("graph_minfi_corrplot") %>% shinycssloaders::withSpinner(),
              DT::DTOutput("table_minfi_corrplot") %>% shinycssloaders::withSpinner()
            )
          )
        )
      )
    ),
    
    tabPanel(
      "DMPs",
      titlePanel("Differentially Methylated Positions"),
      sidebarLayout(
        sidebarPanel(
          width = 3,
          h4("Linear Model Options"),
          
          selectInput("select_limma_voi", "Select Variable of Interest", c()),
          
          pickerInput(
            inputId = "checkbox_limma_covariables",
            label = "Select linear model covariables",
            choices = c(),
            multiple = TRUE,
            options = list(
              `actions-box` = TRUE,
              size = 10,
              `selected-text-format` = "count > 3"
            )
          ),
          
          pickerInput(
            inputId = "checkbox_limma_interactions",
            label = "Select linear model interactions",
            choices = c(),
            multiple = TRUE,
            options = list(
              `actions-box` = TRUE,
              size = 10,
              `selected-text-format` = "count > 3"
            )
          ),
          
          switchInput(
            inputId = "select_limma_weights",
            label = "Array Weights",
            labelWidth = "80px",
            value = FALSE,
          ),
          
          
          shinyjs::disabled(
            actionButton("button_limma_calculatemodel", "Generate Model")
          ),
          tags$br(),
          uiOutput("button_limma_calculatedifs_container")
        ),
        mainPanel(
          width = 9,
          tabsetPanel(
            id = "tabset_limma",
            tabPanel(
              "Model diagnosis",
              value = "model_diagnosis",
              h4("Sigma vs A plot"),
              plotOutput("graph_limma_plotSA") %>% shinycssloaders::withSpinner(),
              h4("Design matrix"),
              DT::DTOutput("table_limma_design") %>% shinycssloaders::withSpinner()
            ),
            tabPanel(
              "Differential CpGs",
              value = "differential_cpgs",
              div(
                style = 'max-width:800px;margin:auto;',
                fluidPage(
                  h4("Heatmap"),
                  textOutput("text_limma_heatmapcount"),
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
                    ),
                  ),
                  
                  column(
                    6,
                    h4("Filtering options"),
                    sliderInput("slider_limma_deltab", "Min. DeltaBeta", 0, 1, 0.2),
                    sliderInput("slider_limma_adjpvalue", "Max. FDR", 0, 1, 0.05),
                    sliderInput("slider_limma_pvalue", "Max. p-value", 0, 1, 1)
                    
                  ),
                  
                  h4("Clustering options", align =
                       "left"),
                  
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
                    3,
                    
                    tags$br(),
                    
                    switchInput(
                      inputId = "select_limma_graphstatic",
                      label = "Static Graph",
                      labelWidth = "100px",
                      value = TRUE,
                    ),
                    
                    switchInput(
                      inputId = "select_limma_colv",
                      label = "Column Dendro.",
                      labelWidth = "100px",
                      value = TRUE,
                    ),
                    
                    switchInput(
                      inputId = "select_limma_colsidecolors",
                      label = "Column Colors",
                      labelWidth = "100px",
                      value = FALSE,
                    )
                    
                  ),
                  
                  column(
                    3,
                    
                    tags$br(),
                    
                    switchInput(
                      inputId = "select_limma_rowsidecolors",
                      label = "Row Colors",
                      labelWidth = "100px",
                      value = FALSE,
                    ),
                    
                    conditionalPanel(
                      "input.select_limma_rowsidecolors",
                      numericInput(
                        "select_limma_knumber",
                        "Clusters number",
                        value = 2,
                        min = 1,
                        max = Inf,
                        step = 1
                      )
                    ),
                    
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
        "Press to download the R objects used for the analysis (RGSet, GenomicRatioSet, Bvalues, Mvalues, etc.)"
      ),
      h3("Download filtered bed files:"),
      selectInput("select_export_bedtype",
                  "Subsetting mode",
                  c("by contrasts", "by heatmap cluster"),
                  selected = "by contrast"),
      downloadButton("download_export_filteredbeds"),
      p(
        "Press to download the created filtered lists of contrasts, or heatmap clusters,
        with the chosen criteria, in bed format (hg19 genome)."
      ),
      h3("Download Workflow Report:"),
      downloadButton("download_export_markdown"),
      p(
        "Press to download the report of all the steps follow and selected in the pipeline, and the results."
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
      img(src = "images/logo.png",
          width = 150), 
      h1("shiny\u00C9PICo") ,
      h3("1.0.0-dev"),
      br(),
      h4(
        tags$a(href = "https://www.gnu.org/licenses/agpl-3.0.html", "GNU Affero GPLv3 License", target="_blank")
      ),
      h4("\u00A9 2020 Octavio Morante-Palacios"),
      h4(
        tags$a(href = "mailto:omorante@carrerasresearch.org?\u00C9PICo!", "omorante@carrerasresearch.org", target="_blank")
      ),
      p(
        "shiny\u00C9PICo is a graphical interface based on Shiny. It is intended for importing and normalizing Illumina DNA methylation arrays (450k or EPIC), exploring DNA methylation data and also for following a statistical analysis to
        detect differentially methylated CpGs, and plotting them in useful and customable heatmaps."
      ),
      p(
        "For suggestions or bug reports, please use the",
        tags$a(href = "https://github.com/omorante/shinyepico/issues", "GitHub", target="_blank"),
        "issues forum."
      )
    )
  )
}
