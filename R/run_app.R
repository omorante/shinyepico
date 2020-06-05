#' Run Shiny ÉPICo!
#'
#' ShinyÉPICo! Interactive minfi and limma pipeline for Illumina methylation arrays
#'
#' @param n_cores Number of cores to be used in parallelized operations in the application. By default, half of your CPU cores. As more cores mean more RAM used, we don't recommend use more than 1 core per 4GB of RAM available
#' @param max_upload_size The limit in MB of the .zip file size to be uploaded. By default, 2000MB. 
#' @return None
#' @examples
#' \donttest{run_shinyepico()}
#' 
#' @importFrom shiny shinyApp
#' @importFrom golem with_golem_options
#' @export
#' 
run_shinyepico <- function(
  n_cores = 1,
  max_upload_size = 2000,
  host = "0.0.0.0",
  port = 80

) {
  with_golem_options(
    app = shinyApp(
      ui = app_ui, 
      server = app_server, options = list(host = host, port=port)
    ), 
    golem_opts = list(n_cores = n_cores, max_upload_size = max_upload_size * 1024^2)
  )
}

