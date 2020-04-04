#' Run Shiny ÉPICo!
#'
#' ShinyÉPICo! Interactive minfi and limma pipeline for Illumina methylation arrays
#'
#' @param n_cores Number of cores to be used in parallelized operations in the application. By default, half of your CPU cores. As more cores mean more RAM used, we don't recommend use more than 1 core per 4GB of RAM available
<<<<<<< HEAD
#' @param max_upload_size The limit in MB of the .zip file size to be uploaded. By default, 2000MB. 
=======
#' @param max_upload_size The limit in MB of the .zip file size to be uploaded. By default, 500MB. 
>>>>>>> 707ab8e45d439906f6ffada0effac4de890932c8
#' @return None
#' @examples
#' \donttest{run_shinyepico()}
#' 
#' @importFrom shiny shinyApp
#' @importFrom golem with_golem_options
#' @export
#' 
run_shinyepico <- function(
<<<<<<< HEAD
  n_cores = if (parallel::detectCores() > 1) {round(parallel::detectCores()/3, digits=0) } else { 1},
  max_upload_size = 2000
=======
  n_cores = if (parallel::detectCores() > 1) {cores = round(parallel::detectCores()/2, digits=0) } else {cores = 1},
  max_upload_size = 500
>>>>>>> 707ab8e45d439906f6ffada0effac4de890932c8
) {
  with_golem_options(
    app = shinyApp(
      ui = app_ui, 
      server = app_server
    ), 
    golem_opts = list(n_cores = n_cores, max_upload_size = max_upload_size * 1024^2)
  )
}

