


#' Run Shiny ÉPICo!
#'
#' ShinyÉPICo! Interactive minfi and limma pipeline for Illumina methylation arrays
#'
#' @param n_cores Number of cores to be used in parallelized operations in the application. By default, half of your CPU cores. Parallelization affects only to mean and differences calculation and it does not suppose a significant memory overhead.
#' @param max_upload_size The limit in MB of the .zip file size to be uploaded. By default, 2000MB.
#' @param host IP used to deploy the server. By default, your local IP (127.0.0.1)
#' @param port Port used to deploy the server.
#' @param seed Seed used to make DMR calculations reproducible. Set it to NULL to not use a seed. DMR calculation is based on permutations, and is subject to uncertainty. Therefore, you can get slightly different results every time you use the application if you do not set a seed.
#' @return None
#' @examples
#' {
#'   if (interactive()) {
#'     run_shinyepico()
#'   }
#' }
#' @importFrom shiny shinyApp
#' @importFrom shiny shinyOptions
#' @export
#'
run_shinyepico <- function(n_cores = parallel::detectCores() / 2,
                           max_upload_size = 2000, host = "127.0.0.1", port = NULL, seed = 123) {

  stopifnot(is.numeric(n_cores) & n_cores >= 1)
  stopifnot(is.numeric(max_upload_size) & max_upload_size > 0)
  stopifnot(is.numeric(seed) | is.null(seed))
  stopifnot(is.numeric(port) | is.null(port))


  # Specified seed to obtain always the same results
  set.seed(seed)

  # Specified n_cores and max_upload_size
  shinyOptions(n_cores = round(n_cores), shiny.maxRequestSize = max_upload_size *
    1024^2)

  shinyApp(ui = app_ui, server = app_server, options = list(
    host = host,
    port = port, launch.browser = TRUE
  ))
}
