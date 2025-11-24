#' Launch the Shiny App
#'
#' This function launches the Shiny application.
#'
#' @return No return value, called for side effects (launches Shiny app)
#'
#' @importFrom shiny runApp
#' @importFrom shinydashboard dashboardPage dashboardHeader dashboardSidebar dashboardBody sidebarMenu menuItem
#' @importFrom ggplot2 ggplot aes geom_line geom_point geom_ribbon labs theme_minimal theme element_text
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr select filter mutate
#' @importFrom DT dataTableOutput renderDataTable datatable
#' @importFrom MASS mvrnorm
#' @importFrom LaplacesDemon rmvt rinvwishart
#' @importFrom shinyjs useShinyjs reset
#' @importFrom stats na.omit sd median quantile rchisq rnorm rt t.test
#' @importFrom graphics grid polygon
#' @importFrom grDevices rgb
#' @export
powerprior_launch_shinyapp <- function() {
  app_dir <- system.file("shiny", "app", package = "powerprior")
  if (app_dir == "") {
    stop("Could not find the Shiny application directory.
              Ensure 'myapp_folder_name' is correctly placed in inst/shiny.")
  }
  shiny::runApp(app_dir, display.mode = "normal")
}
