
# UI ----
ui <- dashboardPage(
  dashboardHeader(title = "Power Prior Bayesian Analysis"),
  
  dashboardSidebar(
    sidebarMenu(
      menuItem("Univariate Analysis", tabName = "univariate", icon = icon("line-chart")),
      menuItem("Multivariate Analysis", tabName = "multivariate", icon = icon("th")),
      menuItem("Help", tabName = "help", icon = icon("question-circle"))
    )
  ),
  
  dashboardBody(
    shinyjs::useShinyjs(),
    tags$head(
      tags$meta(charset = "UTF-8"),
      tags$style(HTML("
        .run-button {
          background-color: #28a745 !important;
          border-color: #28a745 !important;
          color: white !important;
          font-size: 18px !important;
          font-weight: bold !important;
          padding: 15px 40px !important;
          border-radius: 8px !important;
          box-shadow: 0 4px 6px rgba(0,0,0,0.2) !important;
          transition: all 0.3s ease !important;
          cursor: pointer !important;
        }
        .run-button:hover {
          background-color: #218838 !important;
          box-shadow: 0 6px 10px rgba(0,0,0,0.3) !important;
          transform: translateY(-2px) !important;
        }
        .run-button:active {
          box-shadow: 0 2px 4px rgba(0,0,0,0.2) !important;
          transform: translateY(0px) !important;
        }
        .parameter-section {
          background-color: #f8f9fa;
          padding: 15px;
          border-radius: 5px;
          margin-bottom: 15px;
        }
        .data-upload-section {
          background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
          color: white;
          padding: 20px;
          border-radius: 8px;
          margin-bottom: 20px;
        }
        .data-upload-section h4 {
          color: white;
          margin-top: 10px;
          font-weight: bold;
        }
        .tab-content-wrapper {
          padding: 20px;
        }
        .results-section {
          margin-top: 30px;
        }
      "))
    ),
    
    tabItems(
      # Univariate Tab ----
      tabItem(tabName = "univariate",
        div(class = "tab-content-wrapper",
          fluidRow(
            column(12,
              h2("Univariate Power Prior Analysis", style = "margin-bottom: 30px;")
            )
          ),
          
          fluidRow(
            column(6,
              div(class = "data-upload-section",
                h4("UPLOAD AND SELECT DATA"),
                h5("Historical Data", style = "color: white;"),
                fileInput("uni_hist_file", "Choose CSV File", 
                         accept = ".csv", width = "100%"),
                uiOutput("uni_hist_col_ui"),
                h5("Current Data", style = "color: white; margin-top: 20px;"),
                fileInput("uni_curr_file", "Choose CSV File", 
                         accept = ".csv", width = "100%"),
                uiOutput("uni_curr_col_ui"),
                fluidRow(
                  column(6, actionButton("uni_load_data", "LOAD DATA", 
                           class = "btn btn-info", style = "width: 100%; margin-top: 10px; color: white; font-weight: bold;")),
                  column(6, actionButton("uni_clear_data", "CLEAR DATA", 
                           class = "btn btn-danger", style = "width: 100%; margin-top: 10px; color: white; font-weight: bold;"))
                )
              )
            ),
            column(6,
              box(title = "CONFIGURATION", width = 12, status = "primary", solidHeader = TRUE,
                div(class = "parameter-section",
                  h5("Discounting Parameter"),
                  sliderInput("uni_a0", "a0 (0 = no borrowing, 1 = full borrowing)", 
                             min = 0, max = 1, value = 0.5, step = 0.05),
                  helpText("Controls how much historical information to use")
                ),
                
                div(class = "parameter-section",
                  h5("Prior Type"),
                  radioButtons("uni_prior_type", NULL,
                              c("Vague (Non-informative)" = "vague",
                                "Informative" = "informative"),
                              selected = "vague", inline = FALSE),
                  
                  conditionalPanel(
                    condition = "input.uni_prior_type == 'informative'",
                    numericInput("uni_mu0", "Prior Mean (mu0)", value = 0),
                    numericInput("uni_kappa0", "Prior Precision (kappa0)", value = 1, min = 0.01),
                    numericInput("uni_nu0", "Prior DoF (nu0)", value = 3, min = 1),
                    numericInput("uni_sigma2_0", "Prior Variance (sigma2_0)", value = 1, min = 0.01),
                    helpText("Weakly informative: small kappa0 and nu0 (1-3)")
                  )
                )
              )
            )
          ),
          
          fluidRow(
            column(12, style = "text-align: center; margin: 30px 0;",
              actionButton("uni_run_analysis", "RUN ANALYSIS", 
                         class = "run-button", style = "width: 80%; max-width: 500px;")
            )
          ),
          
          div(class = "results-section",
            fluidRow(
              column(6,
                box(title = "HISTORICAL DATA SUMMARY", width = 12, status = "success", solidHeader = TRUE,
                  DT::dataTableOutput("uni_hist_summary")
                )
              ),
              column(6,
                box(title = "CURRENT DATA SUMMARY", width = 12, status = "success", solidHeader = TRUE,
                  DT::dataTableOutput("uni_curr_summary")
                )
              )
            ),
            
            fluidRow(
              column(6,
                box(title = "POWER PRIOR PARAMETERS", width = 12, status = "info", solidHeader = TRUE,
                  DT::dataTableOutput("uni_pp_table")
                )
              ),
              column(6,
                box(title = "POSTERIOR PARAMETERS", width = 12, status = "info", solidHeader = TRUE,
                  DT::dataTableOutput("uni_posterior_table")
                )
              )
            ),
            
            fluidRow(
              column(6,
                box(title = "POSTERIOR DISTRIBUTION OF MEAN (mu)", width = 12, status = "warning", solidHeader = TRUE,
                  plotOutput("uni_posterior_mean_plot", height = "400px")
                )
              ),
              column(6,
                box(title = "POSTERIOR DISTRIBUTION OF VARIANCE (sigma2)", width = 12, status = "warning", solidHeader = TRUE,
                  plotOutput("uni_posterior_var_plot", height = "400px")
                )
              )
            )
          )
        )
      ),
      
      # Multivariate Tab ----
      tabItem(tabName = "multivariate",
        div(class = "tab-content-wrapper",
          fluidRow(
            column(12,
              h2("Multivariate Power Prior Analysis", style = "margin-bottom: 30px;")
            )
          ),
          
          fluidRow(
            column(6,
              div(class = "data-upload-section",
                h4("UPLOAD AND SELECT DATA"),
                h5("Historical Data", style = "color: white;"),
                fileInput("multi_hist_file", "Choose CSV File", 
                         accept = ".csv", width = "100%"),
                uiOutput("multi_hist_cols_ui"),
                h5("Current Data", style = "color: white; margin-top: 20px;"),
                fileInput("multi_curr_file", "Choose CSV File", 
                         accept = ".csv", width = "100%"),
                uiOutput("multi_curr_cols_ui"),
                div(style = "background-color: rgba(255,255,255,0.2); padding: 10px; border-radius: 5px; margin-top: 10px;",
                  p(style = "color: white; font-weight: bold; margin: 0;",
                    "IMPORTANT: Select columns in the SAME order from both datasets. Column matching is based on selection order, not names.")
                ),
                fluidRow(
                  column(6, actionButton("multi_load_data", "LOAD DATA", 
                           class = "btn btn-info", style = "width: 100%; margin-top: 10px; color: white; font-weight: bold;")),
                  column(6, actionButton("multi_clear_data", "CLEAR DATA", 
                           class = "btn btn-danger", style = "width: 100%; margin-top: 10px; color: white; font-weight: bold;"))
                )
              )
            ),
            column(6,
              box(title = "CONFIGURATION", width = 12, status = "primary", solidHeader = TRUE,
                div(class = "parameter-section",
                  h5("Discounting Parameter"),
                  sliderInput("multi_a0", "a0 (0 = no borrowing, 1 = full borrowing)", 
                             min = 0, max = 1, value = 0.5, step = 0.05),
                  helpText("Controls how much historical information to use")
                ),
                
                div(class = "parameter-section",
                  h5("Prior Type"),
                  radioButtons("multi_prior_type", NULL,
                              c("Vague (Non-informative)" = "vague",
                                "Informative" = "informative"),
                              selected = "vague", inline = FALSE),
                  
                  conditionalPanel(
                    condition = "input.multi_prior_type == 'informative'",
                    textInput("multi_mu0", "Prior Mean Vector (e.g., 0,0,0)", value = "0,0,0"),
                    numericInput("multi_kappa0", "Prior Precision (kappa0)", value = 1, min = 0.01),
                    numericInput("multi_nu0", "Prior DoF (nu0)", value = 5, min = 1),
                    textInput("multi_lambda0_diag", "Scale Matrix Diagonal (e.g., 1,1,1)", value = "1,1,1"),
                    helpText("Provide diagonal elements for diagonal Lambda0")
                  )
                )
              )
            )
          ),
          
          fluidRow(
            column(12, style = "text-align: center; margin: 30px 0;",
              actionButton("multi_run_analysis", "RUN ANALYSIS", 
                         class = "run-button", style = "width: 80%; max-width: 500px;")
            )
          ),
          
          div(class = "results-section",
            fluidRow(
              column(6,
                box(title = "HISTORICAL DATA SUMMARY", width = 12, status = "success", solidHeader = TRUE,
                  DT::dataTableOutput("multi_hist_summary")
                )
              ),
              column(6,
                box(title = "CURRENT DATA SUMMARY", width = 12, status = "success", solidHeader = TRUE,
                  DT::dataTableOutput("multi_curr_summary")
                )
              )
            ),
            
            fluidRow(
              column(6,
                box(title = "POWER PRIOR PARAMETERS", width = 12, status = "info", solidHeader = TRUE,
                  verbatimTextOutput("multi_pp_summary")
                )
              ),
              column(6,
                box(title = "POSTERIOR PARAMETERS", width = 12, status = "info", solidHeader = TRUE,
                  verbatimTextOutput("multi_posterior_summary")
                )
              )
            ),
            
            fluidRow(
              column(12,
                box(title = "POSTERIOR MEAN DISTRIBUTIONS", width = 12, status = "warning", solidHeader = TRUE,
                  plotOutput("multi_posterior_plot", height = "500px")
                )
              )
            )
          )
        )
      ),
      
      # Help Tab ----
      tabItem(tabName = "help",
        div(class = "tab-content-wrapper",
          fluidRow(
            column(12,
              h2("POWER PRIOR BAYESIAN ANALYSIS - HELP GUIDE"),
              
              h3("Overview"),
              p("Power priors provide a framework for incorporating historical information in Bayesian analysis 
                while controlling the degree of borrowing through a discounting parameter a0."),
              
              h3("Data Selection"),
              p(strong("Univariate Analysis:"), " After uploading your CSV file, select ONE column from your data for analysis."),
              p(strong("Multivariate Analysis:"), " After uploading your CSV files, select the SAME NUMBER of columns from both historical and current data."),
              
              h3("Univariate Analysis"),
              p("Use univariate analysis when your outcome is a single continuous variable with unknown mean and variance."),
              p(strong("Key Parameters:")),
              tags$ul(
                tags$li(strong("a0 (Discounting Parameter):"), " Controls borrowing from historical data [0,1]"),
                tags$li(strong("mu0 (Prior Mean):"), " Prior belief about the center of the distribution"),
                tags$li(strong("kappa0 (Prior Precision):"), " Confidence in prior mean (higher = more confident)"),
                tags$li(strong("nu0 (Prior DoF):"), " Confidence in variance estimate (1-5 for weakly informative)"),
                tags$li(strong("sigma2_0 (Prior Variance):"), " Prior belief about data spread")
              ),
              
              h3("Multivariate Analysis"),
              p("Use multivariate analysis for multiple correlated outcomes (2+ variables)."),
              p(strong("Key Parameters:")),
              tags$ul(
                tags$li(strong("a0:"), " Same as univariate"),
                tags$li(strong("mu0 (Prior Mean Vector):"), " Prior means for each variable"),
                tags$li(strong("kappa0 (Prior Precision):"), " Same precision applied to all variables"),
                tags$li(strong("nu0 (Prior DoF):"), " Must satisfy nu0 >= p (number of variables)"),
                tags$li(strong("Lambda0 (Scale Matrix):"), " Controls prior covariance structure (diagonal for simplicity)")
              ),
              
              h3("Workflow"),
              tags$ol(
                tags$li("Upload CSV file for historical data"),
                tags$li("Select column(s) to use from historical data"),
                tags$li("Upload CSV file for current data"),
                tags$li("Select the same column(s) from current data"),
                tags$li("Click LOAD DATA to import the selected columns"),
                tags$li("Adjust the Discounting Parameter (a0) with the slider"),
                tags$li("Choose between Vague and Informative priors"),
                tags$li("Click RUN ANALYSIS to compute results"),
                tags$li("Review results in tables and visualizations below"),
                tags$li("Adjust parameters and click RUN ANALYSIS again for sensitivity analysis"),
                tags$li("Click CLEAR DATA to remove files and start over")
              ),
              
              h3("References"),
              p("Huang, Y., Yamaguchi, Y., Homma, G., Maruo, K., & Takeda, K. (2024). 
                'Conjugate Representation of Power Priors for Efficient Bayesian Analysis of Normal Data.' 
                Statistical Science (in press)."),
              p("Ibrahim, J. G., & Chen, M. H. (2000). 
                'Power prior distributions for regression models.' 
                Statistical Science, 15(1), 46-60.")
            )
          )
        )
      )
    )
  )
)

# SERVER ----
server <- function(input, output, session) {
  
  # Format numbers to 3 decimal places
  format_num <- function(x) {
    round(as.numeric(x), 3)
  }
  
  # Univariate reactive data ----
  uni_raw_hist <- reactiveVal(NULL)
  uni_raw_curr <- reactiveVal(NULL)
  
  uni_data <- reactiveValues(
    historical = NULL,
    current = NULL,
    powerprior = NULL,
    posterior = NULL
  )
  
  # Univariate: show column selector for historical data
  output$uni_hist_col_ui <- renderUI({
    req(input$uni_hist_file)
    data <- read.csv(input$uni_hist_file$datapath)
    uni_raw_hist(data)
    selectInput("uni_hist_col", "Select Column from Historical Data:", 
               choices = names(data), width = "100%")
  })
  
  # Univariate: show column selector for current data
  output$uni_curr_col_ui <- renderUI({
    req(input$uni_curr_file)
    data <- read.csv(input$uni_curr_file$datapath)
    uni_raw_curr(data)
    selectInput("uni_curr_col", "Select Column from Current Data:", 
               choices = names(data), width = "100%")
  })
  
  observeEvent(input$uni_load_data, {
    req(uni_raw_hist(), uni_raw_curr(), input$uni_hist_col, input$uni_curr_col)
    tryCatch({
      uni_data$historical <- uni_raw_hist()[[input$uni_hist_col]]
      uni_data$current <- uni_raw_curr()[[input$uni_curr_col]]
      showNotification("Success: Data loaded successfully!", type = "message", duration = 3)
    }, error = function(e) {
      showNotification(paste("Error loading data:", e$message), type = "error", duration = 5)
    })
  })
  
  observeEvent(input$uni_clear_data, {
    uni_data$historical <- NULL
    uni_data$current <- NULL
    uni_data$powerprior <- NULL
    uni_data$posterior <- NULL
    uni_raw_hist(NULL)
    uni_raw_curr(NULL)
    shinyjs::reset("uni_hist_file")
    shinyjs::reset("uni_curr_file")
    showNotification("Data and files cleared successfully!", type = "message", duration = 3)
  })
  
  observeEvent(input$uni_run_analysis, {
    req(uni_data$historical, uni_data$current)
    tryCatch({
      if (input$uni_prior_type == "vague") {
        uni_data$powerprior <- powerprior_univariate(
          uni_data$historical,
          a0 = input$uni_a0
        )
      } else {
        uni_data$powerprior <- powerprior_univariate(
          uni_data$historical,
          a0 = input$uni_a0,
          mu0 = input$uni_mu0,
          kappa0 = input$uni_kappa0,
          nu0 = input$uni_nu0,
          sigma2_0 = input$uni_sigma2_0
        )
      }
      
      uni_data$posterior <- posterior_univariate(
        uni_data$powerprior,
        uni_data$current
      )
      showNotification("Success: Analysis complete!", type = "message", duration = 3)
    }, error = function(e) {
      showNotification(paste("Error in analysis:", e$message), type = "error", duration = 5)
    })
  })
  
  output$uni_hist_summary <- DT::renderDataTable({
    req(uni_data$historical)
    summary_df <- data.frame(
      Metric = c("Sample Size", "Mean", "Std Dev", "Min", "Max"),
      Value = c(
        length(uni_data$historical),
        format_num(mean(uni_data$historical)),
        format_num(sd(uni_data$historical)),
        format_num(min(uni_data$historical)),
        format_num(max(uni_data$historical))
      )
    )
    DT::datatable(summary_df, options = list(dom = 't', paging = FALSE, searching = FALSE))
  })
  
  output$uni_pp_table <- DT::renderDataTable({
    req(uni_data$powerprior)
    pp_df <- data.frame(
      Parameter = c("mu_n", "kappa_n", "nu_n", "sigma2_n"),
      Value = c(
        format_num(uni_data$powerprior$mu_n),
        format_num(uni_data$powerprior$kappa_n),
        format_num(uni_data$powerprior$nu_n),
        format_num(uni_data$powerprior$sigma2_n)
      )
    )
    DT::datatable(pp_df, options = list(dom = 't', paging = FALSE, searching = FALSE))
  })
  
  output$uni_curr_summary <- DT::renderDataTable({
    req(uni_data$current)
    summary_df <- data.frame(
      Metric = c("Sample Size", "Mean", "Std Dev", "Min", "Max"),
      Value = c(
        length(uni_data$current),
        format_num(mean(uni_data$current)),
        format_num(sd(uni_data$current)),
        format_num(min(uni_data$current)),
        format_num(max(uni_data$current))
      )
    )
    DT::datatable(summary_df, options = list(dom = 't', paging = FALSE, searching = FALSE))
  })
  
  output$uni_posterior_table <- DT::renderDataTable({
    req(uni_data$posterior)
    post_df <- data.frame(
      Parameter = c("mu*", "kappa*", "nu*", "sigma2*", "Var(mu)"),
      Value = c(
        format_num(uni_data$posterior$mu_star),
        format_num(uni_data$posterior$kappa_star),
        format_num(uni_data$posterior$nu_star),
        format_num(uni_data$posterior$sigma2_star),
        format_num(uni_data$posterior$sigma2_star / uni_data$posterior$kappa_star)
      )
    )
    DT::datatable(post_df, options = list(dom = 't', paging = FALSE, searching = FALSE))
  })
  
  output$uni_posterior_mean_plot <- renderPlot({
    req(uni_data$posterior)
    mu_star <- uni_data$posterior$mu_star
    kappa_star <- uni_data$posterior$kappa_star
    nu_star <- uni_data$posterior$nu_star
    sigma2_star <- uni_data$posterior$sigma2_star
    
    df <- nu_star
    se <- sqrt(sigma2_star / kappa_star)
    x <- seq(mu_star - 4*se, mu_star + 4*se, length.out = 500)
    y <- dt((x - mu_star) / se, df = df) / se
    
    ggplot(data.frame(x = x, y = y), aes(x = x, y = y)) +
      geom_line(color = "steelblue", size = 1) +
      geom_ribbon(aes(ymin = 0, ymax = y), alpha = 0.3, fill = "steelblue") +
      labs(title = "Posterior Distribution of mu",
           x = "mu", y = "Density") +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
  })
  
  output$uni_posterior_var_plot <- renderPlot({
    req(uni_data$posterior)
    nu_star <- uni_data$posterior$nu_star
    sigma2_star <- uni_data$posterior$sigma2_star
    
    x <- seq(sigma2_star * 0.1, sigma2_star * 3, length.out = 500)
    shape <- nu_star / 2
    scale <- 2 * sigma2_star / nu_star
    y <- (1 / (2^(shape/2) * gamma(shape/2))) * (x^(shape/2 - 1)) * exp(-x / scale)
    
    ggplot(data.frame(x = x, y = y), aes(x = x, y = y)) +
      geom_line(color = "coral", size = 1) +
      geom_ribbon(aes(ymin = 0, ymax = y), alpha = 0.3, fill = "coral") +
      labs(title = "Posterior Distribution of sigma2",
           x = "sigma2", y = "Density") +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
  })
  
  # Multivariate reactive data ----
  multi_raw_hist <- reactiveVal(NULL)
  multi_raw_curr <- reactiveVal(NULL)
  
  multi_data <- reactiveValues(
    historical = NULL,
    current = NULL,
    powerprior = NULL,
    posterior = NULL
  )
  
  # Multivariate: show column selector for historical data
  output$multi_hist_cols_ui <- renderUI({
    req(input$multi_hist_file)
    data <- read.csv(input$multi_hist_file$datapath)
    multi_raw_hist(data)
    selectInput("multi_hist_cols", "Select Columns from Historical Data (hold Ctrl/Cmd):", 
               choices = names(data), multiple = TRUE, width = "100%")
  })
  
  # Multivariate: show column selector for current data
  output$multi_curr_cols_ui <- renderUI({
    req(input$multi_curr_file)
    data <- read.csv(input$multi_curr_file$datapath)
    multi_raw_curr(data)
    selectInput("multi_curr_cols", "Select Columns from Current Data (hold Ctrl/Cmd):", 
               choices = names(data), multiple = TRUE, width = "100%")
  })
  
  observeEvent(input$multi_load_data, {
    req(multi_raw_hist(), multi_raw_curr(), input$multi_hist_cols, input$multi_curr_cols)
    tryCatch({
      if (length(input$multi_hist_cols) != length(input$multi_curr_cols)) {
        showNotification("Error: Must select the same number of columns from both datasets!", 
                        type = "error", duration = 5)
        return()
      }
      
      multi_data$historical <- multi_raw_hist()[, input$multi_hist_cols, drop = FALSE]
      multi_data$current <- multi_raw_curr()[, input$multi_curr_cols, drop = FALSE]
      
      showNotification("Success: Data loaded successfully!", type = "message", duration = 3)
    }, error = function(e) {
      showNotification(paste("Error loading data:", e$message), type = "error", duration = 5)
    })
  })
  
  observeEvent(input$multi_clear_data, {
    multi_data$historical <- NULL
    multi_data$current <- NULL
    multi_data$powerprior <- NULL
    multi_data$posterior <- NULL
    multi_raw_hist(NULL)
    multi_raw_curr(NULL)
    shinyjs::reset("multi_hist_file")
    shinyjs::reset("multi_curr_file")
    showNotification("Data and files cleared successfully!", type = "message", duration = 3)
  })
  
  observeEvent(input$multi_run_analysis, {
    req(multi_data$historical, multi_data$current)
    tryCatch({
      p <- ncol(multi_data$historical)
      
      # Debug: print data info
      cat("Historical data:\n")
      print(head(multi_data$historical))
      cat("\nCurrent data:\n")
      print(head(multi_data$current))
      cat("\nHistorical summary:\n")
      print(colMeans(multi_data$historical))
      cat("\nCurrent summary:\n")
      print(colMeans(multi_data$current))
      
      if (input$multi_prior_type == "vague") {
        multi_data$powerprior <- powerprior_multivariate(
          multi_data$historical,
          a0 = input$multi_a0
        )
      } else {
        mu0_vec <- as.numeric(strsplit(input$multi_mu0, ",")[[1]])
        lambda0_diag <- as.numeric(strsplit(input$multi_lambda0_diag, ",")[[1]])
        Lambda0 <- diag(lambda0_diag)
        
        multi_data$powerprior <- powerprior_multivariate(
          multi_data$historical,
          a0 = input$multi_a0,
          mu0 = mu0_vec,
          kappa0 = input$multi_kappa0,
          nu0 = input$multi_nu0,
          Lambda0 = Lambda0
        )
      }
      
      multi_data$posterior <- posterior_multivariate(
        multi_data$powerprior,
        multi_data$current
      )
      showNotification("Success: Analysis complete!", type = "message", duration = 3)
    }, error = function(e) {
      showNotification(paste("Error in analysis:", e$message), type = "error", duration = 5)
    })
  })
  
  output$multi_hist_summary <- DT::renderDataTable({
    req(multi_data$historical)
    summary_df <- data.frame(
      Variable = names(multi_data$historical),
      Mean = format_num(colMeans(multi_data$historical)),
      SD = format_num(apply(multi_data$historical, 2, sd)),
      Min = format_num(apply(multi_data$historical, 2, min)),
      Max = format_num(apply(multi_data$historical, 2, max))
    )
    DT::datatable(summary_df, options = list(dom = 't', paging = FALSE, searching = FALSE))
  })
  
  output$multi_pp_summary <- renderPrint({
    req(multi_data$powerprior)
    cat("Power Prior Parameters (Multivariate)\n")
    cat("=====================================\n\n")
    cat("Mean vector (mu_n):\n")
    print(round(multi_data$powerprior$mu_n, 3))
    cat("\nPrecision (kappa_n):", round(multi_data$powerprior$kappa_n, 3), "\n")
    cat("Degrees of Freedom (nu_n):", round(multi_data$powerprior$nu_n, 3), "\n")
    cat("\nScale Matrix (Lambda_n):\n")
    print(round(multi_data$powerprior$Lambda_n, 3))
  })
  
  output$multi_curr_summary <- DT::renderDataTable({
    req(multi_data$current)
    summary_df <- data.frame(
      Variable = names(multi_data$current),
      Mean = format_num(colMeans(multi_data$current)),
      SD = format_num(apply(multi_data$current, 2, sd)),
      Min = format_num(apply(multi_data$current, 2, min)),
      Max = format_num(apply(multi_data$current, 2, max))
    )
    DT::datatable(summary_df, options = list(dom = 't', paging = FALSE, searching = FALSE))
  })
  
  output$multi_posterior_summary <- renderPrint({
    req(multi_data$posterior)
    cat("Posterior Parameters (Multivariate)\n")
    cat("====================================\n\n")
    cat("Mean vector (mu*):\n")
    print(round(multi_data$posterior$mu_star, 3))
    cat("\nPrecision (kappa*):", round(multi_data$posterior$kappa_star, 3), "\n")
    cat("Degrees of Freedom (nu*):", round(multi_data$posterior$nu_star, 3), "\n")
    cat("\nScale Matrix (Lambda*):\n")
    print(round(multi_data$posterior$Lambda_star, 3))
  })
  
  output$multi_posterior_plot <- renderPlot({
    req(multi_data$posterior)
    
    mu_star <- multi_data$posterior$mu_star
    kappa_star <- multi_data$posterior$kappa_star
    nu_star <- multi_data$posterior$nu_star
    Lambda_star <- multi_data$posterior$Lambda_star
    p <- multi_data$posterior$p
    
    plot_data <- data.frame()
    
    for (i in 1:p) {
      df <- nu_star - p + 1
      se <- sqrt(Lambda_star[i, i] / (kappa_star * df))
      x <- seq(mu_star[i] - 4*se, mu_star[i] + 4*se, length.out = 300)
      y <- dt((x - mu_star[i]) / se, df = df) / se
      
      plot_data <- rbind(plot_data, data.frame(
        x = x,
        y = y,
        Variable = paste("Variable", i)
      ))
    }
    
    ggplot(plot_data, aes(x = x, y = y, fill = Variable, color = Variable)) +
      geom_line(size = 1) +
      geom_ribbon(aes(ymin = 0, ymax = y), alpha = 0.3) +
      facet_wrap(~Variable, scales = "free") +
      labs(title = "Posterior Distributions of Mean Vector",
           x = "Value", y = "Density") +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
            legend.position = "none")
  })
}

# Run app ----
shinyApp(ui, server)