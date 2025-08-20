# app.R — IV Bias & Weak Instruments (teaching app)
# Place this file in apps/iv_bias_weak_instruments/app.R

suppressWarnings({
  if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
  if (!requireNamespace("DT", quietly = TRUE))       install.packages("DT")
})

library(shiny)
library(ggplot2)
library(DT)

fmt4 <- function(x) sprintf("%.4f", x)

# ---- One simulated dataset + diagnostics -------------------------------------
gen_once <- function(N, beta0, beta1, gamma1, rho, seed = NULL, keep_data = FALSE) {
  if (!is.null(seed)) set.seed(seed)
  Z  <- rnorm(N)
  eS <- rnorm(N)
  eU <- rnorm(N)
  S  <- gamma1 * Z + eS
  u  <- rho * eS + sqrt(pmax(1 - rho^2, 0)) * eU
  Y  <- beta0 + beta1 * S + u
  
  # OLS slope
  b_ols <- coef(lm(Y ~ S))[2]
  
  # IV slope (Wald / covariance ratio)
  cov_ZY <- cov(Z, Y)
  cov_ZS <- cov(Z, S)
  b_iv   <- cov_ZY / cov_ZS
  
  # First stage S ~ Z
  fs <- summary(lm(S ~ Z))
  pi1   <- fs$coefficients["Z", "Estimate"]
  tZ    <- fs$coefficients["Z", "t value"]
  F1    <- as.numeric(tZ^2)
  R2_fs <- fs$r.squared
  
  # Reduced form Y ~ Z
  rf <- summary(lm(Y ~ Z))
  delta  <- rf$coefficients["Z", "Estimate"]
  t_rf   <- rf$coefficients["Z", "t value"]
  R2_rf  <- rf$r.squared
  
  out <- list(
    b_ols = b_ols, b_iv = b_iv,
    pi1 = pi1, F1 = F1, R2_fs = R2_fs,
    delta = delta, t_rf = t_rf, R2_rf = R2_rf
  )
  if (keep_data) out$data <- data.frame(Z = Z, S = S, Y = Y)
  out
}

coef_df <- function(fit) {
  sm <- summary(fit)
  cf <- as.data.frame(sm$coefficients)
  cf$Term <- rownames(cf); rownames(cf) <- NULL
  names(cf) <- c("Estimate","Std. Error","t value","Pr(>|t|)","Term")
  cf <- cf[, c("Term","Estimate","Std. Error","t value","Pr(>|t|)")]
  cf[] <- lapply(cf, function(col) if (is.numeric(col)) round(col,4) else col)
  cf
}
stats_df <- function(fit, first_stage = FALSE) {
  sm <- summary(fit)
  r2 <- sm$r.squared
  if (first_stage) {
    tZ <- sm$coefficients["Z", "t value"]; F1 <- as.numeric(tZ^2)
    data.frame(Statistic = c("R^2","First-stage F"), Value = round(c(r2, F1), 4))
  } else {
    data.frame(Statistic = "R^2", Value = round(r2, 4))
  }
}

# ---- UI ----------------------------------------------------------------------
ui <- fluidPage(
  titlePanel("IV Bias & Weak Instruments"),
  sidebarLayout(
    sidebarPanel(
      h4("DGP"),
      numericInput("beta1", "True β₁ (effect of S on Y)", value = 1.00, step = 0.1),
      numericInput("beta0", "True β₀", value = 0.00, step = 0.1),
      
      h4("Controls"),
      sliderInput("rho", "Endogeneity ρ = Corr(e_S, u)", min = -0.95, max = 0.95, value = 0.50, step = 0.05),
      sliderInput("gamma1", "Instrument strength γ₁ in S = γ₁ Z + e_S", min = 0, max = 2.0, value = 0.30, step = 0.05),
      
      h4("Monte Carlo"),
      numericInput("N", "Sample size per replication", value = 500, min = 50, step = 50),
      numericInput("R", "Replications", value = 1000, min = 100, step = 100),
      numericInput("seed", "Random seed", value = 1234, min = 1, step = 1),
      
      actionButton("run", "Run simulation", class = "btn-primary"),
      hr(),
      helpText("Rule-of-thumb: first-stage F ≲ 10 ⇒ weak instrument concerns.")
    ),
    mainPanel(
      fluidRow(
        column(4, strong("Implied R² (first stage): "), textOutput("r2fs_imp", inline = TRUE)),
        column(4, strong("Avg F (MC): "), textOutput("fstat_mc", inline = TRUE)),
        column(4, strong("Share F<10: "), textOutput("share_weak", inline = TRUE))
      ),
      br(),
      h4("Bias and RMSE (Monte Carlo)"),
      DTOutput("summary_tbl"),
      
      h4("Sampling distributions of β̂ (OLS vs IV)"),
      plotOutput("dens_plot", height = 330),
      
      h4("Bias comparison"),
      plotOutput("bias_plot", height = 250),
      
      hr(),
      h4("First Stage & Reduced Form (single dataset)"),
      h5("First Stage: S ~ Z"),
      strong("Coefficients"), DTOutput("fs_coef_tbl"),
      strong("Regression Statistics"), DTOutput("fs_stats_tbl"),
      br(), br(),
      h5("Reduced Form: Y ~ Z"),
      strong("Coefficients"), DTOutput("rf_coef_tbl"),
      strong("Regression Statistics"), DTOutput("rf_stats_tbl"),
      
      hr(),
      h4("Distribution of first-stage F across Monte Carlo"),
      plotOutput("f_hist", height = 260)
    )
  )
)

# ---- Server ------------------------------------------------------------------
server <- function(input, output, session) {
  
  sims <- eventReactive(input$run, {
    set.seed(input$seed)
    mat <- replicate(input$R, unlist(
      gen_once(input$N, input$beta0, input$beta1, input$gamma1, input$rho)
    ))
    as.data.frame(t(mat))
  })
  
  snapshot <- eventReactive(input$run, {
    gen_once(input$N, input$beta0, input$beta1, input$gamma1, input$rho,
             seed = input$seed + 1, keep_data = TRUE)
  })
  
  # Top metrics
  output$r2fs_imp <- renderText({
    # analytical R^2 for S = γ1 Z + eS with Var(Z)=Var(eS)=1
    fmt4(input$gamma1^2 / (input$gamma1^2 + 1))
  })
  output$fstat_mc <- renderText({ df <- sims(); req(df); fmt4(mean(df$F1, na.rm = TRUE)) })
  output$share_weak <- renderText({ df <- sims(); req(df); fmt4(mean(df$F1 < 10, na.rm = TRUE)) })
  
  # Summary table
  output$summary_tbl <- renderDT({
    df <- sims(); req(df)
    b1 <- input$beta1
    out <- data.frame(
      Estimator = c("OLS","IV (Wald/2SLS)"),
      Mean = c(mean(df$b_ols, na.rm = TRUE), mean(df$b_iv, na.rm = TRUE)),
      Bias = c(mean(df$b_ols, na.rm = TRUE) - b1, mean(df$b_iv, na.rm = TRUE) - b1),
      SD   = c(sd(df$b_ols, na.rm = TRUE), sd(df$b_iv, na.rm = TRUE)),
      RMSE = c(sqrt(mean((df$b_ols - b1)^2, na.rm = TRUE)),
               sqrt(mean((df$b_iv  - b1)^2, na.rm = TRUE)))
    )
    out[] <- lapply(out, function(col) if (is.numeric(col)) round(col, 4) else col)
    datatable(out, rownames = FALSE,
              options = list(dom = "t", paging = FALSE, searching = FALSE, ordering = FALSE,
                             class = "compact stripe hover"))
  })
  
  # First stage / reduced form tables (DT)
  output$fs_coef_tbl <- renderDT({
    snap <- snapshot(); req(snap$data); fit <- lm(S ~ Z, data = snap$data)
    datatable(coef_df(fit), rownames = FALSE,
              options = list(dom="t", paging=FALSE, searching=FALSE, ordering=FALSE,
                             class="compact stripe hover"))
  })
  output$fs_stats_tbl <- renderDT({
    snap <- snapshot(); req(snap$data); fit <- lm(S ~ Z, data = snap$data)
    datatable(stats_df(fit, TRUE), rownames = FALSE,
              options = list(dom="t", paging=FALSE, searching=FALSE, ordering=FALSE,
                             class="compact stripe hover"))
  })
  output$rf_coef_tbl <- renderDT({
    snap <- snapshot(); req(snap$data); fit <- lm(Y ~ Z, data = snap$data)
    datatable(coef_df(fit), rownames = FALSE,
              options = list(dom="t", paging=FALSE, searching=FALSE, ordering=FALSE,
                             class="compact stripe hover"))
  })
  output$rf_stats_tbl <- renderDT({
    snap <- snapshot(); req(snap$data); fit <- lm(Y ~ Z, data = snap$data)
    datatable(stats_df(fit, FALSE), rownames = FALSE,
              options = list(dom="t", paging=FALSE, searching=FALSE, ordering=FALSE,
                             class="compact stripe hover"))
  })
  
  # Density plot
  output$dens_plot <- renderPlot({
    df <- sims(); req(df); b1 <- input$beta1
    long <- rbind(
      data.frame(est = df$b_ols, method = "OLS"),
      data.frame(est = df$b_iv,  method = "IV")
    )
    ggplot(long, aes(x = est, color = method, fill = method)) +
      geom_density(alpha = 0.15, linewidth = 1) +
      geom_vline(xintercept = b1, linetype = "dashed") +
      labs(x = expression(hat(beta)), y = "Density",
           subtitle = paste0("ρ = ", input$rho, ",  γ₁ = ", input$gamma1,
                             "  |  N = ", input$N, ",  R = ", input$R)) +
      theme_minimal() +
      theme(legend.title = element_blank())
  })
  
  # Bias bars
  output$bias_plot <- renderPlot({
    df <- sims(); req(df); b1 <- input$beta1
    bias_ols <- mean(df$b_ols - b1, na.rm = TRUE)
    bias_iv  <- mean(df$b_iv  - b1, na.rm = TRUE)
    bars <- data.frame(Method = c("OLS","IV"), Bias = c(bias_ols, bias_iv))
    ggplot(bars, aes(x = Method, y = Bias, fill = Method)) +
      geom_col() + geom_hline(yintercept = 0, linetype = "dashed") +
      labs(y = "Bias", x = NULL) +
      theme_minimal() + theme(legend.position = "none")
  })
  
  # Histogram of first-stage F
  output$f_hist <- renderPlot({
    df <- sims(); req(df)
    ggplot(df, aes(x = F1)) +
      geom_histogram(bins = 30, color = "white") +
      geom_vline(xintercept = 10, linetype = "dashed", linewidth = 1) +
      labs(x = "First-stage F-statistic", y = "Count",
           subtitle = paste0("Mean F = ", fmt4(mean(df$F1, na.rm = TRUE)),
                             " | Share(F<10) = ", fmt4(mean(df$F1 < 10, na.rm = TRUE)))) +
      theme_minimal()
  })
}

shinyApp(ui, server)