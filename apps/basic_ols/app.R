# app.R — Combined Econometrics Apps (sub-apps in tabs)
# Run: shiny::runApp()

# ---- Packages ----
suppressPackageStartupMessages({
  library(shiny)
  library(ggplot2)
  library(MASS)     # mvrnorm
  library(patchwork)
  library(dplyr)
  library(broom)
})

steel <- "#4682B4"; fire <- "#B22222"

ui <- fluidPage(
  titlePanel("Econometrics Teaching Apps — Combined"),
  tabsetPanel(
    type = "tabs",
    # -------------------------------------------------------------------------
    # 1) Covariance / Variance & Regression
    # -------------------------------------------------------------------------
    tabPanel("Covariance → Regression",
             sidebarLayout(
               sidebarPanel(
                 sliderInput("cv_varx", "Variance of X", min = 0.1, max = 10, value = 1, step = 0.1),
                 numericInput("cv_vary", "Variance of Y", value = 1, min = 0.1, step = 0.1),
                 sliderInput("cv_covxy", "Covariance Cov(X,Y)", min = -10, max = 10, value = 1, step = 0.1),
                 numericInput("cv_n", "Sample size n", value = 100, min = 10),
                 helpText("A valid covariance matrix Σ requires Var(X)>0, Var(Y)>0, and det(Σ)>0.")
               ),
               mainPanel(
                 plotOutput("cv_plot"),
                 verbatimTextOutput("cv_summary"),
                 textOutput("cv_warn")
               )
             )
    ),
    
    # -------------------------------------------------------------------------
    # 2) OLS β1 Convergence (MC) with endogeneity ρ
    # -------------------------------------------------------------------------
    tabPanel("OLS β₁ Convergence (MC)",
             sidebarLayout(
               sidebarPanel(
                 sliderInput("mc_n",  "Sample size n", min = 50, max = 4000, value = 200, step = 50),
                 sliderInput("mc_rho","Corr(x,u) ρ",   min = -0.9, max = 0.9, value = 0,   step = 0.1),
                 numericInput("mc_b1","True β₁", value = 2, step = 0.1),
                 numericInput("mc_R", "Monte Carlo replications", value = 500, min = 100),
                 actionButton("mc_go", "Simulate", class = "btn-primary")
               ),
               mainPanel(
                 plotOutput("mc_hist", height = "420px"),
                 verbatimTextOutput("mc_text")
               )
             )
    ),
    
    # -------------------------------------------------------------------------
    # 3) Interactions & Squared (marginal / total effects)
    # -------------------------------------------------------------------------
    tabPanel("Interactions & Squared",
             sidebarLayout(
               sidebarPanel(
                 h4("True DGP coefficients"),
                 numericInput("ix_b0", "β₀ (intercept)", 0, step = 0.5),
                 numericInput("ix_b1", "β₁ (x effect)", 1, step = 0.5),
                 numericInput("ix_b2", "β₂ (z effect)", 0.5, step = 0.5),
                 numericInput("ix_b3", "β₃ (interaction / curvature)", 0.2, step = 0.1),
                 radioButtons("ix_ztype", "Type of z",
                              c("Continuous" = "cont", "Dummy (0/1)" = "dummy", "z = x (square)" = "square")),
                 conditionalPanel("input.ix_ztype == 'cont'",
                                  sliderInput("ix_zrange", "Range of z (continuous)", -5, 5, c(-3, 3), step = 0.5)
                 ),
                 tags$hr(),
                 h4("Noise / sample"),
                 sliderInput("ix_sigma", "Outcome noise SD σᵤ", min = 0, max = 5, value = 1, step = 0.1),
                 checkboxInput("ix_hetero", "Heteroskedastic σᵤ×(1+|z|)", FALSE),
                 sliderInput("ix_mex", "Measurement error SD in x", min = 0, max = 2, value = 0, step = 0.1),
                 sliderInput("ix_n", "Sample size", 10, 2000, 200, step = 10),
                 numericInput("ix_seed", "Seed", 123, step = 1),
                 actionButton("ix_sim", "Resimulate", class = "btn-secondary")
               ),
               mainPanel(
                 plotOutput("ix_plot", height = "720px"),
                 tags$h4("Regression output (observed vars)"),
                 tableOutput("ix_table")
               )
             )
    ),
    
    # -------------------------------------------------------------------------
    # 4) OVB: omit X2
    # -------------------------------------------------------------------------
    tabPanel("OVB (omit X2)",
             sidebarLayout(
               sidebarPanel(
                 sliderInput("ovb_rho12", "Corr(X1, X2)", min = -1, max = 1, value = 0.5, step = 0.1),
                 sliderInput("ovb_rho2y", "Corr(X2, Y)",  min = -1, max = 1, value = 0.5, step = 0.1),
                 sliderInput("ovb_n",     "Sample size",  min = 50, max = 5000, value = 1000, step = 50),
                 actionButton("ovb_sim", "Simulate", class = "btn-primary")
               ),
               mainPanel(
                 plotOutput("ovb_plot"),
                 verbatimTextOutput("ovb_out")
               )
             )
    ),
    
    # -------------------------------------------------------------------------
    # 5) OVB: 3-variable case
    # -------------------------------------------------------------------------
    tabPanel("OVB (3 variables)",
             sidebarLayout(
               sidebarPanel(
                 sliderInput("ovb3_rho12", "Corr(x1, x2)", min = -1, max = 1, value = 0.5, step = 0.1),
                 sliderInput("ovb3_rho23", "Corr(x2, x3)", min = -1, max = 1, value = 0.5, step = 0.1),
                 sliderInput("ovb3_b1", "β1 (true for x1)", min = -5, max = 5, value = 1, step = 0.1),
                 sliderInput("ovb3_b2", "β2 (true for x2)", min = -5, max = 5, value = 1, step = 0.1),
                 sliderInput("ovb3_b3", "β3 (true for x3)", min = -5, max = 5, value = 1, step = 0.1),
                 numericInput("ovb3_n", "Sample size", value = 1000, min = 10)
               ),
               mainPanel(
                 plotOutput("ovb3_plot"),
                 verbatimTextOutput("ovb3_note")
               )
             )
    ),
    
    # -------------------------------------------------------------------------
    # 6) Attenuation Bias (measurement error)
    # -------------------------------------------------------------------------
    tabPanel("Attenuation Bias",
             sidebarLayout(
               sidebarPanel(
                 numericInput("att_seed", "Seed", value = 123, min = 1, step = 1),
                 sliderInput("att_n", "Sample size", min = 100, max = 20000, value = 1000, step = 100),
                 sliderInput("att_beta", "True slope β", min = -3, max = 3, value = 2, step = 0.1),
                 numericInput("att_alpha", "Intercept α", value = 1, step = 0.1),
                 sliderInput("att_sde", "SD of meas. error σ_e", min = 0, max = 3, value = 0.5, step = 0.05),
                 sliderInput("att_sdu", "SD of outcome noise σ_u", min = 0, max = 3, value = 1, step = 0.05),
                 checkboxInput("att_show_true", "Show true-X scatter/line", TRUE),
                 checkboxInput("att_show_obs",  "Show observed-X scatter/line", TRUE)
               ),
               mainPanel(
                 h4("Coefficient summary"),
                 tableOutput("att_coef"),
                 h4("Reliability & predicted attenuation"),
                 tableOutput("att_rely"),
                 hr(), plotOutput("att_plot", height = 480),
                 hr(), tags$p("Theory: E[β̂_obs] = β × Var(x*) / (Var(x*) + Var(e)).", style = "font-style: italic;")
               )
             )
    )
  )
)

server <- function(input, output, session) {
  
  # ============================================================================
  # 1) Covariance → Regression
  # ============================================================================
  cv_data <- reactive({
    Sigma <- matrix(c(input$cv_varx, input$cv_covxy,
                      input$cv_covxy, input$cv_vary), 2)
    detv <- det(Sigma)
    valid <- (input$cv_varx > 0) && (input$cv_vary > 0) && (detv > 0)
    if (!valid) return(list(valid = FALSE, det = detv))
    
    set.seed(123)
    dat <- as.data.frame(MASS::mvrnorm(n = input$cv_n, mu = c(0, 0), Sigma = Sigma))
    names(dat) <- c("x","y")
    list(valid = TRUE, det = detv, dat = dat)
  })
  
  outputcv <- reactive({
    cv_data()
  })
  
  output$cv_warn <- renderText({
    res <- cv_data()
    if (is.null(res$valid) || res$valid) return("")
    paste0("⚠️ Invalid Σ (not positive definite). Check variances and Cov(X,Y). det(Σ) = ",
           round(res$det, 4))
  })
  
  output$cv_plot <- renderPlot({
    res <- cv_data(); req(res$valid)
    dat <- res$dat
    model <- lm(y ~ x, data = dat)
    ggplot(dat, aes(x, y)) +
      geom_point(color = steel, alpha = 0.6) +
      geom_smooth(method = "lm", se = FALSE, color = fire, linewidth = 1.2) +
      labs(title = "Linear Regression Line",
           subtitle = paste("Slope:", round(coef(model)[2], 3)),
           x = "X", y = "Y") +
      theme_minimal()
  })
  
  output$cv_summary <- renderPrint({
    res <- cv_data(); req(res$valid)
    summary(lm(y ~ x, data = res$dat))
  })
  
  # ============================================================================
  # 2) OLS β1 Convergence (MC)
  # ============================================================================
  mc_res <- eventReactive(input$mc_go, {
    n  <- input$mc_n
    rho <- input$mc_rho
    b1 <- input$mc_b1
    R  <- input$mc_R
    
    betas <- replicate(R, {
      x   <- rnorm(n)
      eps <- rnorm(n)
      u   <- rho * scale(x) + sqrt(1 - rho^2) * eps
      u   <- as.numeric(u) / sd(u)
      y   <- b1 * x + u
      coef(lm(y ~ x))[2]
    })
    list(betas = betas, true = b1)
  })
  
  output$mc_hist <- renderPlot({
    req(mc_res())
    df <- data.frame(beta_hat = mc_res()$betas)
    ggplot(df, aes(beta_hat)) +
      geom_histogram(aes(y = ..density..), bins = 30, fill = "steelblue", alpha = 0.75) +
      geom_vline(xintercept = mc_res()$true, color = "firebrick", linewidth = 1.2) +
      labs(x = expression(hat(beta)[1]), y = "Density",
           title = "Distribution of OLS Estimates of β₁",
           subtitle = paste0("n = ", input$mc_n,
                             ", Corr(x,u) = ", sprintf("%.2f", input$mc_rho),
                             ", True β₁ = ", input$mc_b1)) +
      theme_minimal()
  })
  
  output$mc_text <- renderPrint({
    req(mc_res())
    b <- mc_res()$betas
    cat("Monte Carlo summary of", length(b), "replications:\n\n")
    cat("Mean of β̂₁: ", round(mean(b), 4), "\n")
    cat("SD of β̂₁:   ", round(sd(b), 4), "\n")
    cat("Bias:        ", round(mean(b) - mc_res()$true, 4), "\n")
  })
  
  # ============================================================================
  # 3) Interactions & Squared
  # ============================================================================
  ix_gen <- reactive({
    set.seed(input$ix_seed)
    n <- input$ix_n
    x_true <- runif(n, -5, 5)
    
    if (input$ix_ztype == "cont") {
      z_true <- runif(n, input$ix_zrange[1], input$ix_zrange[2])
      mu <- input$ix_b0 + input$ix_b1 * x_true + input$ix_b2 * z_true + input$ix_b3 * x_true * z_true
    } else if (input$ix_ztype == "dummy") {
      z_true <- rbinom(n, 1, 0.5)
      mu <- input$ix_b0 + input$ix_b1 * x_true + input$ix_b2 * z_true + input$ix_b3 * x_true * z_true
    } else {
      z_true <- x_true
      mu <- input$ix_b0 + (input$ix_b1 + input$ix_b2) * x_true + input$ix_b3 * x_true^2
    }
    
    sigma_vec <- if (isTRUE(input$ix_hetero)) input$ix_sigma * (1 + abs(z_true)) else rep(input$ix_sigma, n)
    y <- mu + rnorm(n, 0, sigma_vec)
    x_obs <- x_true + rnorm(n, 0, input$ix_mex)
    z_obs <- if (input$ix_ztype == "dummy") z_true else z_true
    tibble(y, x = x_obs, z = z_obs, x_true, z_true)
  })
  
  observeEvent(input$ix_sim, { ix_gen() })
  
  ix_fit <- reactive({
    d <- ix_gen(); req(d)
    if (input$ix_ztype == "square") lm(y ~ x + I(x^2), data = d) else lm(y ~ x * z, data = d)
  })
  
  ix_me_x <- reactive({
    b1 <- input$ix_b1; b3 <- input$ix_b3
    if (input$ix_ztype == "cont") {
      z_vals <- seq(input$ix_zrange[1], input$ix_zrange[2], length = 5)
      x_vals <- seq(-5, 5, length = 200)
      do.call(rbind, lapply(z_vals, function(z0)
        data.frame(x = x_vals, me = b1 + b3 * z0, group = paste0("z = ", round(z0, 2)), xlabel = "x")))
    } else if (input$ix_ztype == "dummy") {
      z_vals <- c(0, 1); x_vals <- seq(-5, 5, length = 200)
      do.call(rbind, lapply(z_vals, function(z0)
        data.frame(x = x_vals, me = b1 + b3 * z0, group = paste0("z = ", z0), xlabel = "x")))
    } else {
      x <- seq(-5, 5, length = 200)
      data.frame(x = x, me = input$ix_b1 + input$ix_b2 + 2 * input$ix_b3 * x, group = "∂y/∂x", xlabel = "x")
    }
  })
  
  ix_me_z_ci <- reactive({
    fit <- ix_fit(); V <- vcov(fit); cf <- coef(fit)
    if (input$ix_ztype == "cont") {
      z <- seq(input$ix_zrange[1], input$ix_zrange[2], length = 200)
      b1h <- cf["x"]; b3h <- cf["x:z"]
      var <- V["x","x"] + z^2 * V["x:z","x:z"] + 2 * z * V["x","x:z"]
      data.frame(z, me_hat = b1h + b3h * z, lwr = b1h + b3h * z - 1.96 * sqrt(pmax(var,0)),
                 upr = b1h + b3h * z + 1.96 * sqrt(pmax(var,0)))
    } else if (input$ix_ztype == "dummy") {
      z <- c(0, 1); b1h <- cf["x"]; b3h <- cf["x:z"]
      var <- V["x","x"] + z^2 * V["x:z","x:z"] + 2 * z * V["x","x:z"]
      data.frame(z, me_hat = b1h + b3h * z, lwr = b1h + b3h * z - 1.96 * sqrt(pmax(var,0)),
                 upr = b1h + b3h * z + 1.96 * sqrt(pmax(var,0)))
    } else {
      z <- seq(-5, 5, length = 200)
      b1h <- cf["x"]; b2h <- cf["I(x^2)"]
      var <- V["x","x"] + (2*z)^2 * V["I(x^2)","I(x^2)"] + 2*(2*z) * V["x","I(x^2)"]
      data.frame(z, me_hat = b1h + 2 * b2h * z, lwr = b1h + 2*b2h*z - 1.96 * sqrt(pmax(var,0)),
                 upr = b1h + 2 * b2h * z + 1.96 * sqrt(pmax(var,0)))
    }
  })
  
  ix_pred <- reactive({
    b0 <- input$ix_b0; b1 <- input$ix_b1; b2 <- input$ix_b2; b3 <- input$ix_b3
    x <- seq(-5, 5, length = 200)
    if (input$ix_ztype == "cont") {
      z_vals <- seq(input$ix_zrange[1], input$ix_zrange[2], length = 5)
      df <- expand.grid(x = x, z = z_vals); df$y <- with(df, b0 + b1*x + b2*z + b3*x*z)
      df$group <- factor(round(df$z, 2)); df
    } else if (input$ix_ztype == "dummy") {
      z_vals <- c(0,1); df <- expand.grid(x = x, z = z_vals); df$y <- with(df, b0 + b1*x + b2*z + b3*x*z)
      df$group <- factor(paste0("z = ", df$z)); df
    } else {
      df <- data.frame(x = x); df$y <- b0 + (b1 + b2)*x + b3*x^2; df$group <- "z = x"; df
    }
  })
  
  output$ix_table <- renderTable({
    fit <- ix_fit(); tb <- summary(fit)$coefficients
    if (input$ix_ztype == "square") {
      wanted <- c("(Intercept)", "x", "I(x^2)")
    } else {
      wanted <- c("(Intercept)", "x", "z", "x:z")
    }
    out <- matrix(NA_real_, nrow = length(wanted), ncol = 4,
                  dimnames = list(wanted, colnames(tb)))
    present <- intersect(rownames(tb), wanted)
    out[present, ] <- tb[present, , drop = FALSE]
    round(out, 3)
  }, rownames = TRUE)
  
  output$ix_plot <- renderPlot({
    me_x <- ix_me_x(); me_z <- ix_me_z_ci(); pred <- ix_pred()
    p1 <- ggplot(me_x, aes(x, me, color = group)) +
      geom_line(linewidth = if (input$ix_ztype == "square") 1.2 else 1.1) +
      labs(x = me_x$xlabel[1], y = expression(partialdiff*y/partialdiff*x),
           title = "Marginal Effect of x on y (vs x)") +
      theme_minimal() + theme(legend.title = element_blank())
    p2 <- if (input$ix_ztype == "dummy") {
      ggplot(me_z, aes(z, me_hat)) +
        geom_errorbar(aes(ymin = lwr, ymax = upr), width = 0.1) +
        geom_point(size = 3) + geom_line() + scale_x_continuous(breaks = c(0,1)) +
        labs(x = "z", y = expression(partialdiff*y/partialdiff*x),
             title = "Marginal Effect of x on y (vs z) with 95% CI") + theme_minimal()
    } else {
      ggplot(me_z, aes(z, me_hat)) +
        geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.2) +
        geom_line(linewidth = 1.2) +
        labs(x = if (input$ix_ztype == "square") "z (= x)" else "z",
             y = expression(partialdiff*y/partialdiff*x),
             title = "Marginal Effect of x on y (vs z) with 95% CI") + theme_minimal()
    }
    p3 <- ggplot(pred, aes(x, y, color = group)) +
      geom_line(linewidth = 1.1) +
      labs(x = "x", y = "Predicted y",
           color = if (input$ix_ztype == "square") "z = x" else "z",
           title = "Predicted y vs x (true DGP curves)") + theme_minimal()
    (p1 | p2) / p3 + plot_layout(heights = c(1, 1.2))
  })
  
  # ============================================================================
  # 4) OVB (omit X2)
  # ============================================================================
  ovb_dat <- eventReactive(input$ovb_sim, {
    set.seed(123)
    n <- input$ovb_n; rho12 <- input$ovb_rho12; rho2y <- input$ovb_rho2y
    Sigma <- matrix(c(1, rho12, rho12, 1), 2)
    X <- MASS::mvrnorm(n = n, mu = c(0,0), Sigma = Sigma)
    x1 <- X[,1]; x2 <- X[,2]
    beta1 <- 1; beta2 <- ifelse(abs(rho2y) < 1e-6, 0, rho2y * 1.5)
    y <- beta1 * x1 + beta2 * x2 + rnorm(n)
    mod1 <- lm(y ~ x1); mod2 <- lm(y ~ x1 + x2)
    data.frame(Model = c("Simple (omit X2)","Multiple (with X2)"),
               Estimate = c(coef(mod1)["x1"], coef(mod2)["x1"]))
  })
  
  output$ovb_plot <- renderPlot({
    req(ovb_dat()); df <- ovb_dat()
    ggplot(df, aes(Model, Estimate, fill = Model)) +
      geom_bar(stat = "identity", width = 0.6) +
      geom_hline(yintercept = 1, linetype = "dashed", color = "black") +
      labs(title = "Estimated Coefficient on X1 (omit vs include X2)",
           subtitle = "Dashed line = True effect (1.0)",
           y = "Estimate", x = NULL) +
      theme_minimal() + theme(legend.position = "none")
  })
  
  output$ovb_out <- renderPrint({ req(ovb_dat()); ovb_dat() })
  
  # ============================================================================
  # 5) OVB (3 variables)
  # ============================================================================
  output$ovb3_plot <- renderPlot({
    rho  <- input$ovb3_rho12; rho23 <- input$ovb3_rho23
    Sigma <- matrix(c(1, rho, 0,  rho, 1, rho23, 0, rho23, 1), 3)
    set.seed(123)
    X <- mvrnorm(n = input$ovb3_n, mu = c(0,0,0), Sigma = Sigma)
    x1 <- X[,1]; x2 <- X[,2]; x3 <- X[,3]
    y <- input$ovb3_b1 * x1 + input$ovb3_b2 * x2 + input$ovb3_b3 * x3 + rnorm(input$ovb3_n)
    full <- lm(y ~ x1 + x2 + x3)
    omit2 <- lm(y ~ x1 + x3); omit3 <- lm(y ~ x1 + x2); omit23 <- lm(y ~ x1)
    coefs <- data.frame(
      Model = c("Full (x1,x2,x3)", "Omit x2", "Omit x3", "Omit x2 & x3"),
      Beta1 = c(coef(full)["x1"], coef(omit2)["x1"], coef(omit3)["x1"], coef(omit23)["x1"])
    )
    ggplot(coefs, aes(Model, Beta1)) +
      geom_bar(stat = "identity", fill = "skyblue") +
      geom_hline(yintercept = input$ovb3_b1, linetype = "dashed", color = "red", linewidth = 1) +
      labs(title = "Effect on β1 of Omitting x2/x3",
           y = "Estimated β1", x = NULL) +
      theme_minimal(base_size = 14) + coord_flip()
  })
  output$ovb3_note <- renderText({
    "Dashed red line = true β1. Even when x1 ⟂ x3, omitting x3 can bias β1 via x2."
  })
  
  # ============================================================================
  # 6) Attenuation Bias
  # ============================================================================
  att_dat <- reactive({
    set.seed(input$att_seed)
    n <- input$att_n; beta <- input$att_beta; alpha <- input$att_alpha
    sd_e <- input$att_sde; sd_u <- input$att_sdu
    x_true <- rnorm(n); e <- rnorm(n, 0, sd_e); u <- rnorm(n, 0, sd_u)
    x_obs <- x_true + e; y <- alpha + beta * x_true + u
    tibble(x_true, x_obs, y)
  })
  
  att_fits <- reactive({
    d <- att_dat()
    list(true = lm(y ~ x_true, d), obs = lm(y ~ x_obs, d))
  })
  
  output$att_coef <- renderTable({
    f <- att_fits()
    bind_rows(
      tidy(f$true) %>% mutate(Model = "True (y ~ x_true)"),
      tidy(f$obs)  %>% mutate(Model = "Observed (y ~ x_obs)")
    ) |>
      filter(term != "(Intercept)") |>
      transmute(Model, Term = term,
                Estimate = round(estimate, 3),
                `Std. Error` = round(std.error, 3),
                `t value` = round(statistic, 2),
                `Pr(>|t|)` = signif(p.value, 3))
  })
  
  output$att_rely <- renderTable({
    d <- att_dat()
    var_x_true <- var(d$x_true); var_e <- var(d$x_obs - d$x_true)
    reliability <- var_x_true / (var_x_true + var_e)
    tibble(
      Metric = c("Var(x*)","Var(e)","Reliability Var(x*)/(Var(x*)+Var(e))",
                 "Predicted E[β̂_obs] = β × reliability"),
      Value  = c(round(var_x_true,3), round(var_e,3), round(reliability,3),
                 round(input$att_beta * reliability,3))
    )
  })
  
  output$att_plot <- renderPlot({
    d <- att_dat(); f <- att_fits()
    max_pts <- 4000
    d_plot <- if (nrow(d) > max_pts) d[sample.int(nrow(d), max_pts), ] else d
    g <- ggplot()
    if (isTRUE(input$att_show_obs)) {
      g <- g + geom_point(data = d_plot, aes(x_obs, y), alpha = 0.35) +
        geom_abline(intercept = coef(f$obs)[1], slope = coef(f$obs)[2], color = steel)
    }
    if (isTRUE(input$att_show_true)) {
      g <- g + geom_point(data = d_plot, aes(x_true, y), alpha = 0.25) +
        geom_abline(intercept = coef(f$true)[1], slope = coef(f$true)[2], color = fire)
    }
    g + labs(x = "Regressor (x_true or x_obs)", y = "Outcome (y)",
             title = "Attenuation bias: observed-X slope shrinks toward 0") +
      theme_minimal(base_size = 13)
  })
}

shinyApp(ui, server)