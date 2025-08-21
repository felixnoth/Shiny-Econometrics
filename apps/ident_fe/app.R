# app.R — Combined Shiny app (5 sub-apps)
# Includes ONLY the apps you specified:
# 1) ATE/ATT/ATU vs SDO   2) Randomization Inference (CBT)
# 3) Fixed Effects vs Pooled OLS   4) SE simulation (summary/diagnostics/sweep with decay)
# 5) Clustered SEs vs CCV & TSCB (presets + sweeps)

# ---- Packages ----
library(shiny)
library(ggplot2)
library(DT)
library(dplyr)
library(tidyr)
library(purrr)
library(broom)
library(stringr)
library(MASS)
library(sandwich)
library(reshape2)

fmt4 <- function(x) sprintf("%.4f", x)

# =========================================================
# ============== 1) ATE / ATT / ATU vs SDO ================
# =========================================================
mod_ate_ui <- function(id) {
  ns <- NS(id)
  fluidPage(
    titlePanel("Causal Inference Simulation: ATE, SDO, Biases"),
    sidebarLayout(
      sidebarPanel(
        sliderInput(ns("n"), "Sample Size:", min = 10, max = 500, value = 100, step = 10),
        sliderInput(ns("reps"), "Simulations:", min = 100, max = 5000, value = 1000, step = 100),
        sliderInput(ns("treat_prop"), "Treatment Proportion:", min = 0.1, max = 0.9, value = 0.5, step = 0.1),
        sliderInput(ns("hetero"), "Heterogeneity (SD of Treatment Effect):", min = 0, max = 5, value = 1, step = 0.1),
        checkboxInput(ns("selection"), "Enable Selection Bias", value = FALSE)
      ),
      mainPanel(
        plotOutput(ns("histPlot")),
        verbatimTextOutput(ns("summaryStats"))
      )
    )
  )
}
mod_ate_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    simulateData <- reactive({
      n <- input$n
      reps <- input$reps
      prop_treat <- input$treat_prop
      hetero_sd <- input$hetero
      selection <- input$selection
      
      results <- replicate(reps, {
        y0  <- rnorm(n, mean = 5, sd = 2)
        tau <- rnorm(n, mean = 1, sd = hetero_sd)
        y1  <- y0 + tau
        D   <- if (selection) rbinom(n, 1, plogis(tau)) else rbinom(n, 1, prop_treat)
        Y   <- ifelse(D == 1, y1, y0)
        ATE <- mean(y1 - y0)
        ATT <- mean(y1[D == 1] - y0[D == 1])
        ATU <- mean(y1[D == 0] - y0[D == 0])
        SDO <- mean(Y[D == 1]) - mean(Y[D == 0])
        c(ATE = ATE, ATT = ATT, ATU = ATU, SDO = SDO)
      })
      as.data.frame(t(results))
    })
    
    output$histPlot <- renderPlot({
      df <- simulateData()
      ggplot(df, aes(x = SDO)) +
        geom_histogram(bins = 50, fill = "skyblue", color = "black") +
        geom_vline(aes(xintercept = mean(ATE)),
                   color = "firebrick", linetype = "dashed", linewidth = 1) +
        labs(title = "Distribution of Simple Difference in Means (SDO)",
             subtitle = "Dashed line = True ATE",
             x = "SDO Estimate", y = "Frequency") +
        theme_minimal()
    })
    output$summaryStats <- renderPrint({ summary(simulateData()) })
  })
}

# =========================================================
# ===== 2) Randomization Inference (Fisher — CBT) =========
# =========================================================
ri_base_df <- tibble::tibble(
  name = c("Andy","Ben","Chad","Daniel","Edith","Frank","George","Hank"),
  D_obs = c(1,1,1,1,0,0,0,0),
  Y_obs = c(10,5,16,3,5,7,8,10)
)
ri_stat_sdo <- function(y, d, absolute = TRUE){
  sdo <- mean(y[d==1]) - mean(y[d==0]); if (absolute) abs(sdo) else sdo
}
ri_all_complete_assignments <- function(n = 8, m = 4){
  treated_sets <- combn(n, m); asplit(treated_sets, 2)
}
mod_ri_ui <- function(id) {
  ns <- NS(id)
  fluidPage(
    titlePanel("Randomization Inference (Fisher’s Sharp Null) — CBT Example"),
    sidebarLayout(
      sidebarPanel(
        radioButtons(ns("design"), "Assignment scheme for the null:",
                     choices = c("Complete randomization (exact p-value)" = "complete",
                                 "Bernoulli (Monte Carlo p-value)" = "bernoulli"),
                     selected = "complete"),
        checkboxInput(ns("abs_stat"), "Use |Ȳ_T - Ȳ_C|", TRUE),
        conditionalPanel(
          condition = sprintf("input['%s'] == 'bernoulli'", ns("design")),
          sliderInput(ns("p_treat"), "Bernoulli Pr(Treat)", min = 0, max = 1, value = 0.5, step = 0.05),
          numericInput(ns("B"), "Monte Carlo permutations (B)", value = 5000, min = 100, step = 100)
        ),
        tags$hr(),
        actionButton(ns("shuffle_obs"), "Shuffle observed assignment", class = "btn-primary"),
        helpText("Note: Shuffling always keeps 4 treated. It only changes the 'observed' assignment for teaching.")
      ),
      mainPanel(
        fluidRow(
          column(6, tags$h4("Observed data & assignment"), tableOutput(ns("obs_table"))),
          column(6, tags$h4("Observed test statistic"), verbatimTextOutput(ns("obs_stat")))
        ),
        tags$hr(),
        fluidRow(
          column(7, plotOutput(ns("hist_plot"), height = 360)),
          column(5, tags$h4("Permutation p-value"), verbatimTextOutput(ns("pval_text")),
                 tags$small("p-value = share of permutation statistics ≥ observed."),
                 tags$br(), tags$br(), tags$h5("First few permutations"), tableOutput(ns("perm_head")))
        )
      )
    )
  )
}
mod_ri_server <- function(id) {
  moduleServer(id, function(input, output, session){
    obs <- reactiveVal(ri_base_df)
    observeEvent(input$shuffle_obs, {
      set.seed(sample.int(1e6, 1)); idx <- sample(1:nrow(ri_base_df), 4)
      obs(ri_base_df %>% mutate(D_obs = as.integer(row_number() %in% idx)))
    })
    output$obs_table <- renderTable({ obs() %>% transmute(Name = name, Treated = D_obs, Y_obs) })
    output$obs_stat <- renderPrint({
      d <- obs(); s <- ri_stat_sdo(d$Y_obs, d$D_obs, absolute = isTRUE(input$abs_stat))
      cat(sprintf("%s = %.3f", if (isTRUE(input$abs_stat)) "|SDO|" else "SDO", s))
    })
    perm_results <- reactive({
      d <- obs(); y_null <- d$Y_obs
      if (input$design == "complete"){
        sets <- ri_all_complete_assignments(n = nrow(d), m = sum(d$D_obs))
        stats <- vapply(sets, function(idx){
          D <- integer(nrow(d)); D[idx] <- 1L
          ri_stat_sdo(y_null, D, absolute = isTRUE(input$abs_stat))
        }, numeric(1))
        T_obs <- ri_stat_sdo(y_null, d$D_obs, absolute = isTRUE(input$abs_stat))
        tibble::tibble(T_stat = stats, kind = "perm") %>%
          bind_rows(tibble::tibble(T_stat = T_obs, kind = "observed"))
      } else {
        set.seed(123); B <- as.integer(input$B); p <- input$p_treat
        stats <- replicate(B, {
          D <- rbinom(nrow(d), 1, p); while (length(unique(D)) == 1) D <- rbinom(nrow(d), 1, p)
          ri_stat_sdo(y_null, D, absolute = isTRUE(input$abs_stat))
        })
        T_obs <- ri_stat_sdo(y_null, d$D_obs, absolute = isTRUE(input$abs_stat))
        tibble::tibble(T_stat = stats, kind = "perm") %>%
          bind_rows(tibble::tibble(T_stat = T_obs, kind = "observed"))
      }
    })
    output$pval_text <- renderPrint({
      df <- perm_results(); T_obs <- df %>% filter(kind=="observed") %>% pull(T_stat)
      T_perm <- df %>% filter(kind=="perm") %>% pull(T_stat)
      cat(sprintf("Permutation p-value = %.4f", mean(T_perm >= T_obs)))
    })
    output$perm_head <- renderTable({
      df <- perm_results(); df %>% filter(kind=="perm") %>% head(8) %>%
        transmute(`Test statistic` = round(T_stat, 3))
    })
    output$hist_plot <- renderPlot({
      df <- perm_results(); T_obs <- df %>% filter(kind=="observed") %>% pull(T_stat)
      perms <- df %>% filter(kind=="perm")
      ggplot(perms, aes(x = T_stat)) +
        geom_histogram(bins = 20, fill = "grey70", color = "grey30") +
        geom_vline(xintercept = T_obs, linewidth = 1.1, color = "#B22222") +
        labs(x = "Permutation test statistic", y = "Count",
             title = "Randomization distribution under Fisher’s sharp null",
             subtitle = sprintf("Observed statistic = %.3f", T_obs)) +
        theme_minimal(base_size = 13)
    })
  })
}

# =========================================================
# === 3) Fixed Effects vs Pooled OLS for the slope β1 =====
# =========================================================
gen_panel_fe <- function(n_firms, T, beta0, beta1, sigma_a, sigma_g, s_x_firm, s_x_year, sigma_e = 1) {
  tot <- s_x_firm + s_x_year
  if (tot > 0.999) { s_x_firm <- s_x_firm / tot * 0.999; s_x_year <- s_x_year / tot * 0.999 }
  s_x_idio <- 1 - s_x_firm - s_x_year
  firmid <- rep(seq_len(n_firms), each = T); year <- rep(seq_len(T), times = n_firms)
  a_i <- rnorm(n_firms); g_t <- rnorm(T); a <- a_i[firmid]; g <- g_t[year]
  u <- rnorm(n_firms * T)
  x <- sqrt(s_x_idio) * u + sqrt(s_x_firm) * a + sqrt(s_x_year) * g
  e <- rnorm(n_firms * T, sd = sigma_e)
  y <- beta0 + beta1 * x + sigma_a * a + sigma_g * g + e
  data.frame(y = y, x = x, firm = factor(firmid), year = factor(year))
}
est_all_fe <- function(dat) {
  c(
    `Pooled OLS` = coef(lm(y ~ x, data = dat))[["x"]],
    `Firm FE`    = coef(lm(y ~ x + firm, data = dat))[["x"]],
    `Year FE`    = coef(lm(y ~ x + year, data = dat))[["x"]],
    `Two-way FE` = coef(lm(y ~ x + firm + year, data = dat))[["x"]]
  )
}
summarize_model_fe <- function(bhats, beta1) {
  m  <- mean(bhats); sd <- sd(bhats); bias <- m - beta1; rmse <- sqrt(mean((bhats - beta1)^2))
  c(mean = m, bias = bias, sd = sd, rmse = rmse)
}
mod_fe_ui <- function(id) {
  ns <- NS(id)
  fluidPage(
    titlePanel("Fixed Effects and the Estimated Slope (β₁)"),
    sidebarLayout(
      sidebarPanel(
        h4("Data generating process"),
        numericInput(ns("beta1"), "True β₁", value = 1.0, step = 0.1),
        numericInput(ns("beta0"), "True β₀ (intercept)", value = 0.0, step = 0.1),
        sliderInput(ns("sigma_a"), "Strength of firm effect in y (σ_a)", min = 0, max = 3, value = 1, step = 0.1),
        sliderInput(ns("sigma_g"), "Strength of year effect in y (σ_g)", min = 0, max = 3, value = 1, step = 0.1),
        sliderInput(ns("s_x_firm"), "Share of x from firm FE", min = 0, max = 0.95, value = 0.3, step = 0.01),
        sliderInput(ns("s_x_year"), "Share of x from year FE", min = 0, max = 0.95, value = 0.3, step = 0.01),
        sliderInput(ns("sigma_e"), "Noise sd in y (σ_e)", min = 0.1, max = 3, value = 1, step = 0.1),
        hr(),
        h4("Panel & Monte Carlo"),
        numericInput(ns("n_firms"), "Number of firms", value = 50, min = 5, step = 5),
        numericInput(ns("T"), "Years per firm", value = 20, min = 2, step = 1),
        numericInput(ns("n_iter"), "Monte Carlo iterations", value = 1000, min = 100, step = 100),
        numericInput(ns("seed"), "Random seed", value = 1234, min = 1, step = 1),
        actionButton(ns("run"), "Run Simulation", class = "btn-primary")
      ),
      mainPanel(
        h4("Summary across iterations"),
        DTOutput(ns("summary_table")),
        br(),
        h4("Sampling distributions of β̂₁"),
        plotOutput(ns("dens_plot"), height = 380),
        p("Pooled is biased when x correlates with unobserved firm/year effects.
          One-way FE removes bias in the corresponding dimension; two-way removes both.")
      )
    )
  )
}
mod_fe_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    results <- eventReactive(input$run, {
      set.seed(input$seed); n <- input$n_iter
      bhats <- matrix(NA_real_, nrow = n, ncol = 4)
      colnames(bhats) <- c("Pooled OLS","Firm FE","Year FE","Two-way FE")
      for (j in seq_len(n)) {
        dat <- gen_panel_fe(input$n_firms, input$T, input$beta0, input$beta1,
                            input$sigma_a, input$sigma_g, input$s_x_firm, input$s_x_year, input$sigma_e)
        bhats[j, ] <- est_all_fe(dat)
      }
      sum_list <- lapply(seq_len(ncol(bhats)), function(k) summarize_model_fe(bhats[,k], beta1 = input$beta1))
      list(bhats = bhats, summary = do.call(rbind, sum_list))
    })
    output$summary_table <- renderDT({
      req(results()); tab <- results()$summary
      df <- data.frame(
        Model = rownames(tab),
        Mean = fmt4(tab[, "mean"]), Bias = fmt4(tab[, "bias"]),
        SD = fmt4(tab[, "sd"]), RMSE = fmt4(tab[, "rmse"]), row.names = NULL
      )
      datatable(df, rownames = FALSE, options = list(dom = "t", pageLength = 10))
    })
    output$dens_plot <- renderPlot({
      req(results()); bh <- as.data.frame(results()$bhats)
      bh_long <- melt(bh, variable.name = "Model", value.name = "beta_hat")
      ggplot(bh_long, aes(x = beta_hat, color = Model, fill = Model)) +
        geom_density(alpha = 0.12, linewidth = 0.9) +
        geom_vline(xintercept = input$beta1, linetype = "dashed") +
        labs(x = expression(hat(beta)[1]), y = "Density")
    })
  })
}

# ============================================================================
# 4) Summary-only SE simulation + diagnostics + bias sweep (with decays)
# ============================================================================
compute_cluster_se <- function(fit, data, cluster_by = c("firm","year","both")) {
  cluster_by <- match.arg(cluster_by)
  if (cluster_by == "firm") {
    vc <- sandwich::vcovCL(fit, cluster = data$firmid, type = "HC1")
  } else if (cluster_by == "year") {
    vc <- sandwich::vcovCL(fit, cluster = data$year, type = "HC1")
  } else {
    if (requireNamespace("multiwayvcov", quietly = TRUE)) {
      vc <- multiwayvcov::cluster.vcov(fit, cbind(data$firmid, data$year))
    } else {
      vc <- sandwich::vcovCL(fit, cluster = data$firmid, type = "HC1")
      attr(vc, "fallback_note") <- "two-way clustering not available; fell back to firm clustering"
    }
  }
  se_x <- sqrt(diag(vc))[["x"]]
  attr(se_x, "note") <- attr(vc, "fallback_note")
  se_x
}
build_clustered_decay <- function(n_firms, T, s, mode, decay_firm = 0, decay_year = 0) {
  N <- n_firms * T; firmid <- rep(seq_len(n_firms), each = T); year <- rep(seq_len(T), times = n_firms)
  norm_to_var <- function(v, target_s) {
    if (target_s <= 0) return(rep(0, length(v)))
    sv <- stats::sd(v); if (is.na(sv) || sv == 0) return(rep(0, length(v)))
    v * sqrt(target_s) / sv
  }
  s_firm <- 0; s_year <- 0
  if (mode == "firm") s_firm <- s else if (mode == "year") s_year <- s else if (mode == "both") { s_firm <- s/2; s_year <- s/2 }
  s_idio <- max(0, 1 - (s_firm + s_year))
  v_firm_raw <- numeric(N)
  if (s_firm > 0) {
    z_firm <- rnorm(n_firms)
    w_years <- (1 - decay_firm)^(seq_len(T) - 1)
    w_years <- w_years / sqrt(mean(w_years^2))
    v_firm_raw <- z_firm[firmid] * w_years[year]
    v_firm_raw <- norm_to_var(v_firm_raw, s_firm)
  }
  v_year_raw <- numeric(N)
  if (s_year > 0) {
    z_year <- rnorm(T)
    w_firms <- (1 - decay_year)^(seq_len(n_firms) - 1)
    w_firms <- w_firms / sqrt(mean(w_firms^2))
    v_year_raw <- z_year[year] * w_firms[firmid]
    v_year_raw <- norm_to_var(v_year_raw, s_year)
  }
  v_idio <- if (s_idio > 0) rnorm(N) else numeric(N)
  v_idio <- norm_to_var(v_idio, s_idio)
  v_firm_raw + v_year_raw + v_idio
}
gen_dataset <- function(n_firms, T, share_x, share_r, x_mode, r_mode,
                        decay_firm_x, decay_year_x, decay_firm_r, decay_year_r) {
  firmid <- factor(rep(seq_len(n_firms), each = T))
  year   <- factor(rep(seq_len(T), times = n_firms))
  x <- build_clustered_decay(n_firms, T, s = share_x, mode = x_mode,
                             decay_firm = decay_firm_x, decay_year = decay_year_x)
  r_inner <- build_clustered_decay(n_firms, T, s = share_r, mode = r_mode,
                                   decay_firm = decay_firm_r, decay_year = decay_year_r)
  r <- 2 * r_inner
  y <- 1 * x + r
  data.frame(y = y, x = x, firmid = firmid, year = year)
}
fit_once <- function(dat, fe = c("none","firm","year","twoway"), cluster_by = c("firm","year","both")) {
  fe <- match.arg(fe); cluster_by <- match.arg(cluster_by)
  fmla <- switch(fe,
                 "none"   = y ~ x,
                 "firm"   = y ~ x + firmid,
                 "year"   = y ~ x + year,
                 "twoway" = y ~ x + firmid + year)
  fit <- lm(fmla, data = dat)
  b   <- coef(fit)[["x"]]
  se_classical <- sqrt(diag(vcov(fit)))[["x"]]
  se_cluster   <- compute_cluster_se(fit, dat, cluster_by = cluster_by)
  list(b = b, se_classical = se_classical, se_cluster = se_cluster)
}
mod_se_summary_ui <- function(id) {
  ns <- NS(id)
  fluidPage(
    titlePanel("Summary SE Simulation (decay, FE options) + Diagnostics + Bias Sweep"),
    sidebarLayout(
      sidebarPanel(
        h4("DGP: shares & decay"),
        sliderInput(ns("share_x"), "Share of X variance clustered", min = 0, max = 1, value = 0.5, step = 0.01),
        sliderInput(ns("share_r"), "Share of residual variance clustered", min = 0, max = 1, value = 0.5, step = 0.01),
        selectInput(ns("x_mode"), "X clustering mode", choices = c("none","firm","year","both"), selected = "firm"),
        selectInput(ns("r_mode"), "Residual clustering mode", choices = c("none","firm","year","both"), selected = "firm"),
        sliderInput(ns("decay_firm_x"), "Decay of firm component in X", min = 0, max = 1, value = 0, step = 0.05),
        sliderInput(ns("decay_year_x"), "Decay of year component in X", min = 0, max = 1, value = 0, step = 0.05),
        sliderInput(ns("decay_firm_r"), "Decay of firm component in residual", min = 0, max = 1, value = 0, step = 0.05),
        sliderInput(ns("decay_year_r"), "Decay of year component in residual", min = 0, max = 1, value = 0, step = 0.05),
        
        h4("Panel & simulation"),
        numericInput(ns("n_firms"), "Number of firms (clusters)", value = 50, min = 5, step = 5),
        numericInput(ns("T"), "Years per firm", value = 20, min = 2, step = 1),
        numericInput(ns("n_iter"), "Iterations", value = 1000, min = 10, step = 50),
        numericInput(ns("seed"), "Random seed", value = 1234, min = 1, step = 1),
        
        h4("Regression & SE (summary run)"),
        selectInput(ns("fe_mode"), "Fixed effects", choices = c("none","firm","year","twoway"), selected = "none"),
        selectInput(ns("se_cluster_by"), "SE clustering (summary)", choices = c("firm","year","both"), selected = "firm"),
        actionButton(ns("run"), "Run Simulation", class = "btn-primary"),
        hr(),
        
        h4("Diagnostic panel"),
        numericInput(ns("diag_iter"), "Diagnostic iterations", value = 2000, min = 200, step = 100),
        actionButton(ns("run_diag"), "Run Diagnostic"),
        hr(),
        
        h4("Bias sweep (clusters or years)"),
        selectInput(ns("sweep_dim"), "Sweep dimension", choices = c("firms","years"), selected = "firms"),
        numericInput(ns("sweep_min"), "Sweep min", value = 10, min = 5, step = 5),
        numericInput(ns("sweep_max"), "Sweep max", value = 200, min = 10, step = 10),
        numericInput(ns("sweep_step"), "Sweep step", value = 10, min = 5, step = 5),
        numericInput(ns("sweep_iter"), "Iterations per point", value = 300, min = 50, step = 50),
        selectInput(ns("sweep_fe"), "FE in sweep", choices = c("none","firm","year","twoway"), selected = "none"),
        selectInput(ns("sweep_cluster_by"), "SE clustering in sweep", choices = c("firm","year","both"), selected = "firm"),
        
        checkboxInput(ns("sweep_keepN"), "Keep total N approximately fixed (balanced panel)", value = TRUE),
        numericInput(ns("sweep_N"), "Target total N (for sweep)", value = 10000, min = 100, step = 100),
        
        h5("Sweep-specific decay"),
        sliderInput(ns("sweep_decay_firm_x"), "Decay of firm component in X (sweep)", min = 0, max = 1, value = 0, step = 0.05),
        sliderInput(ns("sweep_decay_year_x"), "Decay of year component in X (sweep)", min = 0, max = 1, value = 0, step = 0.05),
        sliderInput(ns("sweep_decay_firm_r"), "Decay of firm component in residual (sweep)", min = 0, max = 1, value = 0, step = 0.05),
        sliderInput(ns("sweep_decay_year_r"), "Decay of year component in residual (sweep)", min = 0, max = 1, value = 0, step = 0.05),
        
        actionButton(ns("run_sweep"), "Run Sweep")
      ),
      mainPanel(
        h4("Summary (across iterations)"),
        DTOutput(ns("summary_table")),
        hr(),
        h4("Diagnostic: mean(SE) vs sd(β̂)"),
        DTOutput(ns("diag_table")),
        hr(),
        h4("SE Bias vs. Number of Clusters / Years"),
        p("Bias = mean(SE) − SD(β̂). Two curves: classical vs clustered (your chosen clustering)."),
        plotOutput(ns("bias_plot"), height = 380)
      )
    )
  )
}
mod_se_summary_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    
    rv <- reactiveValues(summary = NULL, diag = NULL, bias_df = NULL, sweep_note = NULL)
    
    observeEvent(input$run, {
      set.seed(input$seed); n_iter <- input$n_iter
      b  <- numeric(n_iter); seC <- numeric(n_iter); seK <- numeric(n_iter)
      withProgress(message = "Simulating (summary)...", value = 0, {
        inc <- 1 / n_iter
        for (j in seq_len(n_iter)) {
          dat <- gen_dataset(
            n_firms = input$n_firms, T = input$T,
            share_x = input$share_x, share_r = input$share_r,
            x_mode = input$x_mode, r_mode = input$r_mode,
            decay_firm_x = input$decay_firm_x, decay_year_x = input$decay_year_x,
            decay_firm_r = input$decay_firm_r, decay_year_r = input$decay_year_r
          )
          fit <- fit_once(dat, fe = input$fe_mode, cluster_by = input$se_cluster_by)
          b[j]   <- fit$b; seC[j] <- fit$se_classical; seK[j] <- fit$se_cluster
          incProgress(inc)
        }
      })
      out <- data.frame(
        Metric = c("Avg(β̂)", "SD(β̂)", "Avg(SE classical)", paste0("Avg(SE cluster, ", input$se_cluster_by, ")")),
        Value  = c(mean(b), sd(b), mean(seC), mean(seK))
      )
      out$Value <- fmt4(out$Value); rv$summary <- out
    })
    output$summary_table <- renderDT({
      req(rv$summary); datatable(rv$summary, rownames = FALSE, options = list(dom = "t", pageLength = 10))
    })
    
    observeEvent(input$run_diag, {
      set.seed(input$seed + 777); n_iter <- input$diag_iter
      b <- seC <- seK <- numeric(n_iter)
      withProgress(message = "Running diagnostic...", value = 0, {
        inc <- 1 / n_iter
        for (j in seq_len(n_iter)) {
          dat <- gen_dataset(
            n_firms = input$n_firms, T = input$T,
            share_x = input$share_x, share_r = input$share_r,
            x_mode = input$x_mode, r_mode = input$r_mode,
            decay_firm_x = input$decay_firm_x, decay_year_x = input$decay_year_x,
            decay_firm_r = input$decay_firm_r, decay_year_r = input$decay_year_r
          )
          fit <- fit_once(dat, fe = input$fe_mode, cluster_by = input$se_cluster_by)
          b[j] <- fit$b; seC[j] <- fit$se_classical; seK[j] <- fit$se_cluster
          incProgress(inc)
        }
      })
      sd_b <- sd(b)
      diag_df <- data.frame(
        Metric = c("sd(β̂)", "mean(SE classical)", "mean(SE cluster)",
                   "mean(SE classical) / sd(β̂)", "mean(SE cluster) / sd(β̂)"),
        Value  = c(sd_b, mean(seC), mean(seK), mean(seC)/sd_b, mean(seK)/sd_b)
      )
      diag_df$Value <- fmt4(diag_df$Value); rv$diag <- diag_df
    })
    output$diag_table <- renderDT({
      req(rv$diag); datatable(rv$diag, rownames = FALSE, options = list(dom = "t", pageLength = 10))
    })
    
    observeEvent(input$run_sweep, {
      set.seed(input$seed + 999)
      grid <- seq(input$sweep_min, input$sweep_max, by = input$sweep_step)
      bias_classical <- bias_cluster <- numeric(length(grid)); effN <- integer(length(grid)); note <- NULL
      withProgress(message = "Running sweep...", value = 0, {
        inc <- 1 / length(grid)
        for (i in seq_along(grid)) {
          if (isTRUE(input$sweep_keepN)) {
            if (input$sweep_dim == "firms") {
              n_firms_i <- grid[i]; T_i <- max(2L, as.integer(round(input$sweep_N / n_firms_i)))
            } else {
              T_i <- grid[i]; n_firms_i <- max(2L, as.integer(round(input$sweep_N / T_i)))
            }
          } else {
            if (input$sweep_dim == "firms") { n_firms_i <- grid[i]; T_i <- input$T }
            else { n_firms_i <- input$n_firms; T_i <- grid[i] }
          }
          effN[i] <- n_firms_i * T_i
          b <- seC <- seK <- numeric(input$sweep_iter)
          for (j in seq_len(input$sweep_iter)) {
            dat <- gen_dataset(
              n_firms = n_firms_i, T = T_i,
              share_x = input$share_x, share_r = input$share_r,
              x_mode = input$x_mode, r_mode = input$r_mode,
              decay_firm_x = input$sweep_decay_firm_x, decay_year_x = input$sweep_decay_year_x,
              decay_firm_r = input$sweep_decay_firm_r, decay_year_r = input$sweep_decay_year_r
            )
            fit <- fit_once(dat, fe = input$sweep_fe, cluster_by = input$sweep_cluster_by)
            b[j] <- fit$b; seC[j] <- fit$se_classical; if (is.null(note)) note <- attr(fit$se_cluster, "note"); seK[j] <- fit$se_cluster
          }
          sd_b <- sd(b)
          bias_classical[i] <- mean(seC) - sd_b
          bias_cluster[i]   <- mean(seK) - sd_b
          incProgress(inc, detail = paste0(input$sweep_dim, " = ", grid[i], ", effN = ", effN[i]))
        }
      })
      rv$bias_df <- data.frame(x = grid, classical = bias_classical, cluster = bias_cluster, effN = effN)
      rv$sweep_note <- note
    })
    output$bias_plot <- renderPlot({
      req(rv$bias_df)
      df <- rv$bias_df
      xlab <- if (input$sweep_dim == "firms") "Number of firms" else "Number of years"
      clabel <- if (input$sweep_cluster_by == "both") "Clustered SE (two-way)" else paste0("Clustered SE (", input$sweep_cluster_by, ")")
      subtitle_parts <- c(
        if (!is.null(rv$sweep_note)) rv$sweep_note,
        if (isTRUE(input$sweep_keepN)) paste0("Target N=", input$sweep_N, " | Effective N: min=", min(df$effN), ", max=", max(df$effN)) else NULL,
        paste0("Decay X: firm=", input$sweep_decay_firm_x, ", year=", input$sweep_decay_year_x),
        paste0("Decay r: firm=", input$sweep_decay_firm_r, ", year=", input$sweep_decay_year_r)
      )
      subtitle <- paste(subtitle_parts[!sapply(subtitle_parts, is.null)], collapse = " | ")
      df_long <- data.frame(
        x = rep(df$x, 2),
        bias = c(df$classical, df$cluster),
        type = factor(rep(c("Classical SE", clabel), each = nrow(df)))
      )
      ggplot(df_long, aes(x = x, y = bias, linetype = type)) +
        geom_line() + geom_point() +
        geom_hline(yintercept = 0, linetype = "dashed") +
        labs(title = paste0("SE Bias vs ", if (input$sweep_dim=="firms") "# Firms" else "# Years", " — FE: ", input$sweep_fe),
             subtitle = subtitle, x = xlab, y = "Bias = mean(SE) − SD(β̂)")
    })
  })
}

# ============================================================================
# 5) Clustered SEs vs CCV & TSCB (presets + side-by-side sweeps)
# ============================================================================
needs_sandwich <- function() {
  if (!requireNamespace("sandwich", quietly = TRUE)) {
    stop("Package 'sandwich' is required. Install with install.packages('sandwich').")
  }
}
rbeta_from_mean_var <- function(mean, var) {
  mean <- min(max(mean, 1e-3), 1 - 1e-3)
  var  <- max(min(var, mean*(1-mean) - 1e-6), 1e-6)
  k <- mean*(1-mean)/var - 1
  c(a = mean * k, b = (1-mean) * k)
}
simulate_population <- function(K, N_per_cluster, p_bar, assign_var, tau_mean, tau_sd,
                                alpha_sd = 1, eps_sd = 1, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  if (assign_var <= 0) {
    Wm <- rep(min(max(p_bar, 1e-3), 1-1e-3), K)
  } else {
    ab <- rbeta_from_mean_var(p_bar, assign_var)
    Wm <- rbeta(K, shape1 = ab["a"], shape2 = ab["b"])
    Wm <- pmin(pmax(Wm, 1e-3), 1 - 1e-3)
  }
  alpha_m <- rnorm(K, 0, alpha_sd)
  tau_m   <- rnorm(K, tau_mean, tau_sd)
  cl_id <- rep(seq_len(K), each = N_per_cluster)
  W     <- rbinom(K * N_per_cluster, 1, prob = Wm[cl_id])
  eps   <- rnorm(K * N_per_cluster, 0, eps_sd)
  Y0    <- alpha_m[cl_id] + eps
  Y     <- Y0 + tau_m[cl_id] * W
  data.frame(y = Y, w = W, cluster = factor(cl_id), Wm = Wm[cl_id], tau = tau_m[cl_id])
}
sample_clusters <- function(pop, K_pop, q) {
  if (q >= 1) return(list(dat = pop, sampled_clusters = unique(pop$cluster), qhat = 1, K_pop = K_pop))
  K_s <- max(1L, as.integer(round(q * K_pop)))
  s_clusters <- sample(seq_len(K_pop), size = K_s, replace = FALSE)
  dat <- subset(pop, as.integer(cluster) %in% s_clusters)
  list(dat = droplevels(dat), sampled_clusters = s_clusters, qhat = K_s / K_pop, K_pop = K_pop)
}
fit_model_ccv <- function(dat, with_fe = FALSE) {
  if (with_fe) lm(y ~ w + cluster, data = dat) else lm(y ~ w, data = dat)
}
se_ehw <- function(fit) { needs_sandwich(); as.numeric(sqrt(diag(sandwich::vcovHC(fit, type = "HC1"))["w"])) }
se_cluster <- function(fit, dat) { needs_sandwich(); as.numeric(sqrt(diag(sandwich::vcovCL(fit, cluster = dat$cluster, type = "HC1"))["w"])) }
ccv_variance <- function(v_robust, v_cluster, Wm_by_cluster, qhat) {
  Wm <- Wm_by_cluster
  m1 <- mean(Wm * (1 - Wm)); m2 <- mean((Wm^2) * ((1 - Wm)^2))
  lambda_hat <- if (m2 <= 0) 1 else 1 - qhat * (m1^2) / m2
  lambda_hat <- min(max(lambda_hat, 0), 1)
  V_CCV_1 <- lambda_hat * v_cluster + (1 - lambda_hat) * v_robust
  V_CCV   <- qhat * V_CCV_1 + (1 - qhat) * v_cluster
  V_CCV
}
se_ccv <- function(fit, dat, qhat) {
  needs_sandwich()
  Vrob <- as.numeric(sandwich::vcovHC(fit, type = "HC1")["w","w"])
  Vclu <- as.numeric(sandwich::vcovCL(fit, cluster = dat$cluster, type = "HC1")["w","w"])
  tabs <- with(dat, tapply(w, cluster, function(v) mean(v)))
  Vccv <- ccv_variance(Vrob, Vclu, Wm_by_cluster = as.numeric(tabs), qhat = qhat)
  sqrt(Vccv)
}
se_tscb <- function(dat, with_fe = FALSE, B = 200, seed = NULL, progress_hook = NULL) {
  if (!is.null(seed)) set.seed(seed)
  cl <- droplevels(dat$cluster)
  split_idx <- split(seq_len(nrow(dat)), cl)
  Wm <- sapply(split_idx, function(idx) mean(dat$w[idx]))
  Nk <- sapply(split_idx, length)
  Wm <- pmin(pmax(Wm, 1e-3), 1 - 1e-3)
  pools <- lapply(split_idx, function(idx) {
    list(T = dat[idx[dat$w[idx] == 1], , drop = FALSE],
         C = dat[idx[dat$w[idx] == 0], , drop = FALSE])
  })
  boot_b <- numeric(B)
  for (b in seq_len(B)) {
    if (!is.null(progress_hook)) progress_hook(b, B)
    Wb <- sample(Wm, size = length(Wm), replace = TRUE)
    rows <- vector("list", length(Wm))
    for (j in seq_along(Wm)) {
      n_j <- Nk[j]; nT <- max(0L, as.integer(round(n_j * Wb[j]))); nC <- n_j - nT
      T_pool <- if (nrow(pools[[j]]$T) > 0) pools[[j]]$T else pools[[j]]$C
      C_pool <- if (nrow(pools[[j]]$C) > 0) pools[[j]]$C else pools[[j]]$T
      takeT <- T_pool[sample(seq_len(nrow(T_pool)), size = nT, replace = TRUE), , drop = FALSE]
      takeC <- C_pool[sample(seq_len(nrow(C_pool)), size = nC, replace = TRUE), , drop = FALSE]
      rows[[j]] <- rbind(takeT, takeC)
    }
    boot_dat <- do.call(rbind, rows)
    boot_fit <- fit_model_ccv(boot_dat, with_fe = with_fe)
    boot_b[b] <- coef(boot_fit)[["w"]]
  }
  sd(boot_b)
}
one_run_ccv <- function(K_pop, N_per_cluster, q, p_bar, assign_var, tau_mean, tau_sd,
                        alpha_sd, eps_sd, with_fe, B_boot, use_ccv, use_tscb,
                        seed = NULL, pb = NULL) {
  pop <- simulate_population(K_pop, N_per_cluster, p_bar, assign_var, tau_mean, tau_sd,
                             alpha_sd = alpha_sd, eps_sd = eps_sd, seed = seed)
  samp <- sample_clusters(pop, K_pop = K_pop, q = q)
  dat  <- samp$dat; qhat <- samp$qhat
  fit  <- fit_model_ccv(dat, with_fe = with_fe)
  bhat <- coef(fit)[["w"]]
  seR <- tryCatch(se_ehw(fit), error = function(e) NA_real_)
  seC <- tryCatch(se_cluster(fit, dat), error = function(e) NA_real_)
  seV <- if (use_ccv) tryCatch(se_ccv(fit, dat, qhat = qhat), error = function(e) NA_real_) else NA_real_
  seB <- if (use_tscb) {
    tryCatch(se_tscb(dat, with_fe = with_fe, B = B_boot,
                     progress_hook = function(b, B) { if (!is.null(pb)) setTxtProgressBar(pb, b/B) }),
             error = function(e) NA_real_)
  } else NA_real_
  c(bhat = bhat, se_ehw = seR, se_cluster = seC, se_ccv = seV, se_tscb = seB)
}
mod_se_ccv_ui <- function(id) {
  ns <- NS(id)
  fluidPage(
    titlePanel("Clustered SEs vs CCV & TSCB (presets + sweeps)"),
    sidebarLayout(
      sidebarPanel(
        h4("Preset scenarios"),
        selectInput(ns("preset"), NULL,
                    choices = c("— choose —" = "none",
                                "Homogeneous (no het.)"           = "homog",
                                "Heterogeneous assignment only"   = "het_assign",
                                "Heterogeneous effects only"      = "het_tau",
                                "Both strong (q=1)"               = "both_strong",
                                "Small sample fraction (q=0.3)"   = "small_q"),
                    selected = "homog"),
        helpText("Presets update the controls below. You can tweak further after selecting."),
        h4("Population & Sampling"),
        numericInput(ns("K_pop"), "Population clusters (K)", value = 100, min = 10, step = 10),
        numericInput(ns("Npc"), "Units per cluster", value = 100, min = 10, step = 10),
        sliderInput(ns("q"), "Sample fraction of clusters (q)", min = 0.1, max = 1.0, value = 1.0, step = 0.05),
        
        h4("Assignment & Heterogeneity"),
        sliderInput(ns("pbar"), "Mean treatment probability (p̄)", min = 0.05, max = 0.95, value = 0.5, step = 0.05),
        sliderInput(ns("assign_var"), "Between-cluster assignment variance Var(W_m)", min = 0.000, max = 0.08, value = 0.00, step = 0.002),
        numericInput(ns("tau_mean"), "Mean treatment effect (τ̄)", value = 0.5, step = 0.1),
        sliderInput(ns("tau_sd"), "SD of cluster treatment effects SD(τ_m)", min = 0.0, max = 1.0, value = 0.00, step = 0.05),
        
        h4("Nuisance scales"),
        sliderInput(ns("alpha_sd"), "SD of cluster intercepts (α_m)", min = 0.0, max = 2.0, value = 1.0, step = 0.1),
        sliderInput(ns("eps_sd"), "Idiosyncratic noise SD (ε)", min = 0.2, max = 3.0, value = 1.0, step = 0.1),
        
        h4("Estimator & Monte Carlo"),
        selectInput(ns("fe"), "Regression spec", choices = c("OLS (no FE)" = "ols", "Fixed effects (cluster FE)" = "fe")),
        checkboxInput(ns("use_ccv"), "Include CCV (design-based analytic)", TRUE),
        checkboxInput(ns("use_tscb"), "Include TSCB (two-stage cluster bootstrap)", FALSE),
        numericInput(ns("B_boot"), "Bootstrap reps (TSCB)", value = 80, min = 40, step = 20),
        numericInput(ns("MC"), "Monte Carlo reps", value = 120, min = 20, step = 20),
        numericInput(ns("seed"), "Random seed", value = 1234, min = 1, step = 1),
        actionButton(ns("run"), "Run main simulation", class = "btn-primary"),
        hr(),
        
        h4("Design sweeps"),
        selectInput(ns("sweep_dim"), "Sweep dimension", choices = c("q (sample fraction)" = "q", "K (clusters)" = "K")),
        numericInput(ns("sweep_min"), "Grid min", value = 0.1, step = 0.1),
        numericInput(ns("sweep_max"), "Grid max", value = 1.0, step = 0.1),
        numericInput(ns("sweep_step"), "Grid step", value = 0.1, step = 0.05),
        numericInput(ns("MC_sweep"), "MC per grid point", value = 60, min = 20, step = 10),
        numericInput(ns("B_boot_sweep"), "Bootstrap reps (TSCB) per grid", value = 50, min = 30, step = 10),
        selectInput(ns("fe_sweep"), "Spec in sweep", choices = c("OLS (no FE)" = "ols", "Fixed effects (cluster FE)" = "fe"), selected = "ols"),
        checkboxInput(ns("use_ccv_sweep"), "Include CCV in sweep", TRUE),
        checkboxInput(ns("use_tscb_sweep"), "Include TSCB in sweep", FALSE),
        actionButton(ns("run_sweep"), "Run sweep", class = "btn-secondary")
      ),
      mainPanel(
        uiOutput(ns("pkgWarn")),
        h4("Main simulation: Avg SE vs True SD of β̂"),
        plotOutput(ns("seplot"), height = 340),
        h4("Main simulation: Coverage of 95% CIs for β (uses τ̄)"),
        plotOutput(ns("covplot"), height = 260),
        h4("Summary table"),
        DTOutput(ns("tab")),
        hr(),
        h4("Design sweeps"),
        fluidRow(
          column(6, plotOutput(ns("se_sweep_plot"), height = 320)),
          column(6, plotOutput(ns("cov_sweep_plot"), height = 320))
        )
      )
    )
  )
}
mod_se_ccv_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    
    ok <- requireNamespace("sandwich", quietly = TRUE)
    output$pkgWarn <- renderUI({
      if (!ok) {
        div(style="color:#a00;",
            strong("Missing dependency: "),
            "This tab needs the 'sandwich' package. Install with ",
            code("install.packages('sandwich')"), ".")
      } else NULL
    })
    
    observeEvent(input$preset, {
      if (input$preset == "none") return(NULL)
      if (input$preset == "homog") {
        updateSliderInput(session, "assign_var", value = 0.00)
        updateSliderInput(session, "tau_sd",     value = 0.00)
        updateSliderInput(session, "q",          value = 1.00)
      }
      if (input$preset == "het_assign") {
        updateSliderInput(session, "assign_var", value = 0.02)
        updateSliderInput(session, "tau_sd",     value = 0.00)
        updateSliderInput(session, "q",          value = 1.00)
      }
      if (input$preset == "het_tau") {
        updateSliderInput(session, "assign_var", value = 0.00)
        updateSliderInput(session, "tau_sd",     value = 0.40)
        updateSliderInput(session, "q",          value = 1.00)
      }
      if (input$preset == "both_strong") {
        updateSliderInput(session, "assign_var", value = 0.04)
        updateSliderInput(session, "tau_sd",     value = 0.60)
        updateSliderInput(session, "q",          value = 1.00)
      }
      if (input$preset == "small_q") {
        updateSliderInput(session, "assign_var", value = 0.02)
        updateSliderInput(session, "tau_sd",     value = 0.40)
        updateSliderInput(session, "q",          value = 0.30)
      }
    }, ignoreInit = TRUE)
    
    sims <- eventReactive(input$run, {
      validate(need(ok, "Install the 'sandwich' package and click Run again."))
      set.seed(input$seed)
      MC <- input$MC; with_fe <- (input$fe == "fe")
      mat <- matrix(NA_real_, nrow = MC, ncol = 5)
      colnames(mat) <- c("bhat","se_ehw","se_cluster","se_ccv","se_tscb")
      withProgress(message = "Running main simulation...", value = 0, {
        inc <- 1/MC
        for (r in seq_len(MC)) {
          pb <- NULL
          if (isTRUE(input$use_tscb)) pb <- txtProgressBar(min = 0, max = 1, style = 3)
          mat[r, ] <- tryCatch(
            one_run_ccv(
              K_pop = input$K_pop, N_per_cluster = input$Npc, q = input$q,
              p_bar = input$pbar, assign_var = input$assign_var,
              tau_mean = input$tau_mean, tau_sd = input$tau_sd,
              alpha_sd = input$alpha_sd, eps_sd = input$eps_sd,
              with_fe = with_fe, B_boot = input$B_boot,
              use_ccv = input$use_ccv, use_tscb = input$use_tscb,
              seed = sample.int(1e9, 1), pb = pb
            ),
            error = function(e) { rep(NA_real_, 5) }
          )
          if (!is.null(pb)) close(pb)
          incProgress(inc, detail = paste("rep", r, "of", MC))
        }
      })
      as.data.frame(mat)
    })
    
    output$seplot <- renderPlot({
      df <- sims(); req(df)
      true_sd <- stats::sd(df$bhat, na.rm = TRUE)
      bars <- data.frame(
        Method = c("True SD(β̂)", "EHW", "Cluster",
                   if (input$use_ccv) "CCV" else NULL,
                   if (input$use_tscb) "TSCB" else NULL),
        SE     = c(true_sd,
                   mean(df$se_ehw, na.rm = TRUE),
                   mean(df$se_cluster, na.rm = TRUE),
                   if (input$use_ccv) mean(df$se_ccv, na.rm = TRUE) else NULL,
                   if (input$use_tscb) mean(df$se_tscb, na.rm = TRUE) else NULL)
      )
      ggplot(bars, aes(x = Method, y = SE, fill = Method)) +
        geom_col() +
        labs(y = "Scale", x = NULL,
             subtitle = paste0("K_pop=", input$K_pop, ", q=", input$q,
                               ", Var(W_m)=", input$assign_var,
                               ", SD(τ_m)=", input$tau_sd,
                               ifelse(input$fe=="fe", ", spec=FE", ", spec=OLS"))) +
        theme(legend.position = "none")
    })
    
    output$covplot <- renderPlot({
      df <- sims(); req(df)
      beta_true <- input$tau_mean
      cov_fun <- function(b, se) mean(beta_true >= (b - 1.96*se) & beta_true <= (b + 1.96*se), na.rm = TRUE)
      bars <- data.frame(
        Method = c("EHW","Cluster",
                   if (input$use_ccv) "CCV" else NULL,
                   if (input$use_tscb) "TSCB" else NULL),
        Coverage = c(cov_fun(df$bhat, df$se_ehw),
                     cov_fun(df$bhat, df$se_cluster),
                     if (input$use_ccv) cov_fun(df$bhat, df$se_ccv) else NULL,
                     if (input$use_tscb) cov_fun(df$bhat, df$se_tscb) else NULL)
      )
      ggplot(bars, aes(x = Method, y = Coverage, fill = Method)) +
        geom_col() + geom_hline(yintercept = 0.95, linetype = "dashed") +
        ylim(0, 1) + labs(y = "Empirical coverage", x = NULL) +
        theme(legend.position = "none")
    })
    
    output$tab <- renderDT({
      df <- sims(); req(df)
      out <- data.frame(
        Metric = c("True SD(β̂)", "Avg SE (EHW)", "Avg SE (Cluster)",
                   if (input$use_ccv) "Avg SE (CCV)" else NULL,
                   if (input$use_tscb) "Avg SE (TSCB)" else NULL),
        Value  = fmt4(c(stats::sd(df$bhat, na.rm = TRUE),
                        mean(df$se_ehw, na.rm = TRUE),
                        mean(df$se_cluster, na.rm = TRUE),
                        if (input$use_ccv) mean(df$se_ccv, na.rm = TRUE) else NULL,
                        if (input$use_tscb) mean(df$se_tscb, na.rm = TRUE) else NULL))
      )
      datatable(out, rownames = FALSE, options = list(dom = "t", pageLength = 10))
    })
    
    sweep_res <- eventReactive(input$run_sweep, {
      validate(need(ok, "Install the 'sandwich' package and click Run again."))
      with_fe  <- (input$fe_sweep == "fe")
      dim <- input$sweep_dim
      if (dim == "q") {
        grid <- seq(input$sweep_min, input$sweep_max, by = input$sweep_step)
      } else {
        g <- seq(input$sweep_min, input$sweep_max, by = input$sweep_step)
        grid <- as.integer(round(g[g >= 2])); if (length(grid) == 0) grid <- c(20L, 40L, 60L)
      }
      out <- list()
      withProgress(message = paste("Running sweep over", if (dim=="q") "q" else "K"), value = 0, {
        inc <- 1/length(grid)
        for (i in seq_along(grid)) {
          if (dim == "q") { q_i <- max(0.05, min(1, grid[i])); K_i <- input$K_pop }
          else { q_i <- input$q; K_i <- grid[i] }
          mat <- matrix(NA_real_, nrow = input$MC_sweep, ncol = 5)
          colnames(mat) <- c("bhat","se_ehw","se_cluster","se_ccv","se_tscb")
          for (r in seq_len(input$MC_sweep)) {
            pb <- NULL
            if (isTRUE(input$use_tscb_sweep)) pb <- txtProgressBar(min = 0, max = 1, style = 3)
            mat[r, ] <- tryCatch(
              one_run_ccv(
                K_pop = K_i, N_per_cluster = input$Npc, q = q_i,
                p_bar = input$pbar, assign_var = input$assign_var,
                tau_mean = input$tau_mean, tau_sd = input$tau_sd,
                alpha_sd = input$alpha_sd, eps_sd = input$eps_sd,
                with_fe = with_fe, B_boot = input$B_boot_sweep,
                use_ccv = input$use_ccv_sweep, use_tscb = input$use_tscb_sweep,
                seed = sample.int(1e9, 1), pb = pb
              ),
              error = function(e) rep(NA_real_, 5)
            )
            if (!is.null(pb)) close(pb)
          }
          df <- as.data.frame(mat)
          true_sd <- stats::sd(df$bhat, na.rm = TRUE)
          avg <- c(EHW = mean(df$se_ehw, na.rm = TRUE),
                   Cluster = mean(df$se_cluster, na.rm = TRUE),
                   CCV = if (input$use_ccv_sweep) mean(df$se_ccv, na.rm = TRUE) else NA_real_,
                   TSCB = if (input$use_tscb_sweep) mean(df$se_tscb, na.rm = TRUE) else NA_real_)
          cov_fun <- function(b, se) mean(input$tau_mean >= (b - 1.96*se) & input$tau_mean <= (b + 1.96*se), na.rm = TRUE)
          cov <- c(EHW = cov_fun(df$bhat, df$se_ehw),
                   Cluster = cov_fun(df$bhat, df$se_cluster),
                   CCV = if (input$use_ccv_sweep) cov_fun(df$bhat, df$se_ccv) else NA_real_,
                   TSCB = if (input$use_tscb_sweep) cov_fun(df$bhat, df$se_tscb) else NA_real_)
          out[[i]] <- list(grid = grid[i], true_sd = true_sd, avg = avg, cov = cov)
          incProgress(inc, detail = paste(if (dim=="q") "q" else "K", "=", grid[i]))
        }
      })
      gv   <- sapply(out, function(z) z$grid)
      t_sd <- sapply(out, function(z) z$true_sd)
      avgM <- do.call(rbind, lapply(out, function(z) z$avg))
      covM <- do.call(rbind, lapply(out, function(z) z$cov))
      list(grid = gv, true_sd = t_sd, avg = avgM, cov = covM, dim = dim)
    })
    
    output$se_sweep_plot <- renderPlot({
      SR <- sweep_res(); req(SR)
      grid <- SR$grid; dim <- SR$dim
      df <- data.frame(
        x = rep(grid, times = 1 + ncol(SR$avg)),
        Method = factor(c(rep("True SD(β̂)", length(grid)),
                          rep(colnames(SR$avg), each = length(grid))),
                        levels = c("True SD(β̂)","EHW","Cluster","CCV","TSCB")),
        Value = c(SR$true_sd, as.numeric(t(SR$avg)))
      )
      df <- df[!is.na(df$Value), ]
      xlab <- if (dim == "q") "Sample fraction q" else "Number of clusters K"
      ggplot(df, aes(x = x, y = Value, linetype = Method)) +
        geom_line() + geom_point() +
        labs(title = "Avg SE vs True SD across grid", x = xlab, y = "Scale") +
        theme(legend.position = "bottom")
    })
    
    output$cov_sweep_plot <- renderPlot({
      SR <- sweep_res(); req(SR)
      grid <- SR$grid; dim <- SR$dim
      covdf <- data.frame(
        x = rep(grid, times = ncol(SR$cov)),
        Method = factor(rep(colnames(SR$cov), each = length(grid)),
                        levels = c("EHW","Cluster","CCV","TSCB")),
        Coverage = as.numeric(t(SR$cov))
      )
      covdf <- covdf[!is.na(covdf$Coverage), ]
      xlab <- if (dim == "q") "Sample fraction q" else "Number of clusters K"
      ggplot(covdf, aes(x = x, y = Coverage, linetype = Method)) +
        geom_line() + geom_point() +
        geom_hline(yintercept = 0.95, linetype = "dashed") +
        ylim(0, 1) +
        labs(title = "95% CI coverage across grid", x = xlab, y = "Coverage") +
        theme(legend.position = "bottom")
    })
  })
}

# ===========================
# App Shell (5 tabs)
# ===========================
ui <- navbarPage(
  "Econometrics Teaching — Combined Apps",
  tabPanel("ATE vs SDO",         mod_ate_ui("ate")),
  tabPanel("Randomization Inference", mod_ri_ui("ri")),
  tabPanel("Fixed Effects vs OLS",    mod_fe_ui("fe")),
  tabPanel("SE Sim (summary/diag/sweep)", mod_se_summary_ui("sesum")),
  tabPanel("Clustered SEs vs CCV & TSCB", mod_se_ccv_ui("ccv"))
)

server <- function(input, output, session) {
  mod_ate_server("ate")
  mod_ri_server("ri")
  mod_fe_server("fe")
  mod_se_summary_server("sesum")
  mod_se_ccv_server("ccv")
}

shinyApp(ui, server)
