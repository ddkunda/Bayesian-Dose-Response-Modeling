# app.R

library(shiny)
library(rstan)
library(ggplot2)
library(dplyr)
library(parallel)

# For faster Stan compilation (caching) and parallel sampling:
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# ---------------------------------------------------------------------------
# STAN MODEL:
#   We model the mean response at dose d:
#     E(d) = E0 + (Emax * d^h) / (EC50^h + d^h)
#
#   We assume a Normal likelihood for observed data with known or estimated SD.
#
#   NOTE: We use the new array syntax required by Stan v2.32+.
# ---------------------------------------------------------------------------
stan_code <- "
data {
  int<lower=1> N;                // number of dose levels
  array[N] real dose;            // dose values
  array[N] real response;        // observed mean percent reduction
  real<lower=0> sigma;           // known SD (if not known, you can turn this into a parameter)

  // Prior parameters (means and SDs)
  real prior_E0_mu;
  real<lower=0> prior_E0_sd;

  real prior_Emax_mu;
  real<lower=0> prior_Emax_sd;

  real prior_EC50_mu;
  real<lower=0> prior_EC50_sd;

  real prior_h_mu;
  real<lower=0> prior_h_sd;
}
parameters {
  real E0;
  real Emax;
  real<lower=0> EC50;
  real<lower=0> h;
}
model {
  // Priors
  E0   ~ normal(prior_E0_mu,   prior_E0_sd);
  Emax ~ normal(prior_Emax_mu, prior_Emax_sd);
  EC50 ~ normal(prior_EC50_mu, prior_EC50_sd);
  h    ~ normal(prior_h_mu,    prior_h_sd);

  // Likelihood
  for (i in 1:N) {
    real mu_i = E0 + (Emax * pow(dose[i], h)) / (pow(EC50, h) + pow(dose[i], h));
    response[i] ~ normal(mu_i, sigma);
  }
}
"

# Compile the model once
stan_model_obj <- stan_model(model_code = stan_code)


# ---------------------------------------------------------------------------
# SHINY UI
# ---------------------------------------------------------------------------
ui <- fluidPage(
  titlePanel("Bayesian Dose-Response Modeling"),
  
  sidebarLayout(
    sidebarPanel(
      h4("1. Input Data"),
      helpText("Enter dose levels and observed percent reduction means (comma-separated)."),
      textAreaInput("dose_levels", "Dose Levels (e.g., mg/kg)", "0, 5, 10, 20, 40", rows = 2),
      textAreaInput("observed_reduction", "Observed Percent Reduction (%)", "5, 15, 30, 45, 50", rows = 2),
      numericInput("std_dev", "Known SD of Observations", value = 5, step = 1),
      
      h4("2. Priors for Emax Model"),
      fluidRow(
        column(6, numericInput("E0_mu",  "E0 prior mean",  value = 5, step = 1)),
        column(6, numericInput("E0_sd",  "E0 prior SD",    value = 5, step = 1))
      ),
      fluidRow(
        column(6, numericInput("Emax_mu", "Emax prior mean", value = 50, step = 5)),
        column(6, numericInput("Emax_sd", "Emax prior SD",   value = 10, step = 1))
      ),
      fluidRow(
        column(6, numericInput("EC50_mu", "EC50 prior mean", value = 15, step = 1)),
        column(6, numericInput("EC50_sd", "EC50 prior SD",   value = 5, step = 1))
      ),
      fluidRow(
        column(6, numericInput("h_mu",    "h prior mean",    value = 1, step = 0.1)),
        column(6, numericInput("h_sd",    "h prior SD",      value = 0.5, step = 0.1))
      ),
      
      h4("3. Target Effect Probability"),
      numericInput("target_effect", "Target Reduction (%)", value = 40, step = 5),
      helpText("We'll compute P(E(d) > target). For each dose, we can see which dose crosses the threshold."),
      
      actionButton("runAnalysis", "Run Bayesian Analysis")
    ),
    
    mainPanel(
      h4(HTML("<u>Posterior Summaries</u>")),
      verbatimTextOutput("posteriorSummary"),
      
      br(),
      
      h4(HTML("<u>Dose-Response Plot</u>")),
      plotOutput("doseResponsePlot", height = "400px"),
      
      br(),
      
      h4(HTML("<u>Target Effect Probability per Dose</u>")),
      tableOutput("targetEffectTable"),
      
      h5(
        span("All questions can be sent to Christian Dide-Agossou, PhD:", style = "color:black; display: inline-block;"),
        span(htmlOutput("uicmt14"), style = "display: inline-block;")
      )
    )
  )
)

# ---------------------------------------------------------------------------
# SHINY SERVER
# ---------------------------------------------------------------------------
server <- function(input, output, session) {
  
  # Run the Bayesian dose-response model upon user click
  runModel <- eventReactive(input$runAnalysis, {
    # 1. Parse user data
    dose_vec <- as.numeric(strsplit(input$dose_levels, "[,\\s]+")[[1]])
    resp_vec <- as.numeric(strsplit(input$observed_reduction, "[,\\s]+")[[1]])
    
    # Basic checks
    validate(
      need(length(dose_vec) == length(resp_vec), "Dose and response vectors must have the same length!")
    )
    
    # 2. Build Stan data
    stan_data_list2 <- list(
      N = length(dose_vec),
      dose = dose_vec,
      response = resp_vec,
      sigma = input$std_dev,    # Known SD
      prior_E0_mu    = input$E0_mu,
      prior_E0_sd    = input$E0_sd,
      prior_Emax_mu  = input$Emax_mu,
      prior_Emax_sd  = input$Emax_sd,
      prior_EC50_mu  = input$EC50_mu,
      prior_EC50_sd  = input$EC50_sd,
      prior_h_mu     = input$h_mu,
      prior_h_sd     = input$h_sd
    )
    
    # 3. Fit model
    fit <- sampling(
      stan_model_obj,
      data = stan_data_list2,
      iter = 3000,
      warmup = 1000,
      chains = 4,
      seed = 42
    )
    
    return(list(fit = fit, dose_vec = dose_vec, resp_vec = resp_vec))
  })
  
  # Show posterior summary
  output$posteriorSummary <- renderPrint({
    req(runModel())
    fit <- runModel()$fit
    print(summary(fit, probs = c(0.025, 0.5, 0.975))$summary)
  })
  
  # Plot the posterior dose-response curve
  # We'll sample from the posterior and compute the mean E(d) for a sequence of doses
  output$doseResponsePlot <- renderPlot({
    req(runModel())
    fit      <- runModel()$fit
    dose_vec <- runModel()$dose_vec
    resp_vec <- runModel()$resp_vec
    
    posterior_draws <- rstan::extract(fit)
    
    # 1. Create a grid of doses for plotting
    dose_grid <- seq(0, max(dose_vec)*1.2, length.out = 100)
    
    # 2. For each dose in dose_grid, compute the entire posterior distribution of E(d)
    #    We'll store draws in a matrix: rows = posterior draws, columns = doses
    ndraws <- length(posterior_draws$E0)  # total MCMC draws
    E_draws_mat <- sapply(dose_grid, function(d) {
      E0_draws   <- posterior_draws$E0
      Emax_draws <- posterior_draws$Emax
      EC50_draws <- posterior_draws$EC50
      h_draws    <- posterior_draws$h
      
      # Vector of E(d) for each posterior draw
      E0_draws + (Emax_draws * d^h_draws) / (EC50_draws^h_draws + d^h_draws)
    })
    # E_draws_mat: dimensions ndraws x length(dose_grid)
    
    # 3. Compute mean, lower CI, and upper CI for each dose
    E_mean  <- apply(E_draws_mat, 2, mean)
    E_lower <- apply(E_draws_mat, 2, quantile, probs = 0.025)
    E_upper <- apply(E_draws_mat, 2, quantile, probs = 0.975)
    
    df_plot <- data.frame(
      dose     = dose_grid,
      mean_est = E_mean,
      lower    = E_lower,
      upper    = E_upper
    )
    
    # 4. Plot
    ggplot() +
      # Ribbon for 95% credible interval
      geom_ribbon(
        data = df_plot,
        aes(x = dose, ymin = lower, ymax = upper),
        fill = "blue", alpha = 0.2
      ) +
      # Mean posterior dose-response curve
      geom_line(
        data = df_plot,
        aes(x = dose, y = mean_est),
        color = "blue", linewidth = 1
      ) +
      # Observed points
      geom_point(
        data = data.frame(dose = dose_vec, resp = resp_vec),
        aes(x = dose, y = resp),
        color = "red", size = 4
      ) +
      theme_minimal() +
      labs(
        x = "Dose (e.g., mg/kg)",
        y = "Percent Reduction (%)",
        title = "Posterior Mean Dose-Response Curve (with 95% Credible Interval)"
      )+
      theme_minimal(base_size = 14) +     # Increase the overall base font size
      theme(
        plot.title   = element_text(size = 16, face = "bold"), # Title size
        axis.title.x = element_text(size = 16, face = "bold"), # X-axis label size
        axis.title.y = element_text(size = 16, face = "bold"), # Y-axis label size
        axis.text.x  = element_text(size = 14),                # X tick label size
        axis.text.y  = element_text(size = 14)                 # Y tick label size
      )
  })
  
  
  # Probability that E(d) > target, for each dose in the original data
  output$targetEffectTable <- renderTable({
    req(runModel())
    fit      <- runModel()$fit
    dose_vec <- runModel()$dose_vec
    posterior_draws <- rstan::extract(fit)
    
    target_val <- input$target_effect
    
    # For each dose in dose_vec, compute the posterior distribution of E(d),
    # then compute the fraction of draws that exceed target_val.
    results <- lapply(seq_along(dose_vec), function(i){
      d <- dose_vec[i]
      E0_draws   <- posterior_draws$E0
      Emax_draws <- posterior_draws$Emax
      EC50_draws <- posterior_draws$EC50
      h_draws    <- posterior_draws$h
      
      Ed_draws <- E0_draws + (Emax_draws * d^h_draws) / (EC50_draws^h_draws + d^h_draws)
      prob     <- mean(Ed_draws > target_val)
      
      data.frame(
        Dose = d,
        "Probability E(d) > Target" = sprintf("%.1f%%", 100 * prob),
        check.names = FALSE
      )
    })
    
    do.call(rbind, results)
  }, digits = 0)
  
  output$uicmt14 <- renderUI({
    email <- "christian.dideagossou@gmail.com"
    link <- paste0("mailto:", email)
    tags$a(href = link, email)  # Use tags$a to create the hyperlink
  })
}

shinyApp(ui = ui, server = server)
