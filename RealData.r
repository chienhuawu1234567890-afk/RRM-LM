############################################################
## Real Data Illustration: Boston Housing Dataset
## RRM + Existing Privacy-Preserving Methods
############################################################

rm(list = ls())
set.seed(12345)

############################################################
## 1. Packages
############################################################

required_pkgs <- c("MASS", "dplyr")

for (pkg in required_pkgs) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

############################################################
## 2. Load Boston Housing Data
############################################################

data("Boston", package = "MASS")

dat <- Boston %>%
  dplyr::select(
    MEDV = medv,
    RM = rm,
    LSTAT = lstat,
    PTRATIO = ptratio,
    NOX = nox,
    DIS = dis
  )

n <- nrow(dat)

formula_main <- MEDV ~ RM + LSTAT + PTRATIO + NOX + DIS

############################################################
## 3. Original Model
############################################################

fit_original <- lm(formula_main, data = dat)

extract_results <- function(fit, data, method_name,
                            Delta2 = NA_real_,
                            sigma_DP = NA_real_,
                            epsilon = NA_real_,
                            delta = NA_real_) {
  
  sm <- summary(fit)
  
  beta_RM <- sm$coefficients["RM", "Estimate"]
  se_RM   <- sm$coefficients["RM", "Std. Error"]
  p_RM    <- sm$coefficients["RM", "Pr(>|t|)"]
  
  pred <- predict(fit, newdata = data)
  rmse <- sqrt(mean((data$MEDV - pred)^2))
  
  data.frame(
    Method = method_name,
    Beta_RM = beta_RM,
    SE_RM = se_RM,
    p_value_RM = p_RM,
    Prediction_RMSE = rmse,
    Delta2 = Delta2,
    sigma_DP = sigma_DP,
    epsilon = epsilon,
    delta = delta,
    stringsAsFactors = FALSE
  )
}

############################################################
## 4. Proposed Moment-Invariant RRM
############################################################

moment_invariant_rrm_lm <- function(data, formula, gamma = 0.50) {
  
  fit <- lm(formula, data = data)
  
  y <- model.response(model.frame(fit))
  mu_hat <- fitted(fit)
  sigma_hat <- summary(fit)$sigma
  
  eps_sd <- sqrt(1 - gamma^2) * sigma_hat
  
  y_star <- mu_hat +
    gamma * (y - mu_hat) +
    rnorm(length(y), mean = 0, sd = eps_sd)
  
  data_star <- data
  response_name <- all.vars(formula)[1]
  data_star[[response_name]] <- y_star
  
  data_star
}

gamma <- 0.50

dat_rrm <- moment_invariant_rrm_lm(
  data = dat,
  formula = formula_main,
  gamma = gamma
)

fit_rrm <- lm(formula_main, data = dat_rrm)

############################################################
## 5. Existing Method 1: Random Swap
############################################################

random_swap <- function(data, swap_rate = 0.30) {
  
  data_star <- data
  n <- nrow(data)
  n_swap <- floor(swap_rate * n)
  
  swap_ids <- sample(seq_len(n), n_swap)
  data_star$MEDV[swap_ids] <- sample(data_star$MEDV[swap_ids])
  
  data_star
}

dat_swap <- random_swap(dat, swap_rate = 0.30)
fit_swap <- lm(formula_main, data = dat_swap)

############################################################
## 6. Existing Method 2: PRAM-like Perturbation
############################################################

pram_response <- function(data, perturb_rate = 0.30) {
  
  data_star <- data
  n <- nrow(data)
  
  ids <- sample(seq_len(n), floor(perturb_rate * n))
  
  y <- data$MEDV
  y_sd <- sd(y)
  
  data_star$MEDV[ids] <- data_star$MEDV[ids] +
    rnorm(length(ids), mean = 0, sd = 0.50 * y_sd)
  
  data_star
}

dat_pram <- pram_response(dat, perturb_rate = 0.30)
fit_pram <- lm(formula_main, data = dat_pram)

############################################################
## 7. Existing Method 3: Randomized Response
############################################################

randomized_response <- function(data, response_rate = 0.70) {
  
  data_star <- data
  n <- nrow(data)
  
  replace_ids <- which(runif(n) > response_rate)
  
  data_star$MEDV[replace_ids] <- sample(
    data$MEDV,
    length(replace_ids),
    replace = TRUE
  )
  
  data_star
}

dat_rr <- randomized_response(dat, response_rate = 0.70)
fit_rr <- lm(formula_main, data = dat_rr)

############################################################
## 8. Existing Method 4: Synthetic Data
############################################################

synthetic_data <- function(data, formula) {
  
  fit <- lm(formula, data = data)
  mu_hat <- fitted(fit)
  sigma_hat <- summary(fit)$sigma
  
  data_star <- data
  data_star$MEDV <- rnorm(
    nrow(data),
    mean = mu_hat,
    sd = sigma_hat
  )
  
  data_star
}

dat_syn <- synthetic_data(dat, formula_main)
fit_syn <- lm(formula_main, data = dat_syn)

############################################################
## 9. Existing Method 5: Standard Gaussian DP
############################################################

standard_gaussian_dp <- function(y,
                                 B = 5,
                                 epsilon = 1,
                                 delta = 1e-5) {
  
  y_scaled <- as.numeric(scale(y))
  y_clipped <- pmax(pmin(y_scaled, B), -B)
  
  Delta2 <- 2 * B
  
  sigma_DP <- Delta2 *
    sqrt(2 * log(1.25 / delta)) / epsilon
  
  y_dp_scaled <- y_clipped +
    rnorm(length(y), mean = 0, sd = sigma_DP)
  
  y_dp <- as.numeric(y_dp_scaled * sd(y) + mean(y))
  
  list(
    y_dp = y_dp,
    Delta2 = Delta2,
    sigma_DP = sigma_DP,
    epsilon = epsilon,
    delta = delta
  )
}

############################################################
## 10. Existing Method 6: DP-RRM
############################################################

privacy_calibrated_dp_rrm <- function(y,
                                      gamma = 0.50,
                                      B = 5,
                                      epsilon = 1,
                                      delta = 1e-5) {
  
  y_scaled <- as.numeric(scale(y))
  y_clipped <- pmax(pmin(y_scaled, B), -B)
  
  Delta2 <- 2 * B * gamma
  
  sigma_DP <- Delta2 *
    sqrt(2 * log(1.25 / delta)) / epsilon
  
  y_dp_scaled <- gamma * y_clipped +
    rnorm(length(y), mean = 0, sd = sigma_DP)
  
  y_dp <- as.numeric(y_dp_scaled * sd(y) + mean(y))
  
  list(
    y_dp = y_dp,
    Delta2 = Delta2,
    sigma_DP = sigma_DP,
    epsilon = epsilon,
    delta = delta
  )
}

epsilon <- 1
delta <- 1e-5
B_clip <- 5

dp_obj <- standard_gaussian_dp(
  y = dat$MEDV,
  B = B_clip,
  epsilon = epsilon,
  delta = delta
)

dat_dp <- dat
dat_dp$MEDV <- dp_obj$y_dp
fit_dp <- lm(formula_main, data = dat_dp)

dp_rrm_obj <- privacy_calibrated_dp_rrm(
  y = dat$MEDV,
  gamma = gamma,
  B = B_clip,
  epsilon = epsilon,
  delta = delta
)

dat_dp_rrm <- dat
dat_dp_rrm$MEDV <- dp_rrm_obj$y_dp
fit_dp_rrm <- lm(formula_main, data = dat_dp_rrm)

############################################################
## 11. Geometric Diagnostics
############################################################

mahal_diag <- function(data_original, data_release) {
  
  X0 <- as.matrix(data_original)
  X1 <- as.matrix(data_release)
  
  D0 <- mahalanobis(X0, colMeans(X0), cov(X0))
  D1 <- mahalanobis(X1, colMeans(X1), cov(X1))
  
  cutoff <- qchisq(0.975, df = ncol(X0))
  
  O0 <- as.numeric(D0 > cutoff)
  O1 <- as.numeric(D1 > cutoff)
  
  ROP <- ifelse(
    sum(O0) == 0,
    NA_real_,
    100 * sum(O0 * O1) / sum(O0)
  )
  
  OOR <- ifelse(
    sum((O0 + O1) > 0) == 0,
    NA_real_,
    sum(O0 * O1) / sum((O0 + O1) > 0)
  )
  
  rho_M <- suppressWarnings(
    cor(rank(D0), rank(D1), method = "spearman")
  )
  
  Delta_D <- mean(abs(D0 - D1))
  
  list(
    ROP = ROP,
    OOR = OOR,
    rho_Mahalanobis = rho_M,
    Delta_D = Delta_D
  )
}

cook_diag <- function(fit_original, fit_release) {
  
  CD0 <- cooks.distance(fit_original)
  CD1 <- cooks.distance(fit_release)
  
  rho_Cook <- suppressWarnings(
    cor(rank(CD0), rank(CD1), method = "spearman")
  )
  
  Delta_CD <- mean(abs(CD0 - CD1))
  
  list(
    rho_Cook = rho_Cook,
    Delta_CD = Delta_CD
  )
}

############################################################
## 12. Disclosure-Risk Diagnostics
############################################################

nearest_neighbor_risk <- function(data_original,
                                  data_release,
                                  k_values = c(1, 5, 10)) {
  
  X0 <- as.matrix(scale(data_original))
  X1 <- as.matrix(scale(data_release))
  
  n <- nrow(X0)
  ranks <- numeric(n)
  
  for (i in seq_len(n)) {
    dists <- rowSums(
      (X0 - matrix(
        X1[i, ],
        nrow = n,
        ncol = ncol(X0),
        byrow = TRUE
      ))^2
    )
    
    ranks[i] <- rank(dists, ties.method = "first")[i]
  }
  
  out <- list()
  out$R_NN <- mean(ranks == 1)
  
  for (k in k_values) {
    out[[paste0("R_", k)]] <- mean(ranks <= k)
  }
  
  out
}

reconstruction_error <- function(data_original, data_release) {
  
  X0 <- as.matrix(scale(data_original))
  X1 <- as.matrix(scale(data_release))
  
  mean(sqrt(rowSums((X0 - X1)^2)))
}

release_variability_lm <- function(data,
                                   formula,
                                   gamma = 0.50,
                                   B_release = 200) {
  
  releases <- replicate(
    B_release,
    {
      dat_star <- moment_invariant_rrm_lm(
        data = data,
        formula = formula,
        gamma = gamma
      )
      dat_star$MEDV
    }
  )
  
  V_i <- apply(releases, 1, var)
  
  mean(V_i)
}

############################################################
## 13. Helper Function for Diagnostics
############################################################

get_all_diagnostics <- function(data_release, fit_release) {
  
  geo <- mahal_diag(dat, data_release)
  cook <- cook_diag(fit_original, fit_release)
  disc <- nearest_neighbor_risk(dat, data_release)
  RE <- reconstruction_error(dat, data_release)
  
  list(
    ROP = geo$ROP,
    OOR = geo$OOR,
    rho_Cook = cook$rho_Cook,
    rho_Mahalanobis = geo$rho_Mahalanobis,
    Delta_D = geo$Delta_D,
    Delta_CD = cook$Delta_CD,
    R_NN = disc$R_NN,
    R_1 = disc$R_1,
    R_5 = disc$R_5,
    R_10 = disc$R_10,
    Reconstruction_Error = RE
  )
}

diag_rrm    <- get_all_diagnostics(dat_rrm, fit_rrm)
diag_swap   <- get_all_diagnostics(dat_swap, fit_swap)
diag_pram   <- get_all_diagnostics(dat_pram, fit_pram)
diag_rr     <- get_all_diagnostics(dat_rr, fit_rr)
diag_syn    <- get_all_diagnostics(dat_syn, fit_syn)
diag_dp     <- get_all_diagnostics(dat_dp, fit_dp)
diag_dp_rrm <- get_all_diagnostics(dat_dp_rrm, fit_dp_rrm)

V_rrm <- release_variability_lm(
  data = dat,
  formula = formula_main,
  gamma = gamma,
  B_release = 200
)

############################################################
## 14. Combine Estimation Results
############################################################

results <- rbind(
  extract_results(fit_original, dat, "Original"),
  extract_results(fit_rrm, dat_rrm, "Proposed RRM"),
  extract_results(fit_swap, dat_swap, "Random Swap"),
  extract_results(fit_pram, dat_pram, "PRAM"),
  extract_results(fit_rr, dat_rr, "Randomized Response"),
  extract_results(fit_syn, dat_syn, "Synthetic Data"),
  extract_results(
    fit_dp,
    dat_dp,
    "Differential Privacy",
    Delta2 = dp_obj$Delta2,
    sigma_DP = dp_obj$sigma_DP,
    epsilon = dp_obj$epsilon,
    delta = dp_obj$delta
  ),
  extract_results(
    fit_dp_rrm,
    dat_dp_rrm,
    "DP-RRM",
    Delta2 = dp_rrm_obj$Delta2,
    sigma_DP = dp_rrm_obj$sigma_DP,
    epsilon = dp_rrm_obj$epsilon,
    delta = dp_rrm_obj$delta
  )
)

############################################################
## 15. Add Diagnostics
############################################################

diag_list <- list(
  NULL,
  diag_rrm,
  diag_swap,
  diag_pram,
  diag_rr,
  diag_syn,
  diag_dp,
  diag_dp_rrm
)

add_diag_col <- function(name) {
  sapply(diag_list, function(x) {
    if (is.null(x)) NA_real_ else x[[name]]
  })
}

results$ROP <- add_diag_col("ROP")
results$OOR <- add_diag_col("OOR")
results$rho_Cook <- add_diag_col("rho_Cook")
results$rho_Mahalanobis <- add_diag_col("rho_Mahalanobis")
results$Delta_D <- add_diag_col("Delta_D")
results$Delta_CD <- add_diag_col("Delta_CD")

results$R_NN <- add_diag_col("R_NN")
results$R_1 <- add_diag_col("R_1")
results$R_5 <- add_diag_col("R_5")
results$R_10 <- add_diag_col("R_10")
results$Reconstruction_Error <- add_diag_col("Reconstruction_Error")

results$Repeated_Release_Var <- c(
  NA_real_,
  V_rrm,
  NA_real_,
  NA_real_,
  NA_real_,
  NA_real_,
  NA_real_,
  NA_real_
)

############################################################
## 16. Print Full Results
############################################################

results_print <- results

num_cols <- sapply(results_print, is.numeric)

results_print[num_cols] <- lapply(
  results_print[num_cols],
  round,
  digits = 4
)

print(results_print)

############################################################
## 17. Separate Tables
############################################################

estimation_table <- results %>%
  dplyr::select(
    Method,
    Beta_RM,
    SE_RM,
    p_value_RM,
    Prediction_RMSE
  )

geometric_table <- results %>%
  dplyr::select(
    Method,
    ROP,
    OOR,
    rho_Cook,
    rho_Mahalanobis,
    Delta_D,
    Delta_CD
  )

disclosure_table <- results %>%
  dplyr::select(
    Method,
    R_NN,
    R_1,
    R_5,
    R_10,
    Reconstruction_Error,
    Repeated_Release_Var
  )

privacy_table <- results %>%
  dplyr::select(
    Method,
    Delta2,
    sigma_DP,
    epsilon,
    delta
  )

comparison_table <- data.frame(
  Method = results$Method,
  Utility_Preservation = c(
    "Excellent",
    "Moderate--High",
    "Moderate",
    "Moderate",
    "Low--Moderate",
    "High",
    "Very Low",
    "Low"
  ),
  Disclosure_Protection = c(
    "None",
    "Strong",
    "Moderate",
    "Moderate",
    "Moderate",
    "Moderate",
    "Very Strong",
    "Strong"
  ),
  Formal_DP = c(
    "No",
    "No",
    "No",
    "No",
    "No",
    "No",
    "Yes",
    "Yes"
  )
)

############################################################
## 18. Print Separate Tables
############################################################

print(estimation_table %>% mutate(across(where(is.numeric), ~ round(.x, 4))))
print(geometric_table %>% mutate(across(where(is.numeric), ~ round(.x, 4))))
print(disclosure_table %>% mutate(across(where(is.numeric), ~ round(.x, 4))))
print(privacy_table %>% mutate(across(where(is.numeric), ~ round(.x, 4))))
print(comparison_table)

############################################################
## 19. Save Results
############################################################

write.csv(results, "real_data_rrm_all_methods_results.csv", row.names = FALSE)

write.csv(estimation_table, "real_data_estimation_results.csv", row.names = FALSE)

write.csv(geometric_table, "real_data_geometric_diagnostics.csv", row.names = FALSE)

write.csv(disclosure_table, "real_data_disclosure_diagnostics.csv", row.names = FALSE)

write.csv(privacy_table, "real_data_privacy_calibration.csv", row.names = FALSE)

write.csv(comparison_table, "real_data_method_comparison.csv", row.names = FALSE)

############################################################
## End
############################################################
