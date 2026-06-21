############################################################
## Randomized Refraction Model Simulation Codes
## Studies I--V: Inference, Robustness, Diagnostics, Privacy
############################################################

rm(list = ls())
set.seed(12345)

############################################################
## 0. Helper Functions
############################################################

ols_summary <- function(x, y, beta_true = 1, alpha = 0.05) {
  fit <- lm(y ~ x)
  est <- coef(summary(fit))["x", "Estimate"]
  se  <- coef(summary(fit))["x", "Std. Error"]
  pval <- coef(summary(fit))["x", "Pr(>|t|)"]

  ci_l <- est + qt(alpha / 2, df = length(y) - 2) * se
  ci_u <- est + qt(1 - alpha / 2, df = length(y) - 2) * se

  c(
    beta_hat = est,
    se = se,
    cover = as.numeric(ci_l <= beta_true & beta_true <= ci_u),
    reject = as.numeric(pval < alpha)
  )
}

moment_rrm <- function(y, mu, gamma, sigma = 1) {
  n <- length(y)
  e_star <- rnorm(n, mean = 0, sd = sqrt(1 - gamma^2) * sigma)
  y_star <- mu + gamma * (y - mu) + e_star
  y_star
}

cook_diag <- function(x, y) {
  fit <- lm(y ~ x)
  cooks.distance(fit)
}

mahal_diag <- function(x, y) {
  z <- cbind(x, y)
  center <- colMeans(z)
  S <- cov(z)
  mahalanobis(z, center = center, cov = S)
}

safe_cor <- function(a, b) {
  if (sd(a) == 0 || sd(b) == 0) return(NA_real_)
  cor(a, b, method = "spearman")
}

gaussian_dp_noise_sd <- function(B, gamma, epsilon, delta) {
  delta_sens <- 2 * B * gamma
  delta_sens / epsilon * sqrt(2 * log(1.25 / delta))
}

clip_vec <- function(y, B) {
  pmax(pmin(y, B), -B)
}

############################################################
## 1. Main Simulation Settings
############################################################

R <- 2000
n_vec <- c(100, 250, 500)
gamma_vec <- c(0.2, 0.5, 0.8, 1.0)

beta0 <- 0
beta1 <- 1
sigma <- 1
alpha <- 0.05

############################################################
## Study I: Inferential Preservation
############################################################

study1_results <- data.frame()

for (n in n_vec) {
  for (gamma in gamma_vec) {

    store_original <- matrix(NA, R, 4)
    store_rrm <- matrix(NA, R, 4)

    for (r in 1:R) {
      x <- rnorm(n)
      mu <- beta0 + beta1 * x
      y <- mu + rnorm(n, 0, sigma)

      y_star <- moment_rrm(y, mu, gamma, sigma)

      store_original[r, ] <- ols_summary(x, y, beta_true = beta1)
      store_rrm[r, ] <- ols_summary(x, y_star, beta_true = beta1)
    }

    colnames(store_original) <- c("beta_hat", "se", "cover", "reject")
    colnames(store_rrm) <- c("beta_hat", "se", "cover", "reject")

    tmp <- rbind(
      data.frame(
        Study = "Inferential preservation",
        Method = "Original",
        n = n,
        gamma = gamma,
        Bias = mean(store_original[, "beta_hat"]) - beta1,
        SD = sd(store_original[, "beta_hat"]),
        MSE = mean((store_original[, "beta_hat"] - beta1)^2),
        Mean_SE = mean(store_original[, "se"]),
        Coverage = mean(store_original[, "cover"]),
        Power = mean(store_original[, "reject"])
      ),
      data.frame(
        Study = "Inferential preservation",
        Method = "Moment-invariant RRM",
        n = n,
        gamma = gamma,
        Bias = mean(store_rrm[, "beta_hat"]) - beta1,
        SD = sd(store_rrm[, "beta_hat"]),
        MSE = mean((store_rrm[, "beta_hat"] - beta1)^2),
        Mean_SE = mean(store_rrm[, "se"]),
        Coverage = mean(store_rrm[, "cover"]),
        Power = mean(store_rrm[, "reject"])
      )
    )

    study1_results <- rbind(study1_results, tmp)
  }
}

############################################################
## Study II: Robustness and Influence Attenuation
############################################################

study2_results <- data.frame()

for (n in n_vec) {
  for (gamma in gamma_vec) {

    max_cd_original <- max_cd_rrm <- numeric(R)
    mean_cd_original <- mean_cd_rrm <- numeric(R)
    rho_cook <- numeric(R)
    delta_cd <- numeric(R)

    for (r in 1:R) {
      x <- rnorm(n)
      mu <- beta0 + beta1 * x
      y <- mu + rnorm(n, 0, sigma)

      y_star <- moment_rrm(y, mu, gamma, sigma)

      cd <- cook_diag(x, y)
      cd_star <- cook_diag(x, y_star)

      max_cd_original[r] <- max(cd)
      max_cd_rrm[r] <- max(cd_star)

      mean_cd_original[r] <- mean(cd)
      mean_cd_rrm[r] <- mean(cd_star)

      rho_cook[r] <- safe_cor(rank(cd), rank(cd_star))
      delta_cd[r] <- mean(abs(cd - cd_star))
    }

    study2_results <- rbind(
      study2_results,
      data.frame(
        Study = "Robustness and influence",
        n = n,
        gamma = gamma,
        Max_CD_Original = mean(max_cd_original),
        Max_CD_RRM = mean(max_cd_rrm),
        Mean_CD_Original = mean(mean_cd_original),
        Mean_CD_RRM = mean(mean_cd_rrm),
        Delta_CD = mean(delta_cd),
        Rho_Cook = mean(rho_cook, na.rm = TRUE)
      )
    )
  }
}

############################################################
## Study III: Geometric Diagnostics
############################################################

study3_results <- data.frame()

for (n in n_vec) {
  for (gamma in gamma_vec) {

    n_pre <- n_post <- numeric(R)
    ROP <- OOR <- rho_mahal <- delta_D <- numeric(R)

    threshold <- qchisq(0.975, df = 2)

    for (r in 1:R) {
      x <- rnorm(n)
      mu <- beta0 + beta1 * x
      y <- mu + rnorm(n, 0, sigma)

      y_star <- moment_rrm(y, mu, gamma, sigma)

      D <- mahal_diag(x, y)
      D_star <- mahal_diag(x, y_star)

      O <- D > threshold
      O_star <- D_star > threshold

      n_pre[r] <- sum(O)
      n_post[r] <- sum(O_star)

      ROP[r] <- ifelse(sum(O) > 0, 100 * sum(O & O_star) / sum(O), NA)
      OOR[r] <- ifelse(sum(O | O_star) > 0, sum(O & O_star) / sum(O | O_star), NA)

      rho_mahal[r] <- safe_cor(rank(D), rank(D_star))
      delta_D[r] <- mean(abs(D - D_star))
    }

    study3_results <- rbind(
      study3_results,
      data.frame(
        Study = "Geometric diagnostics",
        n = n,
        gamma = gamma,
        N_Pre = mean(n_pre),
        N_Post = mean(n_post),
        ROP = mean(ROP, na.rm = TRUE),
        OOR = mean(OOR, na.rm = TRUE),
        Rho_Mahalanobis = mean(rho_mahal, na.rm = TRUE),
        Delta_D = mean(delta_D)
      )
    )
  }
}

############################################################
## Study IV: Privacy Calibration and Benchmark Comparisons
############################################################

epsilon_vec <- c(0.5, 1, 2)
delta <- 1e-5
B <- 5

study4_results <- data.frame()

for (n in n_vec) {
  for (gamma in gamma_vec) {
    for (epsilon in epsilon_vec) {

      store_original <- matrix(NA, R, 4)
      store_gaussian_dp <- matrix(NA, R, 4)
      store_rrm_dp <- matrix(NA, R, 4)

      sd_gaussian_dp <- gaussian_dp_noise_sd(B, gamma = 1, epsilon, delta)
      sd_rrm_dp <- gaussian_dp_noise_sd(B, gamma = gamma, epsilon, delta)

      for (r in 1:R) {
        x <- rnorm(n)
        mu <- beta0 + beta1 * x
        y <- mu + rnorm(n, 0, sigma)

        y_clip <- clip_vec(y, B)

        y_gaussian_dp <- y_clip + rnorm(n, 0, sd_gaussian_dp)

        y_rrm_dp <- gamma * y_clip + rnorm(n, 0, sd_rrm_dp)

        store_original[r, ] <- ols_summary(x, y, beta_true = beta1)
        store_gaussian_dp[r, ] <- ols_summary(x, y_gaussian_dp, beta_true = beta1)
        store_rrm_dp[r, ] <- ols_summary(x, y_rrm_dp, beta_true = beta1)
      }

      colnames(store_original) <- c("beta_hat", "se", "cover", "reject")
      colnames(store_gaussian_dp) <- c("beta_hat", "se", "cover", "reject")
      colnames(store_rrm_dp) <- c("beta_hat", "se", "cover", "reject")

      tmp <- rbind(
        data.frame(
          Study = "Privacy calibration",
          Method = "Original",
          n = n,
          gamma = gamma,
          epsilon = epsilon,
          delta = delta,
          Noise_SD = 0,
          Bias = mean(store_original[, "beta_hat"]) - beta1,
          SD = sd(store_original[, "beta_hat"]),
          MSE = mean((store_original[, "beta_hat"] - beta1)^2),
          Coverage = mean(store_original[, "cover"]),
          Power = mean(store_original[, "reject"])
        ),
        data.frame(
          Study = "Privacy calibration",
          Method = "Gaussian DP",
          n = n,
          gamma = gamma,
          epsilon = epsilon,
          delta = delta,
          Noise_SD = sd_gaussian_dp,
          Bias = mean(store_gaussian_dp[, "beta_hat"]) - beta1,
          SD = sd(store_gaussian_dp[, "beta_hat"]),
          MSE = mean((store_gaussian_dp[, "beta_hat"] - beta1)^2),
          Coverage = mean(store_gaussian_dp[, "cover"]),
          Power = mean(store_gaussian_dp[, "reject"])
        ),
        data.frame(
          Study = "Privacy calibration",
          Method = "Privacy-calibrated RRM",
          n = n,
          gamma = gamma,
          epsilon = epsilon,
          delta = delta,
          Noise_SD = sd_rrm_dp,
          Bias = mean(store_rrm_dp[, "beta_hat"]) - beta1,
          SD = sd(store_rrm_dp[, "beta_hat"]),
          MSE = mean((store_rrm_dp[, "beta_hat"] - beta1)^2),
          Coverage = mean(store_rrm_dp[, "cover"]),
          Power = mean(store_rrm_dp[, "reject"])
        )
      )

      study4_results <- rbind(study4_results, tmp)
    }
  }
}

############################################################
## Study V: High-Dimensional and Non-Gaussian Settings
############################################################

generate_error <- function(n, type = "normal", eta = 0.1) {
  if (type == "normal") {
    rnorm(n)
  } else if (type == "t3") {
    rt(n, df = 3) / sqrt(3)
  } else if (type == "skewed") {
    rchisq(n, df = 3) - 3
  } else if (type == "contaminated") {
    z <- rbinom(n, 1, eta)
    rnorm(n, 0, ifelse(z == 1, 5, 1))
  } else {
    stop("Unknown error type")
  }
}

study5_results <- data.frame()

p_vec <- c(10, 50, 100)
error_types <- c("normal", "t3", "skewed", "contaminated")

for (n in c(250, 500)) {
  for (p in p_vec) {
    for (gamma in gamma_vec) {
      for (etype in error_types) {

        beta <- rep(0, p)
        beta[1:5] <- c(1, 0.8, 0.6, 0.4, 0.2)

        mse_original <- mse_rrm <- numeric(R)
        pred_original <- pred_rrm <- numeric(R)

        for (r in 1:R) {
          X <- matrix(rnorm(n * p), n, p)
          mu <- as.vector(X %*% beta)
          e <- generate_error(n, etype)
          y <- mu + e

          y_star <- moment_rrm(y, mu, gamma, sigma = sd(e))

          fit_original <- lm(y ~ X)
          fit_rrm <- lm(y_star ~ X)

          beta_hat_original <- coef(fit_original)[-1]
          beta_hat_rrm <- coef(fit_rrm)[-1]

          mse_original[r] <- mean((beta_hat_original - beta)^2)
          mse_rrm[r] <- mean((beta_hat_rrm - beta)^2)

          pred_original[r] <- mean((y - fitted(fit_original))^2)
          pred_rrm[r] <- mean((y_star - fitted(fit_rrm))^2)
        }

        study5_results <- rbind(
          study5_results,
          data.frame(
            Study = "High-dimensional non-Gaussian",
            n = n,
            p = p,
            gamma = gamma,
            Error = etype,
            MSE_Original = mean(mse_original),
            MSE_RRM = mean(mse_rrm),
            Prediction_Error_Original = mean(pred_original),
            Prediction_Error_RRM = mean(pred_rrm)
          )
        )
      }
    }
  }
}

############################################################
## Save Results
############################################################

write.csv(study1_results, "rrm_study1_inferential_preservation.csv", row.names = FALSE)
write.csv(study2_results, "rrm_study2_robustness_influence.csv", row.names = FALSE)
write.csv(study3_results, "rrm_study3_geometric_diagnostics.csv", row.names = FALSE)
write.csv(study4_results, "rrm_study4_privacy_calibration.csv", row.names = FALSE)
write.csv(study5_results, "rrm_study5_highdim_nongaussian.csv", row.names = FALSE)

############################################################
## Print Key Results
############################################################

cat("\nStudy I: Inferential Preservation\n")
print(head(study1_results, 12))

cat("\nStudy II: Robustness and Influence\n")
print(head(study2_results, 12))

cat("\nStudy III: Geometric Diagnostics\n")
print(head(study3_results, 12))

cat("\nStudy IV: Privacy Calibration\n")
print(head(study4_results, 12))

cat("\nStudy V: High-Dimensional and Non-Gaussian\n")
print(head(study5_results, 12))
