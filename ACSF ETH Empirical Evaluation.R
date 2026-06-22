# ==============================================================================
# SETUP & DATA 
# ==============================================================================
library(quantmod)
library(urca)
library(forecast)
# Ethereum Dataset
getSymbols("ETH-USD", src = "yahoo", from = "2020-01-01", to = "2026-06-19")
eth_adj <- na.omit(Cl(get("ETH-USD")))
colnames(eth_adj) <- "ETH_Price"
# Augmented Dickey-Fuller validation
adf_joint <- ur.df(eth_adj, type = "trend", selectlags = "AIC")
test_stats <- adf_joint@teststat
crit_vals  <- adf_joint@cval
tau3_stat   <- test_stats[1]
tau3_crit5  <- crit_vals[1, "5pct"] 
phi3_stat   <- test_stats[3]
phi3_crit5  <- crit_vals[3, "5pct"]
# Output ADF Summary
cat("Tau3 Unit Root Test Statistic :", round(tau3_stat, 4), " (5% Crit:", tau3_crit5, ")\n")
cat("Phi3 Joint Trend Test Statistic:", round(phi3_stat, 4), " (5% Crit:", phi3_crit5, ")\n")
if (tau3_stat > tau3_crit5) {
  cat("Non-Stationary. Log-Difference Transformation.\n")
  eth_transformed <- diff(log(eth_adj))[-1]
} else {
  if (phi3_stat > phi3_crit5) {
    cat("Trend-Stationary. Residual Detrending.\n")
    time_idx <- 1:length(eth_adj)
    eth_transformed <- residuals(lm(as.numeric(eth_adj) ~ time_idx))
  } else {
    cat("Stationary Level Series. Base Adjusted Prices.\n")
    eth_transformed <- eth_adj
  }
}
eth_transformed <- na.omit(eth_transformed)
colnames(eth_transformed) <- "Stationary Data"
# Stationary verification ADF test
adf_ret  <- ur.df(eth_transformed, type = "drift", selectlags = "AIC")
tau2_ret <- adf_ret@teststat[1]
# ==============================================================================
# Table 1 / Plot 1 - Descriptive Stats and Stationarity Tests
# ==============================================================================
eth_report_table <- data.frame(
  Metric = c("Observations (N)", "Mean", "Variance", "ADF Unit Root Test", "Joint Trend Test"),
  Raw_Prices = c(
    length(eth_adj),
    round(mean(eth_adj), 4),
    round(var(eth_adj), 4),
    paste0(round(tau3_stat, 3), " [CV 5%: -3.43]"),
    paste0(round(phi3_stat, 3), " [CV 5%:  6.25]")
  ),
  Log_Returns = c(
    length(eth_transformed),
    round(mean(eth_transformed * 100), 6),
    round(var(eth_transformed * 100), 6),
    paste0(round(tau2_ret, 3), " [CV 5%: -2.86]"),
    "N/A"
  )
)
print(eth_report_table, row.names = FALSE)
#Price & Return Graphs
par(mfrow = c(2, 1), mar = c(4, 4, 2, 2))
plot(eth_adj, main = "Ethereum Closing Price", col = "darkblue", xlab = "Date", major.ticks = "years", grid.col = "lightgray")
plot(eth_transformed, main = "Ethereum Daily Log Returns", col = "darkgreen", xlab = "Date", major.ticks = "years", grid.col = "lightgray")
# ==============================================================================
# GAS FRAMEWORK
# ==============================================================================
r <- as.numeric(eth_transformed) * 100 
N_total    <- length(r)
train_size <- floor(0.70 * N_total)
r_in  <- r[1:train_size]
r_out <- r[(train_size + 1):N_total]
N_test  <- N_total - train_size
c_90 <- as.numeric(quantile(abs(r_in), 0.90))
c_95 <- as.numeric(quantile(abs(r_in), 0.95))
c_97 <- as.numeric(quantile(abs(r_in), 0.97))
c_99 <- as.numeric(quantile(abs(r_in), 0.99))
thresholds <- c(c_90, c_95, c_97, c_99)
threshold_labels <- c("Q90", "Q95","Q97","Q99")
diks_censored_score <- function(r_t, mu_t, sig2_t, c_val) {
  sigma_t <- sqrt(sig2_t)
  res     <- r_t - mu_t
  if (abs(r_t) <= c_val) {
    s_m <- res
    s_v <- (res^2 / sig2_t) - 1
  } else {
    z_plus  <- (c_val - mu_t) / sigma_t
    z_minus <- (c_val + mu_t) / sigma_t
    G_f <- 2 - pnorm(z_plus) - pnorm(z_minus)
    G_f <- max(G_f, 1e-10)
    s_m <- sigma_t * (dnorm(z_plus) - dnorm(z_minus)) / G_f
    s_v <- (z_plus * dnorm(z_plus) + z_minus * dnorm(z_minus)) / G_f
  }
  list(s_m = s_m, s_v = s_v)
}
#GAS update equation
run_gas_filter <- function(theta, r_data, type = "gaussian", c_val = NULL, alpha_val = NULL) {
  omega_m <- theta[1]; A_m <- theta[2]; B_m <- theta[3]
  omega_v <- theta[4]; A_v <- theta[5]; B_v <- theta[6]
  T_len   <- length(r_data)
  mu <- numeric(T_len); lambda <- numeric(T_len) 
  mu[1] <- mean(r_data); lambda[1] <- log(var(r_data))
  if (is.null(alpha_val)) alpha_val <- 0
  for(t in 1:(T_len - 1)) {
    sig2 <- exp(lambda[t]); res <- r_data[t] - mu[t]
    s_m_gauss <- res; s_v_gauss <- (res^2 / sig2) - 1
    if (type == "gaussian") {
      s_m <- s_m_gauss; s_v <- s_v_gauss
    } else if (type == "censored") {
      cens <- diks_censored_score(r_data[t], mu[t], sig2, c_val)
      s_m <- cens$s_m; s_v <- cens$s_v
    } else if (type == "convex") {
      cens <- diks_censored_score(r_data[t], mu[t], sig2, c_val)
      s_m <- alpha_val * s_m_gauss + (1 - alpha_val) * cens$s_m
      s_v <- alpha_val * s_v_gauss + (1 - alpha_val) * cens$s_v
    }
    mu[t+1]     <- omega_m + A_m * s_m + B_m * mu[t]
    lambda[t+1] <- omega_v + A_v * s_v + B_v * lambda[t]
  }
  return(list(mu = mu, lambda = lambda))
}
gas_objective <- function(theta, r_data, type = "gaussian", c_val = NULL, alpha_input = NULL) {
  if(abs(theta[3]) >= 0.999 || abs(theta[6]) >= 0.999) return(1e10)
  alpha_val <- if (type == "convex") 1 / (1 + exp(-alpha_input)) else NULL
  paths     <- run_gas_filter(theta, r_data, type, c_val, alpha_val)
  sig2      <- exp(paths$lambda)
  nll       <- 0.5 * sum(log(2 * pi) + paths$lambda + ((r_data - paths$mu)^2 / sig2))
  if(is.nan(nll) || is.na(nll) || is.infinite(nll)) return(1e10)
  return(nll)
}
calc_crps_vec     <- function(actual, mu, sigma) { z <- (actual - mu) / sigma; return(sigma * (z * (2 * pnorm(z) - 1) + 2 * dnorm(z) - (1 / sqrt(pi)))) }
calc_logscore_vec <- function(actual, mu, sigma) { return(-0.5 * (log(2 * pi) + 2 * log(sigma) + ((actual - mu)^2 / (sigma^2)))) }
add_stars <- function(stat, pval) {
  if (is.na(pval)) return(as.character(stat))
  if (pval < 0.01) return(paste0(stat, "***"))
  if (pval < 0.05) return(paste0(stat, "**"))
  if (pval < 0.10) return(paste0(stat, "*"))
  return(as.character(stat))
}
# ==============================================================================
# IN-SAMPLE TRAINING
# ==============================================================================
results_table <- data.frame(Architecture = character(), NLL = numeric(), AIC = numeric(), BIC = numeric(), stringsAsFactors = FALSE)
initial_theta <- c(mean(r_in)*0.1, 0.05, 0.90, log(var(r_in))*0.1, 0.05, 0.9)
fit_g <- optim(initial_theta, gas_objective, r_data = r_in, type = "gaussian", method = "BFGS")
results_table <- rbind(results_table, data.frame(Architecture = "Log Likelihood", 
                                                 NLL = round(-fit_g$value, 4),
                                                 AIC = round(2*6 + 2*fit_g$value, 4), BIC = round(6*log(train_size) + 2*fit_g$value, 4)
))
for(i in 1:length(thresholds)) {
  current_c <- thresholds[i]
  lbl   <- threshold_labels[i]
  fit_cens <- optim(initial_theta, gas_objective, r_data = r_in, type = "censored", 
                    c_val = current_c, method = "BFGS", control = list(maxit = 500))
  nll_c <- fit_cens$value
  k_c   <- length(initial_theta)
  aic_c <- 2 * k_c + 2 * nll_c
  bic_c <- k_c * log(train_size) + 2 * nll_c
  results_table <- rbind(results_table, data.frame(
    Architecture = paste0("Censored @", lbl),
    NLL = -nll_c, AIC = aic_c, BIC = bic_c
  ))
  #alpha = 1 / 1 + exp(-1)
  initial_convex_theta <- c(initial_theta, 1)
  # Utilize L-BFGS-B
  fit_conv <- optim(
    par = initial_convex_theta, 
    fn = function(p) {
      gas_objective(theta = p[1:6], r_data = r_in, type = "convex", c_val = current_c, alpha_input = p[7])
    }, 
    method = "L-BFGS-B",
    lower = c(-Inf, 0.0001, -0.999, -Inf, 0.0001, -0.999, -Inf),
    upper = c( Inf, 0.9999,  0.999,  Inf, 0.9999,  0.999,  Inf),
    control = list(maxit = 1000)
  )
  alpha_star <- 1 / (1 + exp(-fit_conv$par[7]))
  nll_v <- fit_conv$value
  k_v   <- length(initial_convex_theta)
  aic_v <- 2 * k_v + 2 * nll_v
  bic_v <- k_v * log(train_size) + 2 * nll_v
  results_table <- rbind(results_table, data.frame(Architecture = paste0("ACSF(",sprintf("%.4f", alpha_star) ,") @", lbl), 
                                                   NLL = -nll_v, AIC = aic_v, BIC = bic_v
  ))
}
print(results_table, row.names = FALSE)
# ==============================================================================
# SCORE IMPACT FUNCTIONS 
# ==============================================================================
library(ggplot2)
library(tidyr)
u_grid <- seq(-4, 4, length.out = 1000)
alphas_est <- numeric(length(thresholds))
for(j in 1:length(thresholds)) {
  cc <- thresholds[j]
  initial_convex_theta <- c(initial_theta, 1)
  fit_conv <- optim(
    par = initial_convex_theta, 
    fn = function(p) { gas_objective(theta = p[1:6], r_data = r_in, type = "convex", c_val = cc, alpha_input = p[7]) }, 
    method = "L-BFGS-B",
    lower = c(-Inf, 0.0001, -0.999, -Inf, 0.0001, -0.999, -Inf),
    upper = c( Inf, 0.9999,  0.999,  Inf, 0.9999,  0.999,  Inf),
    control = list(maxit = 1000)
  )
  alphas_est[j] <- 1 / (1 + exp(-fit_conv$par[7]))
}
df_list_mean <- list()
df_list_var  <- list()
df_list_vlines <- list()
dynamic_labels <- character(length(thresholds))
for(j in 1:length(thresholds)) {
  cc <- thresholds[j]
  alpha_val <- alphas_est[j]
  lbl <- paste0(threshold_labels[j], " (\u03b1 = ", round(alpha_val, 3), ")")
  dynamic_labels[j] <- lbl
  std_C <- cc / sd(r_in)
  s_m_gauss <- u_grid
  s_m_cens  <- ifelse(abs(u_grid) <= std_C, u_grid, 0)
  s_m_conv  <- alpha_val * s_m_gauss + (1 - alpha_val) * s_m_cens
  s_v_gauss <- u_grid^2 - 1
  mills_ratio <- std_C * dnorm(std_C) / (1 - pnorm(std_C))
  s_v_cens  <- ifelse(abs(u_grid) <= std_C, u_grid^2 - 1, mills_ratio)
  s_v_conv  <- alpha_val * s_v_gauss + (1 - alpha_val) * s_v_cens
  df_list_mean[[j]] <- data.frame(u = rep(u_grid, 3), Score = c(s_m_gauss, s_m_cens, s_m_conv),
                                  Model = rep(c("Log-likelihood", "Censored", "ACSF"), each = 1000), Quantile = lbl)
  df_list_var[[j]]  <- data.frame(u = rep(u_grid, 3), Score = c(s_v_gauss, s_v_cens, s_v_conv),
                                  Model = rep(c("Log-likelihood", "Censored", "ACSF"), each = 1000), Quantile = lbl)
  df_list_vlines[[j]] <- data.frame(Quantile = lbl, std_C = std_C)
}
df_mean_all   <- do.call(rbind, df_list_mean)
df_var_all    <- do.call(rbind, df_list_var)
df_vlines_all <- do.call(rbind, df_list_vlines)
df_mean_all$Quantile   <- factor(df_mean_all$Quantile, levels = dynamic_labels)
df_var_all$Quantile    <- factor(df_var_all$Quantile, levels = dynamic_labels)
df_vlines_all$Quantile <- factor(df_vlines_all$Quantile, levels = dynamic_labels)
palette_colors <- c("Log-likelihood" = "#555555", "Censored" = "#E69F00", "ACSF" = "#0072B2")
# ==============================================================================
# Plot 2: LOCATION IMPACT 
# ==============================================================================
ggplot(df_mean_all, aes(x = u, y = Score, color = Model, linetype = Model)) +
  geom_vline(data = df_vlines_all, aes(xintercept = std_C), linetype = "dotted", color = "gray60") +
  geom_vline(data = df_vlines_all, aes(xintercept = -std_C), linetype = "dotted", color = "gray60") +
  geom_line(linewidth = 0.6) +
  facet_wrap(~ Quantile, nrow = 2) +
  scale_color_manual(values = palette_colors) +
  scale_linetype_manual(values = c("Log-likelihood" = "solid", "Censored" = "solid", "ACSF" = "solid")) +
  labs(
    x = "Standardized Innovation",
    y = "Score Weight"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "gray95"),
    strip.text = element_text(face = "bold", size = 10, color = "#333333")
  )
# ==============================================================================
# Plot 3: SCALE IMPACT
# ==============================================================================
ggplot(df_var_all, aes(x = u, y = Score, color = Model, linetype = Model)) +
  geom_vline(data = df_vlines_all, aes(xintercept = std_C), linetype = "dotted", color = "gray60") +
  geom_vline(data = df_vlines_all, aes(xintercept = -std_C), linetype = "dotted", color = "gray60") +
  geom_line(linewidth = 0.5) +
  facet_wrap(~ Quantile, nrow = 2) +
  scale_color_manual(values = palette_colors) +
  scale_linetype_manual(values = c("Log-likelihood" = "solid", "Censored" = "solid", "ACSF" = "solid")) +
  labs(
    x = "Standardized Innovation",
    y = "Score Weight"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "gray95"),
    strip.text = element_text(face = "bold", size = 10, color = "#333333")
  )
# ==============================================================================
# FIXED-PARAMETER ROLLING WINDOW 
# ==============================================================================
rolling_labels <- c("90", "95", "97", "99")
fit_params <- list()
for(i in 1:length(thresholds)) {
  lbl <- rolling_labels[i]
  current_c <- thresholds[i]
  fit_cens <- optim(initial_theta, gas_objective, r_data = r_in, type = "censored", 
                    c_val = current_c, method = "BFGS", control = list(maxit = 500))
  fit_params[[paste0("Censored @", lbl)]] <- fit_cens$par
  fit_conv <- optim(
    par = c(initial_theta, 1), 
    fn = function(p) {
      gas_objective(theta = p[1:6], r_data = r_in, type = "convex", c_val = current_c, alpha_input = p[7])
    }, 
    method = "L-BFGS-B",
    lower = c(-Inf, 0.0001, -0.999, -Inf, 0.0001, -0.999, -Inf),
    upper = c( Inf, 0.9999,  0.999,  Inf, 0.9999,  0.999,  Inf),
    control = list(maxit = 1000)
  )
  fit_params[[paste0("ACSF @", lbl)]] <- fit_conv$par
}
models <- c("Log-likelihood", "Censored @90", "ACSF @90","Censored @95", "ACSF @95","Censored @97",
            "ACSF @97","Censored @99", "ACSF @99")
pred_mu    <- matrix(NA, nrow = N_test, ncol = length(models), dimnames = list(NULL, models))
pred_sigma <- matrix(NA, nrow = N_test, ncol = length(models), dimnames = list(NULL, models))
actual_r     <- numeric(N_test)
realized_var <- numeric(N_test)

for(i in 1:N_test) {
  start_idx <- i; end_idx <- i + train_size - 1
  r_window  <- r[start_idx:end_idx]
  target_r <- r[end_idx + 1]; actual_r[i] <- target_r
  realized_var[i] <- (target_r - mean(r_window))^2
  # Track 1
  path_g <- run_gas_filter(fit_g$par, r_window, type = "gaussian")
  res_g  <- r_window[train_size] - path_g$mu[train_size]
  s_m_g  <- res_g; s_v_g <- (res_g^2 / exp(path_g$lambda[train_size])) - 1
  pred_mu[i, "Log-likelihood"]    <- fit_g$par[1] + fit_g$par[2]*s_m_g + fit_g$par[3]*path_g$mu[train_size]
  pred_sigma[i, "Log-likelihood"] <- sqrt(exp(fit_g$par[4] + fit_g$par[5]*s_v_g + fit_g$par[6]*path_g$lambda[train_size]))
  # Track 2
  for(j in 1:length(thresholds)) {
    cc <- thresholds[j]; lbl <- rolling_labels[j]
    p_c    <- fit_params[[paste0("Censored @", lbl)]]
    path_c <- run_gas_filter(p_c[1:6], r_window, type = "censored", c_val = cc) 
    sig2_c <- exp(path_c$lambda[train_size])
    cens_c <- diks_censored_score(r_window[train_size], path_c$mu[train_size], sig2_c, cc)
    pred_mu[i, paste0("Censored @", lbl)]    <- p_c[1] + p_c[2]*cens_c$s_m + p_c[3]*path_c$mu[train_size]
    pred_sigma[i, paste0("Censored @", lbl)] <- sqrt(exp(p_c[4] + p_c[5]*cens_c$s_v + p_c[6]*path_c$lambda[train_size]))
    p_v    <- fit_params[[paste0("ACSF @", lbl)]]
    a_star <- 1 / (1 + exp(-p_v[7])) 
    path_v <- run_gas_filter(p_v[1:6], r_window, type = "convex", c_val = cc, alpha_val = a_star)
    sig2_v <- exp(path_v$lambda[train_size]); res_v <- r_window[train_size] - path_v$mu[train_size]
    s_m_gau <- res_v; s_v_gau <- (res_v^2 / sig2_v) - 1
    cens_v  <- diks_censored_score(r_window[train_size], path_v$mu[train_size], sig2_v, cc)
    s_m_acsf <- a_star*s_m_gau + (1-a_star)*cens_v$s_m
    s_v_acsf <- a_star*s_v_gau + (1-a_star)*cens_v$s_v
    pred_mu[i, paste0("ACSF @", lbl)]    <- p_v[1] + p_v[2]*s_m_acsf + p_v[3]*path_v$mu[train_size]
    pred_sigma[i, paste0("ACSF @", lbl)] <- sqrt(exp(p_v[4] + p_v[5]*s_v_acsf + p_v[6]*path_v$lambda[train_size]))
  }
}
# ==============================================================================
# TABLE 3: OUT-OF-SAMPLE DENSITY VALIDATION
# ==============================================================================
crps_base     <- calc_crps_vec(actual_r, pred_mu[, "Log-likelihood"], pred_sigma[, "Log-likelihood"])
logscore_base <- calc_logscore_vec(actual_r, pred_mu[, "Log-likelihood"], pred_sigma[, "Log-likelihood"])
density_table <- data.frame(Architecture = models, Mean_CRPS = NA, DM_CRPS_Stat = NA, Mean_LogScore = NA, DM_LogScore_Stat = NA)
for(m in models) {
  c_crps <- calc_crps_vec(actual_r, pred_mu[, m], pred_sigma[, m])
  c_ls   <- calc_logscore_vec(actual_r, pred_mu[, m], pred_sigma[, m])
  density_table[density_table$Architecture == m, "Mean_CRPS"]     <- round(mean(c_crps), 5)
  density_table[density_table$Architecture == m, "Mean_LogScore"] <- round(mean(c_ls), 5)
  if (m != "Log-likelihood") {
    dm_c <- dm.test(crps_base, c_crps, alternative = "greater", h = 1, power = 1)
    dm_l <- dm.test(-logscore_base, -c_ls, alternative = "greater", h = 1, power = 1)
    density_table[density_table$Architecture == m, "DM_CRPS_Stat"]     <- add_stars(round(dm_c$statistic, 3), dm_c$p.value)
    density_table[density_table$Architecture == m, "DM_LogScore_Stat"] <- add_stars(round(dm_l$statistic, 3), dm_l$p.value)
  } else {
    density_table[density_table$Architecture == m, c("DM_CRPS_Stat", "DM_LogScore_Stat")] <- "Ref"
  }
}
print(density_table, row.names = FALSE)
# ==============================================================================
# TABLE 4: VARIANCE ACCURACY & VALUE-AT-RISK
# ==============================================================================
tail_table <- data.frame(Architecture = models, Variance_MSE = NA, VaR97_Rate = NA, Kupiec97_Stat = NA, VaR99_Rate = NA, Kupiec99_Stat = NA)
kupiec_star_calc <- function(actual, mu, sigma, alpha_risk) {
  var_line <- mu + qnorm(alpha_risk) * sigma; violations <- sum(actual < var_line)
  N <- length(actual); p_hat <- violations / N; p <- alpha_risk
  if (violations == 0 || violations == N) return("—")
  lr_stat <- -2 * ((N - violations) * log(1 - p) + violations * log(p) - (N - violations) * log(1 - p_hat) - violations * log(p_hat))
  return(add_stars(round(lr_stat, 3), 1 - pchisq(lr_stat, df = 1)))
}
for(m in models) {
  tail_table[tail_table$Architecture == m, "Variance_MSE"] <- round(mean((realized_var - (pred_sigma[, m]^2))^2), 5)
  # Q95 Boundaries
  v97 <- pred_mu[, m] + qnorm(0.03) * pred_sigma[, m]
  tail_table[tail_table$Architecture == m, "VaR97_Rate"]   <- round(sum(actual_r < v97) / N_test, 4)
  tail_table[tail_table$Architecture == m, "Kupiec97_Stat"] <- kupiec_star_calc(actual_r, pred_mu[, m], pred_sigma[, m], 0.03)
  # Q99 Boundaries
  v99 <- pred_mu[, m] + qnorm(0.01) * pred_sigma[, m]
  tail_table[tail_table$Architecture == m, "VaR99_Rate"]   <- round(sum(actual_r < v99) / N_test, 4)
  tail_table[tail_table$Architecture == m, "Kupiec99_Stat"] <- kupiec_star_calc(actual_r, pred_mu[, m], pred_sigma[, m], 0.01)
}
print(tail_table, row.names = FALSE)
# ==============================================================================
# Plot 4: VaR Price Tracking
# ==============================================================================
library(ggplot2)
library(tidyr)
oos_price_actual <- as.numeric(eth_adj[(train_size + 1):N_total])
oos_price_lagged <- as.numeric(eth_adj[train_size: (N_total - 1)])
floor_loglik <- pred_mu[, "Log-likelihood"] + qnorm(0.01) * pred_sigma[, "Log-likelihood"]
floor_cens97 <- pred_mu[, "Censored @97"]         + qnorm(0.01) * pred_sigma[, "Censored @97"]
floor_acsf97 <- pred_mu[, "ACSF @97"]         + qnorm(0.01) * pred_sigma[, "ACSF @97"]
df_tracking <- data.frame(
  Day            = 1:length(oos_price_actual),
  Actual_Price   = oos_price_actual,
  Log_Likelihood = oos_price_lagged * exp(floor_loglik / 100),
  ACSF_97        = oos_price_lagged * exp(floor_acsf97 / 100),
  Censored_97    = oos_price_lagged * exp(floor_cens97 / 100)
)
df_tracking_long <- tidyr::pivot_longer(
  df_tracking, 
  cols     = c(Actual_Price, Log_Likelihood, ACSF_97, Censored_97),
  names_to = "Series", 
  values_to = "Price"
)
df_tracking_long$Series <- factor(
  df_tracking_long$Series,
  levels = c("Actual_Price", "Log_Likelihood", "ACSF_97", "Censored_97"),
  labels = c("ETH Price", "Log-likelihood", "ACSF @97", "Censored @97")
)
ggplot(df_tracking_long, aes(x = Day, y = Price, color = Series, linetype = Series, size = Series)) +
  geom_line() +
  scale_color_manual(values = c(
    "ETH Price"      = "black",
    "Log-likelihood" = "red",
    "ACSF @97"       = "#0072B2",
    "Censored @97"   = "gold"
  )) +
  scale_linetype_manual(values = c(
    "ETH Price"      = "solid",
    "Log-likelihood" = "solid",
    "ACSF @97"       = "solid",
    "Censored @97"   = "solid"
  )) +
  scale_size_manual(values = c(
    "ETH Price"      = 0.7,
    "Log-likelihood" = 0.8,
    "ACSF @97"       = 0.6,
    "Censored @97"   = 0.4
  )) +
  labs(
    x = "Trading Days",
    y = "Price (USD)"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title        = element_text(face = "bold", size = 13 , hjust = 0.5),
    legend.title      = element_blank(),
    legend.position   = "top",
    legend.box        = "horizontal",
    panel.grid.minor  = element_blank(),
    panel.grid.major  = element_line(color = "gray93")
  )
identical(pred_mu[,"ACSF @99"], pred_mu[,"Log-likelihood"])

