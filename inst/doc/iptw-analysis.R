## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 5,
  warning = FALSE,
  message = FALSE
)

library(riskdiff)

## ----basic_example------------------------------------------------------------
# Load the data
data(cachar_sample)

# Quick look at the data
head(cachar_sample)
table(cachar_sample$areca_nut, cachar_sample$abnormal_screen)

## ----calculate_weights--------------------------------------------------------
# Calculate ATE weights for areca nut use
iptw_result <- calc_iptw_weights(
  data = cachar_sample,
  treatment = "areca_nut", 
  covariates = c("age", "sex", "residence", "smoking", "tobacco_chewing"),
  weight_type = "ATE",
  verbose = TRUE
)

# Examine the results
print(iptw_result)

## ----check_assumptions--------------------------------------------------------
# Comprehensive assumption checking
assumptions <- check_iptw_assumptions(iptw_result, verbose = TRUE)

## ----balance_plots, fig.width=8, fig.height=6---------------------------------
# Create balance plots (requires ggplot2)
if (requireNamespace("ggplot2", quietly = TRUE)) {
  plots <- create_balance_plots(iptw_result, plot_type = "both")
  print(plots$love_plot)
  print(plots$ps_plot)
}

## ----causal_effect------------------------------------------------------------
# Estimate ATE using IPTW
rd_causal <- calc_risk_diff_iptw(
  data = cachar_sample,
  outcome = "abnormal_screen",
  treatment = "areca_nut",
  covariates = c("age", "sex", "residence", "smoking", "tobacco_chewing"),
  weight_type = "ATE",
  verbose = TRUE
)

print(rd_causal)
summary(rd_causal)

## ----ate_example--------------------------------------------------------------
rd_ate <- calc_risk_diff_iptw(
  data = cachar_sample,
  outcome = "abnormal_screen", 
  treatment = "areca_nut",
  covariates = c("age", "sex", "residence", "smoking"),
  weight_type = "ATE"
)

cat("ATE: The average causal effect of areca nut use in the population\n")
cat("Risk Difference:", scales::percent(rd_ate$rd_iptw, accuracy = 0.01), "\n")

## ----att_example--------------------------------------------------------------
rd_att <- calc_risk_diff_iptw(
  data = cachar_sample,
  outcome = "abnormal_screen",
  treatment = "areca_nut", 
  covariates = c("age", "sex", "residence", "smoking"),
  weight_type = "ATT"
)

cat("ATT: The average causal effect among areca nut users\n")
cat("Risk Difference:", scales::percent(rd_att$rd_iptw, accuracy = 0.01), "\n")

## ----atc_example--------------------------------------------------------------
rd_atc <- calc_risk_diff_iptw(
  data = cachar_sample,
  outcome = "abnormal_screen",
  treatment = "areca_nut",
  covariates = c("age", "sex", "residence", "smoking"), 
  weight_type = "ATC"
)

cat("ATC: The average causal effect among non-users of areca nut\n")
cat("Risk Difference:", scales::percent(rd_atc$rd_iptw, accuracy = 0.01), "\n")

## ----bootstrap_example--------------------------------------------------------
rd_bootstrap <- calc_risk_diff_iptw(
  data = cachar_sample,
  outcome = "head_neck_abnormal",
  treatment = "tobacco_chewing",
  covariates = c("age", "sex", "residence", "areca_nut"),
  bootstrap_ci = TRUE,
  boot_n = 500,  # Use more in practice (1000+)
  verbose = FALSE
)

print(rd_bootstrap)

## ----ps_models----------------------------------------------------------------
# Logistic regression (default)
ps_logit <- calc_iptw_weights(
  data = cachar_sample,
  treatment = "tobacco_chewing",
  covariates = c("age", "sex", "residence", "areca_nut"),
  method = "logistic"
)

# Probit regression
ps_probit <- calc_iptw_weights(
  data = cachar_sample,
  treatment = "tobacco_chewing", 
  covariates = c("age", "sex", "residence", "areca_nut"),
  method = "probit"
)

# Compare propensity score distributions
cat("Logistic PS range:", round(range(ps_logit$ps), 3), "\n")
cat("Probit PS range:", round(range(ps_probit$ps), 3), "\n")

## ----stabilization------------------------------------------------------------
# Unstabilized weights
ps_unstab <- calc_iptw_weights(
  data = cachar_sample,
  treatment = "areca_nut",
  covariates = c("age", "sex", "residence"),
  stabilize = FALSE
)

# Stabilized weights (default)
ps_stab <- calc_iptw_weights(
  data = cachar_sample,
  treatment = "areca_nut",
  covariates = c("age", "sex", "residence"), 
  stabilize = TRUE
)

cat("Unstabilized weight variance:", round(var(ps_unstab$weights), 2), "\n")
cat("Stabilized weight variance:", round(var(ps_stab$weights), 2), "\n")

## ----trimming-----------------------------------------------------------------
# Check for extreme weights
summary(ps_stab$weights)

# Trim at 1st and 99th percentiles
ps_trimmed <- calc_iptw_weights(
  data = cachar_sample,
  treatment = "areca_nut",
  covariates = c("age", "sex", "residence"),
  trim_weights = TRUE,
  trim_quantiles = c(0.01, 0.99)
)

cat("Original weight range:", round(range(ps_stab$weights), 2), "\n")
cat("Trimmed weight range:", round(range(ps_trimmed$weights), 2), "\n")

## ----comparison---------------------------------------------------------------
# Traditional regression-based risk difference
rd_regression <- calc_risk_diff(
  data = cachar_sample,
  outcome = "abnormal_screen",
  exposure = "areca_nut",
  adjust_vars = c("age", "sex", "residence", "smoking"),
  link = "auto"
)

# IPTW-based causal risk difference  
rd_iptw <- calc_risk_diff_iptw(
  data = cachar_sample,
  outcome = "abnormal_screen",
  treatment = "areca_nut",
  covariates = c("age", "sex", "residence", "smoking"),
  weight_type = "ATE"
)

# Compare results
comparison_table <- data.frame(
  Method = c("Regression Adjustment", "IPTW (ATE)"),
  Risk_Difference = scales::percent(c(rd_regression$rd, rd_iptw$rd_iptw), accuracy = 0.01),
  CI_Lower = scales::percent(c(rd_regression$ci_lower, rd_iptw$ci_lower), accuracy = 0.01),
  CI_Upper = scales::percent(c(rd_regression$ci_upper, rd_iptw$ci_upper), accuracy = 0.01),
  P_Value = sprintf("%.3f", c(rd_regression$p_value, rd_iptw$p_value))
)

print(comparison_table)

## ----poor_balance, eval=FALSE-------------------------------------------------
# # Check which variables have poor balance
# assumptions <- check_iptw_assumptions(iptw_result)
# poor_balance_vars <- assumptions$balance$poor_balance_vars
# 
# if (length(poor_balance_vars) > 0) {
#   cat("Variables with poor balance:", paste(poor_balance_vars, collapse = ", "), "\n")
# 
#   # Try including interactions or polynomial terms
#   iptw_improved <- calc_iptw_weights(
#     data = cachar_sample,
#     treatment = "areca_nut",
#     covariates = c("age", "I(age^2)", "sex", "residence",
#                    "smoking", "age:sex"),  # Add interactions
#     weight_type = "ATE"
#   )
# }

## ----extreme_ps, eval=FALSE---------------------------------------------------
# # Check propensity score distribution
# iptw_result <- calc_iptw_weights(
#   data = cachar_sample,
#   treatment = "areca_nut",
#   covariates = c("age", "sex", "residence")
# )
# 
# # Identify subjects with extreme scores
# extreme_low <- which(iptw_result$ps < 0.05)
# extreme_high <- which(iptw_result$ps > 0.95)
# 
# if (length(extreme_low) > 0 || length(extreme_high) > 0) {
#   cat("Consider trimming sample to region of common support\n")
# 
#   # Restrict to common support
#   common_support <- iptw_result$ps >= 0.05 & iptw_result$ps <= 0.95
#   data_restricted <- cachar_sample[common_support, ]
# 
#   # Re-analyze with restricted sample
#   rd_restricted <- calc_risk_diff_iptw(
#     data = data_restricted,
#     outcome = "abnormal_screen",
#     treatment = "areca_nut",
#     covariates = c("age", "sex", "residence")
#   )
# }

## ----model_spec, eval=FALSE---------------------------------------------------
# # Simple model
# ps_simple <- calc_iptw_weights(
#   data = cachar_sample,
#   treatment = "areca_nut",
#   covariates = c("age", "sex")
# )
# 
# # Complex model with interactions
# ps_complex <- calc_iptw_weights(
#   data = cachar_sample,
#   treatment = "areca_nut",
#   covariates = c("age", "I(age^2)", "sex", "residence",
#                  "smoking", "tobacco_chewing", "age:sex")
# )
# 
# # Compare balance
# check_iptw_assumptions(ps_simple, verbose = FALSE)
# check_iptw_assumptions(ps_complex, verbose = FALSE)

## ----sensitivity, eval=FALSE--------------------------------------------------
# # Simulate an unmeasured confounder
# set.seed(123)
# cachar_sample$unmeasured_confounder <- rbinom(nrow(cachar_sample), 1, 0.3)
# 
# # Compare results with and without the unmeasured confounder
# rd_without_u <- calc_risk_diff_iptw(
#   data = cachar_sample,
#   outcome = "abnormal_screen",
#   treatment = "areca_nut",
#   covariates = c("age", "sex", "residence")
# )
# 
# rd_with_u <- calc_risk_diff_iptw(
#   data = cachar_sample,
#   outcome = "abnormal_screen",
#   treatment = "areca_nut",
#   covariates = c("age", "sex", "residence", "unmeasured_confounder")
# )
# 
# cat("Without unmeasured confounder:", scales::percent(rd_without_u$rd_iptw), "\n")
# cat("With unmeasured confounder:", scales::percent(rd_with_u$rd_iptw), "\n")
# cat("Difference:", scales::percent(abs(rd_without_u$rd_iptw - rd_with_u$rd_iptw)), "\n")

