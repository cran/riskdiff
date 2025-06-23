## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
 collapse = TRUE,
 comment = "#>",
 fig.width = 7,
 fig.height = 4,
 fig.align = "center"
)

## ----load libraries-----------------------------------------------------------
library(riskdiff)
library(dplyr)
library(ggplot2)

## ----load example data--------------------------------------------------------
data(cachar_sample)

# Basic structure
glimpse(cachar_sample)

# Key variables for our analysis
cachar_sample %>%
  select(abnormal_screen, areca_nut, smoking, alcohol, age, sex, residence) %>%
  summary()

## ----simple unadjusted analysis-----------------------------------------------
# Simple unadjusted analysis
rd_simple <- calc_risk_diff(
  data = cachar_sample,
  outcome = "abnormal_screen",
  exposure = "areca_nut"
)

print(rd_simple)

## ----visualisation------------------------------------------------------------
# Let's visualize what this risk difference means
exposure_summary <- cachar_sample %>%
  group_by(areca_nut) %>%
  summarise(
    n = n(),
    cases = sum(abnormal_screen),
    risk = mean(abnormal_screen),
    se = sqrt(risk * (1 - risk) / n)
  ) %>%
  mutate(
    risk_percent = risk * 100,
    se_percent = se * 100
  )

print(exposure_summary)

# The risk difference is simply:
rd_value <- diff(exposure_summary$risk)
cat("Risk difference:", round(rd_value * 100, 1), "percentage points\n")

## ----adjusted analysis--------------------------------------------------------
# Age and sex adjusted analysis
rd_adjusted <- calc_risk_diff(
  data = cachar_sample,
  outcome = "abnormal_screen",
  exposure = "areca_nut",
  adjust_vars = c("age", "sex")
)
print(rd_adjusted)

# Show comparison in a readable format
cat("\n=== COMPARISON: UNADJUSTED vs ADJUSTED ===\n")
cat("Unadjusted Risk Difference:", sprintf("%.2f%%", rd_simple$rd * 100), 
    sprintf("(%.2f%%, %.2f%%)", rd_simple$ci_lower * 100, rd_simple$ci_upper * 100), "\n")
cat("Adjusted Risk Difference:  ", sprintf("%.2f%%", rd_adjusted$rd * 100), 
    sprintf("(%.2f%%, %.2f%%)", rd_adjusted$ci_lower * 100, rd_adjusted$ci_upper * 100), "\n")
cat("Difference in estimates:   ", sprintf("%.2f%%", (rd_adjusted$rd - rd_simple$rd) * 100), "percentage points\n")

## ----stratified analysis------------------------------------------------------
# Stratified by residence (urban vs rural)
rd_stratified <- calc_risk_diff(
  data = cachar_sample,
  outcome = "abnormal_screen",
  exposure = "areca_nut",
  strata = "residence"
)

print(rd_stratified)

## ----visualise stratified results---------------------------------------------
# Create a forest plot of risk differences
plot_data <- rd_stratified %>%
  mutate(
    label = paste0(residence, "\n(n=", n_obs, ")"),
    rd_percent = rd * 100,
    ci_lower_percent = ci_lower * 100,
    ci_upper_percent = ci_upper * 100
  )

ggplot(plot_data, aes(x = rd_percent, y = label)) +
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin = ci_lower_percent, xmax = ci_upper_percent), 
                 height = 0.2) +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
  labs(
    title = "Risk Difference for Cancer by Betel Nut Chewing",
    subtitle = "Stratified by Urban/Rural Residence",
    x = "Risk Difference (percentage points)",
    y = ""
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.y = element_text(size = 11)
  )

## ----convergence handling-----------------------------------------------------
# Force different link functions
rd_identity <- calc_risk_diff(
  cachar_sample, "abnormal_screen", "areca_nut", 
  link = "identity"
)

rd_log <- calc_risk_diff(
  cachar_sample, "abnormal_screen", "areca_nut", 
  link = "log"
)

# Compare model types used
cat("Identity link model type:", rd_identity$model_type, "\n")
cat("Log link model type:", rd_log$model_type, "\n")

## ----rare outcomes------------------------------------------------------------
# Create a rare outcome (1% prevalence)
cachar_sample$rare_outcome <- rbinom(nrow(cachar_sample), 1, 0.01)

rd_rare <- calc_risk_diff(
  cachar_sample, 
  "rare_outcome", 
  "areca_nut"
)

print(rd_rare)

## ----missing data-------------------------------------------------------------
# Create a copy with some missing data for demonstration
set.seed(123)  # For reproducibility
cachar_with_missing <- cachar_sample %>%
  mutate(
    # Introduce more modest missing data (~3% in age, ~2% in alcohol)
    age_with_missing = ifelse(runif(n()) < 0.03, NA, age),
    alcohol_with_missing = ifelse(runif(n()) < 0.02, NA, alcohol)
  )

# Check the missing data patterns
missing_summary <- cachar_with_missing %>%
  summarise(
    total_observations = n(),
    age_missing = sum(is.na(age_with_missing)),
    alcohol_missing = sum(is.na(alcohol_with_missing)),
    total_missing_any = sum(!complete.cases(select(., age_with_missing, alcohol_with_missing, abnormal_screen, areca_nut))),
    complete_cases = sum(complete.cases(select(., age_with_missing, alcohol_with_missing, abnormal_screen, areca_nut)))
  )

print(missing_summary)

# Analysis with variables that have missing data
rd_missing <- calc_risk_diff(
  cachar_with_missing,
  "abnormal_screen",
  "areca_nut",
  adjust_vars = c("age_with_missing", "alcohol_with_missing")
)

# Compare with complete case analysis
rd_complete <- calc_risk_diff(
  cachar_sample,
  "abnormal_screen", 
  "areca_nut",
  adjust_vars = c("age", "alcohol")
)

cat("\n=== IMPACT OF MISSING DATA ===\n")
cat("Complete data analysis (n=", rd_complete$n_obs, "): ", sprintf("%.2f%%", rd_complete$rd * 100), 
    sprintf(" (%.2f%%, %.2f%%)", rd_complete$ci_lower * 100, rd_complete$ci_upper * 100), "\n")

# Check if missing data analysis succeeded
if (!is.na(rd_missing$rd)) {
  cat("Missing data analysis (n=", rd_missing$n_obs, "):  ", sprintf("%.2f%%", rd_missing$rd * 100), 
      sprintf(" (%.2f%%, %.2f%%)", rd_missing$ci_lower * 100, rd_missing$ci_upper * 100), "\n")
  cat("Cases lost to missing data: ", rd_complete$n_obs - rd_missing$n_obs, "\n")
} else {
  cat("Missing data analysis: FAILED (insufficient data or convergence issues)\n")
  cat("Attempted to use n =", rd_missing$n_obs, "complete cases\n")
  cat("Cases lost to missing data: ", nrow(cachar_with_missing) - rd_missing$n_obs, "\n\n")
  
  cat("📚 LESSON: This demonstrates why missing data can be problematic:\n")
  cat("   • Listwise deletion can dramatically reduce sample size\n")
  cat("   • Small samples may cause model convergence failures\n") 
  cat("   • Consider multiple imputation for better missing data handling\n")
  cat("   • The riskdiff package gracefully handles these failures\n")
}

## ----basic syntax example, eval=FALSE-----------------------------------------
# # Example usage:
# result <- calc_risk_diff(
#   data = cachar_sample,           # Your dataset
#   outcome = "abnormal_screen",    # Binary outcome variable (0/1)
#   exposure = "areca_nut",         # Exposure of interest
#   adjust_vars = c("age", "sex"),  # Variables to adjust for
#   strata = "residence",           # Stratification variables
#   link = "auto",                  # Link function: "auto", "identity", "log", "logit"
#   alpha = 0.05,                   # Significance level (0.05 = 95% CI)
#   verbose = FALSE                 # Print diagnostic messages if TRUE
# )

