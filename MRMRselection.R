# ============================================================
# MRMR Feature Selection
# Top 100 donor/image features for eGFR outcomes and DGF
# ============================================================

library(tidyverse)
library(caret)
library(mRMRe)

# ------------------------------------------------------------
# Load data
# ------------------------------------------------------------

df <- Reduced_Features_NO_COLOR

# Keep only transplanted kidneys for outcomes
df_transplanted <- df %>%
  filter(Allocated_Discarded == 1)

# ------------------------------------------------------------
# Outcomes
# ------------------------------------------------------------

egfr_outcomes <- c(
  "eGFR_CKD_EPI_1M",
  "eGFR_CKD_EPI_3M",
  "eGFR_CKD_EPI_6M",
  "eGFR_CKD_EPI_12M",
  "eGFR_CKD_EPI_3years",
  "eGFR_CKD_EPI_5years")

dgf_outcome <- "Delayed_Graft_Function"

# ------------------------------------------------------------
# Columns to exclude from predictors
# ------------------------------------------------------------

cols_to_exclude <- c(
  # IDs/admin
  "Allocated_Discarded",
  "Slide_number",
  "Number_of_Sections_Per_Slide",
  "Transplant_date_OR_Harvest_date",
  "Recipient_code",
  "RECIPIENT_DATA",
  
  # recipient variables
  "Recipient_birth_date",
  "Recipient_age",
  "Recipient_Sex",
  "Recipient_sex",
  "Recipient_weight",
  "Recipient_height",
  "Recipient_BMI",
  "Recipient_CKD_ethiology",
  "Country_of_birth",
  
  # follow-up/outcome timing
  "Last_follow_up_date",
  "Transplant_admission_discharge_date",
  "Return_to_dialysis_date",
  "Graft_loss_date",
  "Death_date",
  "Cause_of_death",
  "FollowUp_Status",
  "Follow_up_time_Months",
  
  # all outcomes / post-transplant variables
  "Cr1M", "eGFR_CKD_EPI_1M",
  "Cr3M", "eGFR_CKD_EPI_3M",
  "Cr6M", "eGFR_CKD_EPI_6M",
  "Cr12M", "eGFR_CKD_EPI_12M",
  "Cr3years", "eGFR_CKD_EPI_3years",
  "Cr5years", "eGFR_CKD_EPI_5years",
  "Cr10years", "eGFRMDRD_10years",
  "eGFR_CKD_EPI_last_follow_up",
  "Cr_at_last_Followup",
  "Delayed_Graft_Function",
  
  # derived outcome variables
  "eGFR_Variation_1Y_3Y_total",
  "eGFR_Variation_1Y_5Y_total",
  "eGFR_less30_OR_Dialysis__3Y",
  "eGFR_less30_OR_Dialysis__5Y",
  "eGFR_less_30_dialysis_Decline_sup_5ml__3Y",
  "eGFR_less_30_dialysis_Decline_sup_5ml__5Y",
  "eGFR_12M_ABOVE_50_CKD_EPI",
  "eGFR_12M_ABOVE_55_CKD_EPI",
  "eGFR_12M_ABOVE_60_CKD_EPI",
  
  # graft/death outcomes
  "Graft_loss_Returned_to_Dialysis",
  "Graft_loss",
  "Graft_loss_code",
  "Graft_loss_death_censored",
  "Death_With_a_Functioning_Graft",
  "Death",
  
  # section headers / non-feature labels
  "DONOR_VARIABLES",
  "REMUZZI_CALCULATIONS",
  "AI_FEATURES",
  "COMPUTATIONAL_MORPHOMETRIC",
  "GLOMERULAR_GRANULAR_DATA_REVIEWED_MF01",
  "CT_GRANULAR_FEATURES_REV",
  "GRANULAR_FEATURES_VASCULAR_MEAN_VALUES_ALL_SIZES_rev",
  "GRANULAR_FEATURES_VASCULAR_ARTERY_WITH_MAX_AREA_R",
  "GRANULAR_FEATURES_VASCULAR_ARTERIOLAR_AREA_INF_35000_PX_R",
  "GRANULAR_FEATURES_VASCULAR_ARTERY_AREA_SUP_35000_PX_R")

# ------------------------------------------------------------
# Build clean numeric feature matrix
# ------------------------------------------------------------

features_numeric <- df_transplanted %>%
  select(-any_of(cols_to_exclude)) %>%
  select(where(is.numeric))

# Remove near-zero variance features
nzv <- nearZeroVar(features_numeric)

if (length(nzv) > 0) {
  features_numeric <- features_numeric[, -nzv, drop = FALSE]}

# Remove highly correlated features
cor_matrix <- cor(features_numeric, use = "pairwise.complete.obs")
cor_matrix[is.na(cor_matrix)] <- 0

high_cor <- findCorrelation(cor_matrix, cutoff = 0.80)

if (length(high_cor) > 0) {
  features_numeric <- features_numeric[, -high_cor, drop = FALSE]}

cat("Features remaining after cleaning:", ncol(features_numeric), "\n")

# ------------------------------------------------------------
# MRMR helper function
# ------------------------------------------------------------

run_mrmr <- function(outcome_vector, features, outcome_name, top_n = 100) {
  
  outcome_vector <- as.numeric(outcome_vector)
  valid_rows <- !is.na(outcome_vector)
  
  y <- outcome_vector[valid_rows]
  x <- features[valid_rows, , drop = FALSE]
  
  complete_rows <- complete.cases(x)
  y <- y[complete_rows]
  x <- x[complete_rows, , drop = FALSE]
  
  if (length(y) < 20) {
    warning(paste("Skipping", outcome_name, "- too few complete cases."))
    return(NULL)}
  
  if (ncol(x) < 1) {
    warning(paste("Skipping", outcome_name, "- no usable predictors."))
    return(NULL)}
  
  x <- x %>%
    mutate(across(everything(), ~ suppressWarnings(as.numeric(.))))
  
  x <- x[, colSums(!is.na(x)) > 0, drop = FALSE]
  
  complete_rows <- complete.cases(x) & !is.na(y)
  y <- y[complete_rows]
  x <- x[complete_rows, , drop = FALSE]
  
  if (length(y) < 20 || ncol(x) < 1) {
    warning(paste("Skipping", outcome_name, "- too few complete cases or predictors."))
    return(NULL)}
  
  mrmr_df <- data.frame(
    outcome = as.numeric(y),
    x,
    check.names = FALSE)
  
  mrmr_df[] <- lapply(mrmr_df, function(col) suppressWarnings(as.numeric(col)))
  dd <- mRMR.data(data = mrmr_df)
  
  feature_count <- min(top_n, ncol(x))
  
  results <- mRMR.classic(
    data = dd,
    target_indices = 1,
    feature_count = feature_count)
  
  selected <- featureNames(results)
  
  # remove target column if mRMRe returns it
  selected <- selected[selected != "outcome"]
  
  tibble(
    outcome = outcome_name,
    rank = seq_along(selected),
    feature = selected)}

# ------------------------------------------------------------
# Run MRMR for each eGFR outcome
# ------------------------------------------------------------

egfr_feature_results <- map_dfr(egfr_outcomes, function(outcome_name) {
  run_mrmr(
    outcome_vector = df_transplanted[[outcome_name]],
    features = features_numeric,
    outcome_name = outcome_name,
    top_n = 100)})

# ------------------------------------------------------------
# Run MRMR for DGF
# ------------------------------------------------------------

dgf_feature_results <- run_mrmr(
  outcome_vector = df_transplanted[[dgf_outcome]],
  features = features_numeric,
  outcome_name = dgf_outcome,
  top_n = 100)

# ------------------------------------------------------------
# Save results
# ------------------------------------------------------------

write.csv(
  egfr_feature_results,
  "top100_features_eGFR_all_outcomes.csv",
  row.names = FALSE)

write.csv(
  dgf_feature_results,
  "top100_features_DGF.csv",
  row.names = FALSE)

# Save one combined file
combined_feature_results <- bind_rows(
  egfr_feature_results,
  dgf_feature_results)

write.csv(combined_feature_results,
  "top100_features_all_outcomes.csv",
  row.names = FALSE)

cat("Done. Saved MRMR top feature files.\n")