# ============================================================
# CLEAN MRMR FEATURE SELECTION
# ============================================================

library(tidyverse)
library(caret)
library(mRMRe)

# LOAD DATA

df <- Reduced_Features_NO_COLOR

# transplanted only for outcomes
df_transplanted <- df %>%
  filter(Allocated_Discarded == 1)


# OUTCOMES

egfr_outcomes <- c(
  "eGFR_CKD_EPI_1M",
  "eGFR_CKD_EPI_3M",
  "eGFR_CKD_EPI_6M",
  "eGFR_CKD_EPI_12M",
  "eGFR_CKD_EPI_3years",
  "eGFR_CKD_EPI_5years")

dgf_outcome <- "Delayed_Graft_Function"


# EXACT COLUMNS TO REMOVE

cols_to_exclude <- c(
  "Allocated_Discarded",
  "Slide_number",
  "Recipient_code",
  "RECIPIENT_DATA",
  "Transplant_date_OR_Harvest_date",
  
  # outcomes
  egfr_outcomes,
  dgf_outcome,
  
  "Cr1M","Cr3M","Cr6M","Cr12M",
  "Cr3years","Cr5years","Cr10years",
  "Cr_at_last_Followup",
  
  # admin
  "Last_follow_up_date",
  "Death_date",
  "Graft_loss_date")


# MORE PATTERN FILTERS 

bad_name_patterns <- c(
  "^Recipient",
  "recipient",
  "RECIPIENT",
  
  "^eGFR",
  "eGFR_",
  
  "^Cr",
  "Creatinine",
  
  "Delayed_Graft_Function",
  "^DGF$",
  
  "Graft",
  "Death",
  "Dialysis",
  "Follow",
  "Transplant",
  "Return_to",
  
  "Recipient_",
  "Slide_number",
  "Allocated_Discarded")

bad_regex <- paste(bad_name_patterns, collapse = "|")


# BUILD CLEAN FEATURE MATRIX


all_bad_cols <- names(df_transplanted)[
  names(df_transplanted) %in% cols_to_exclude |
    stringr::str_detect(names(df_transplanted), bad_regex)]

features_numeric <- df_transplanted %>%
  select(-any_of(all_bad_cols)) %>%
  select(where(is.numeric))

cat("Numeric candidate features after leakage removal:",
    ncol(features_numeric), "\n")

# REMOVE LOW INFORMATION FEATURES

nzv <- nearZeroVar(features_numeric)

if(length(nzv) > 0){
  features_numeric <- features_numeric[, -nzv, drop = FALSE]}


# REMOVE HIGH CORRELATION

cor_matrix <- cor(features_numeric, use = "pairwise.complete.obs")
cor_matrix[is.na(cor_matrix)] <- 0

high_cor <- findCorrelation(cor_matrix, cutoff = 0.80)

if(length(high_cor) > 0){
  features_numeric <- features_numeric[, -high_cor, drop = FALSE]}

cat("Features remaining after cleanup:",
    ncol(features_numeric), "\n")

# MRMR FUNCTION


run_mrmr <- function(outcome_vector,
                     features,
                     outcome_name,
                     top_n = 100){
  
  y <- suppressWarnings(as.numeric(outcome_vector))
  
  keep <- !is.na(y)
  
  y <- y[keep]
  x <- features[keep, , drop = FALSE]
  
  # force numeric
  x <- x %>% mutate(across(everything(),
                  ~ suppressWarnings(as.numeric(.))))
  
  # remove all-NA cols
  x <- x[, colSums(!is.na(x)) > 0, drop = FALSE]
  
  # complete rows only
  keep2 <- complete.cases(x) & !is.na(y)
  
  y <- y[keep2]
  x <- x[keep2, , drop = FALSE]
  
  if(length(y) < 20 || ncol(x) < 5){
    warning(paste("Skipping", outcome_name))
    return(NULL)}
  
  feature_count <- min(top_n + 30, ncol(x))
  
  mrmr_df <- data.frame(
    outcome = y,
    x,
    check.names = FALSE)
  
  mrmr_df[] <- lapply(
    mrmr_df,
    function(z) suppressWarnings(as.numeric(z)))
  
  dd <- mRMR.data(data = mrmr_df)
  
  results <- mRMR.classic(
    data = dd,
    target_indices = 1,
    feature_count = feature_count)
  
  selected <- featureNames(results)
  
  selected <- selected[selected != "outcome"]
  
  selected <- selected[
    !stringr::str_detect(selected, bad_regex)]
  
  selected <- selected[
    !(selected %in% cols_to_exclude)]
  
  selected <- unique(selected)
  
  selected <- head(selected, top_n)
  
  tibble(
    outcome = outcome_name,
    rank = seq_along(selected),
    feature = selected)}


# RUN EGFR MRMR

egfr_feature_results <- map_dfr(
  egfr_outcomes,
  function(outcome_name){
    
    run_mrmr(
      outcome_vector = df_transplanted[[outcome_name]],
      features = features_numeric,
      outcome_name = outcome_name,
      top_n = 100)})

# RUN DGF MRMR

dgf_feature_results <- run_mrmr(
  outcome_vector = df_transplanted[[dgf_outcome]],
  features = features_numeric,
  outcome_name = dgf_outcome,
  top_n = 100)

# CHECK COUNTS

cat("\nCounts by eGFR outcome:\n")
print(egfr_feature_results %>% count(outcome))

cat("\nCounts for DGF:\n")
print(dgf_feature_results %>% count(outcome))

# ------------------------------------------------------------
# SAVE FILES
# ------------------------------------------------------------

write.csv(
  egfr_feature_results,
  "top100_features_eGFR_all_outcomes.csv",
  row.names = FALSE)

write.csv(
  dgf_feature_results,
  "top100_features_DGF.csv",
  row.names = FALSE)

combined_feature_results <- bind_rows(
  egfr_feature_results,
  dgf_feature_results)

write.csv(
  combined_feature_results,
  "top100_features_all_outcomes.csv",
  row.names = FALSE)

cat("\nDone. Clean MRMR files saved.\n")
