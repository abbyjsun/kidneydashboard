# app.R
# Kidney Transplant Outcomes Dashboard
# Abby Sun

library(shiny)
library(tidyverse)
library(pROC)
library(glmnet)
library(ranger)
library(xgboost)

# ============================================================
# CONFIG
# ============================================================

eGFR_outcomes <- c(
  "eGFR_CKD_EPI_1M",
  "eGFR_CKD_EPI_3M",
  "eGFR_CKD_EPI_6M",
  "eGFR_CKD_EPI_12M",
  "eGFR_CKD_EPI_3years",
  "eGFR_CKD_EPI_5years")

predictors_full <- c(
  "Donor_age",
  "Donor_Sex_0_Male_1_Female",
  "Donor_BMI",
  "Donor_Hypertension_History",
  "Donor_Diabetes_History",
  "Donor_final_creatinine",
  "KDPI_2024")

DGF_col <- "Delayed_Graft_Function"
dgf_predictors <- predictors_full
kdpi_predictors <- c("KDPI_2024")

stage_levels <- c(
  "G1 (≥90)", "G2 (60–89)", "G3a (45–59)",
  "G3b (30–44)", "G4 (15–29)", "G5 (<15)")

DATA_PATH <- "Reduced_Features_NO_COLOR 2.csv"

# MRMR feature-ranking files generated from your MRMR script
MRMR_EGFR_PATH <- "top100_features_eGFR_all_outcomes.csv"
MRMR_DGF_PATH <- "top100_features_DGF.csv"

ml_model_names <- c(
  "Lasso",
  "Ridge",
  "Elastic Net",
  "Random Forest",
  "XGBoost")

ml_predictor_set_names <- c(
  "Donor clinical/demographic only",
  "Joint MRMR-selected donor clinical + image features")

# ============================================================
# HELPER FUNCTIONS
# ============================================================

ckd_stage <- function(e) {
  dplyr::case_when(
    is.na(e) ~ NA_character_,
    e >= 90 ~ "G1 (≥90)",
    e >= 60 ~ "G2 (60–89)",
    e >= 45 ~ "G3a (45–59)",
    e >= 30 ~ "G3b (30–44)",
    e >= 15 ~ "G4 (15–29)",
    TRUE ~ "G5 (<15)")}

rmse_from_mse <- function(mse) sqrt(mse)

safe_as_numeric <- function(x) suppressWarnings(as.numeric(x))

safe_binary_01 <- function(x) {
  if (is.factor(x)) x <- as.character(x)
  
  if (is.character(x)) {
    x <- tolower(trimws(x))
    x <- ifelse(
      x %in% c("yes", "y", "true", "1"), 1,
      ifelse(x %in% c("no", "n", "false", "0"), 0, NA))}
  
  y <- suppressWarnings(as.numeric(x))
  y[!(y %in% c(0, 1))] <- NA_real_
  y}

predict_with_intervals <- function(fit, newdata_row) {
  ci <- as.data.frame(predict(fit, newdata = newdata_row, interval = "confidence"))
  pi <- as.data.frame(predict(fit, newdata = newdata_row, interval = "prediction"))
  
  tibble(
    Predicted_eGFR = as.numeric(ci$fit),
    Pred_CI_L = as.numeric(ci$lwr),
    Pred_CI_U = as.numeric(ci$upr),
    Pred_PI_L = as.numeric(pi$lwr),
    Pred_PI_U = as.numeric(pi$upr))}

pick_id_col <- function(d) {
  candidates <- c("Slide_number", "Donor_ID", "DonorID", "Donor_id", "Recipient_code", "ID", "Id", "id")
  hit <- candidates[candidates %in% names(d)]
  if (length(hit) > 0) return(hit[1])
  
  for (nm in names(d)) {
    if (!all(is.na(d[[nm]]))) return(nm)}
  names(d)[1]}

calc_auc_safe <- function(truth, score) {
  keep <- !is.na(truth) & !is.na(score)
  truth <- truth[keep]
  score <- score[keep]
  
  if (length(truth) == 0 || length(unique(truth)) < 2) return(NA_real_)
  as.numeric(pROC::auc(truth, score))}

roc_df_safe <- function(truth, score) {
  keep <- !is.na(truth) & !is.na(score)
  truth <- truth[keep]
  score <- score[keep]
  
  if (length(truth) == 0 || length(unique(truth)) < 2) return(NULL)
  
  roc_obj <- pROC::roc(response = truth, predictor = score, quiet = TRUE)
  
  tibble(
    specificity = roc_obj$specificities,
    sensitivity = roc_obj$sensitivities,
    fpr = 1 - specificity)}

top_coef_table <- function(fit, top_n = 5) {
  co <- coef(fit)
  co <- co[names(co) != "(Intercept)"]
  
  tibble(
    Predictor = names(co),
    Coefficient = as.numeric(co),
    Abs_Coefficient = abs(as.numeric(co))) %>%
    arrange(desc(Abs_Coefficient)) %>%
    slice_head(n = top_n) %>%
    select(Predictor, Coefficient)}

# ============================================================
# MACHINE LEARNING HELPER FUNCTIONS
# ============================================================

is_leakage_or_recipient_feature <- function(x) {
  bad_patterns <- c(
    "^Recipient", "recipient", "RECIPIENT",
    "eGFR", "Cr[0-9]", "Creatinine_", "Delayed_Graft_Function", "DGF",
    "Graft", "Death", "Dialysis", "Follow", "Transplant_date", "Return_to",
    "Allocated_Discarded", "Slide_number")
  str_detect(x, paste(bad_patterns, collapse = "|"))}

read_mrmr_features <- function(path, outcome_name, data_names, top_n = 100) {
  if (!file.exists(path)) return(character(0))
  
  mrmr <- read.csv(path, stringsAsFactors = FALSE)
  if (!("feature" %in% names(mrmr))) return(character(0))
  
  if ("outcome" %in% names(mrmr)) {
    mrmr <- mrmr %>% filter(.data$outcome == outcome_name)}
  
  if ("rank" %in% names(mrmr)) {
    mrmr <- mrmr %>% arrange(.data$rank)}
  
  feats <- mrmr$feature
  feats <- feats[!is.na(feats)]
  feats <- feats[!(feats %in% c("outcome", outcome_name, DGF_col, eGFR_outcomes))]
  feats <- feats[!is_leakage_or_recipient_feature(feats)]
  feats <- feats[feats %in% data_names]
  unique(head(feats, top_n))}

mrmr_features_for_outcome <- function(outcome_name, data_names, top_n = 100) {
  if (outcome_name == DGF_col) {
    read_mrmr_features(MRMR_DGF_PATH, outcome_name, data_names, top_n)} 
  else {read_mrmr_features(MRMR_EGFR_PATH, outcome_name, data_names, top_n)}}

ml_predictor_sets_for_outcome <- function(data_names, outcome_name, top_n = 100) {
  clinical <- predictors_full[predictors_full %in% data_names]
  mrmr_feats <- mrmr_features_for_outcome(outcome_name, data_names, top_n)
  
  list(
    "Donor clinical/demographic only" = clinical,
    "Joint MRMR-selected donor clinical + image features" = mrmr_feats)}

prepare_ml_data <- function(dat, outcome, predictors, outcome_type) {
  needed <- unique(c(outcome, predictors))
  needed <- needed[needed %in% names(dat)]
  if (!(outcome %in% needed)) return(NULL)
  
  sub <- dat %>% select(all_of(needed))
  y <- sub[[outcome]]
  
  if (outcome_type == "Binary") {
    y <- safe_binary_01(y)} 
  else {
    y <- safe_as_numeric(y)}
  
  x_df <- sub %>% select(-all_of(outcome))
  x_df <- x_df %>% mutate(across(everything(), ~ suppressWarnings(as.numeric(.))))
  x_df <- x_df[, colSums(!is.na(x_df)) > 0, drop = FALSE]
  
  keep <- !is.na(y) & complete.cases(x_df)
  y <- y[keep]
  x_df <- x_df[keep, , drop = FALSE]
  
  if (length(y) < 20 || ncol(x_df) < 1) return(NULL)
  if (outcome_type == "Binary" && length(unique(y)) < 2) return(NULL)
  
  x <- model.matrix(~ ., data = x_df)[, -1, drop = FALSE]
  if (ncol(x) < 1) return(NULL)
  
  list(y = y, x_df = x_df, x = x, predictors = names(x_df))}

align_ml_matrix <- function(new_x_df, train_x_cols) {
  new_x <- model.matrix(~ ., data = new_x_df)[, -1, drop = FALSE]
  colnames(new_x) <- make.names(colnames(new_x), unique = TRUE)  # ADD THIS
  
  missing_cols <- setdiff(train_x_cols, colnames(new_x))
  if (length(missing_cols) > 0) {
    add <- matrix(0, nrow = nrow(new_x), ncol = length(missing_cols))
    colnames(add) <- missing_cols
    new_x <- cbind(new_x, add)}
  extra_cols <- setdiff(colnames(new_x), train_x_cols)
  if (length(extra_cols) > 0) {
    new_x <- new_x[, !(colnames(new_x) %in% extra_cols), drop = FALSE]}
  new_x[, train_x_cols, drop = FALSE]}

fit_ml_core <- function(x, x_df, y, model_name, outcome_type) {
  family_type <- ifelse(outcome_type == "Binary", "binomial", "gaussian")
  
  tryCatch({
    if (model_name %in% c("Lasso", "Ridge", "Elastic Net")) {
      alpha_val <- switch(
        model_name,
        "Lasso" = 1,
        "Ridge" = 0,
        "Elastic Net" = 0.5)
      
      cvfit <- glmnet::cv.glmnet(
        x = x,
        y = y,
        family = family_type,
        alpha = alpha_val,
        nfolds = min(5, nrow(x)),
        standardize = TRUE)
      
      list(type = "glmnet", model = cvfit, x_cols = colnames(x))} else if (model_name == "Random Forest") {
      train_df <- data.frame(y = y, x_df, check.names = FALSE)
      
      if (outcome_type == "Binary") {
        train_df$y <- factor(train_df$y, levels = c(0, 1))
        rf <- ranger::ranger(
          y ~ .,
          data = train_df,
          probability = TRUE,
          num.trees = 500
        )} else {
        rf <- ranger::ranger(
          y ~ .,
          data = train_df,
          num.trees = 500)}
      
      list(type = "ranger", model = rf, predictors = names(x_df), outcome_type = outcome_type)} else if (model_name == "XGBoost") {
        objective <- ifelse(outcome_type == "Binary", "reg:logistic", "reg:squarederror")
        y_xgb <- as.numeric(y)
        dtrain <- xgboost::xgb.DMatrix(data = x)
        xgboost::setinfo(dtrain, "label", y_xgb)
        xgb <- xgboost::xgb.train(
          params = list(
            objective        = objective,
            max_depth        = 3,
            learning_rate    = 0.05,
            subsample        = 0.8,
            colsample_bytree = 0.8
          ),
          data    = dtrain,
          nrounds = 100)
        list(type = "xgboost", model = xgb, x_cols = colnames(x))} else {
      NULL}}, error = function(e) {
        message("fit_ml_core ERROR [", model_name, "]: ", conditionMessage(e))
        NULL
      })}

predict_ml_core <- function(fit_obj, x, x_df, outcome_type) {
  tryCatch({
    if (fit_obj$type == "glmnet") {
      as.numeric(predict(
        fit_obj$model,
        newx = x,
        s = "lambda.min",
        type = ifelse(outcome_type == "Binary", "response", "link")
      ))
    } else if (fit_obj$type == "ranger") {
      p <- predict(fit_obj$model, data = x_df)$predictions
      if (outcome_type == "Binary") as.numeric(p[, "1"]) else as.numeric(p)
    } else if (fit_obj$type == "xgboost") {
      as.numeric(predict(fit_obj$model, newdata = x))
    } else {
      rep(NA_real_, nrow(x_df))}}, error = function(e) rep(NA_real_, nrow(x_df)))}

make_cv_folds <- function(y, k = 5, seed = 123) {
  set.seed(seed)
  n <- length(y)
  k <- min(k, n)
  
  if (length(unique(y)) == 2) {
    folds <- rep(NA_integer_, n)
    for (cls in unique(y)) {
      idx <- which(y == cls)
      folds[idx] <- sample(rep(seq_len(k), length.out = length(idx)))
    }
    folds} else {
    sample(rep(seq_len(k), length.out = n))}}

fit_ml_model <- function(dat, outcome, predictors, predictor_set, model_name, outcome_type) {
  ml_dat <- prepare_ml_data(dat, outcome, predictors, outcome_type)
  if (is.null(ml_dat)) return(NULL)
  
  y <- ml_dat$y
  x <- ml_dat$x
  x_df <- ml_dat$x_df
  
  n <- length(y)
  k <- ifelse(n < 40, 3, 5)
  folds <- make_cv_folds(y, k = k, seed = 123)
  
  cv_pred <- rep(NA_real_, n)
  
  for (fold in sort(unique(folds))) {
    train_idx <- which(folds != fold)
    test_idx <- which(folds == fold)
    
    # Binary outcomes need both outcome classes in training.
    if (outcome_type == "Binary" && length(unique(y[train_idx])) < 2) next
    
    fit_fold <- fit_ml_core(
      x = x[train_idx, , drop = FALSE],
      x_df = x_df[train_idx, , drop = FALSE],
      y = y[train_idx],
      model_name = model_name,
      outcome_type = outcome_type)
    
    if (is.null(fit_fold)) next
    
    # Align XGBoost/glmnet test matrices to training columns if needed.
    if (fit_fold$type %in% c("glmnet", "xgboost")) {
      x_test <- x[test_idx, , drop = FALSE]
      missing_cols <- setdiff(fit_fold$x_cols, colnames(x_test))
      if (length(missing_cols) > 0) {
        add <- matrix(0, nrow = nrow(x_test), ncol = length(missing_cols))
        colnames(add) <- missing_cols
        x_test <- cbind(x_test, add)
      }
      extra_cols <- setdiff(colnames(x_test), fit_fold$x_cols)
      if (length(extra_cols) > 0) {
        x_test <- x_test[, !(colnames(x_test) %in% extra_cols), drop = FALSE]}
      x_test <- x_test[, fit_fold$x_cols, drop = FALSE]} else {
      x_test <- x[test_idx, , drop = FALSE]}
    
    cv_pred[test_idx] <- predict_ml_core(
      fit_obj = fit_fold,
      x = x_test,
      x_df = x_df[test_idx, , drop = FALSE],
      outcome_type = outcome_type)}
  
  keep_pred <- !is.na(cv_pred)
  if (sum(keep_pred) < 10) return(NULL)
  
  if (outcome_type == "Continuous") {
    mse <- mean((y[keep_pred] - cv_pred[keep_pred])^2, na.rm = TRUE)
    auc <- NA_real_
    score <- mse} else {
    mse <- NA_real_
    auc <- calc_auc_safe(y[keep_pred], cv_pred[keep_pred])
    score <- auc}
  
  # Final refit on all available complete TX rows for prediction/recommendation
  final_fit <- fit_ml_core(
    x = x,
    x_df = x_df,
    y = y,
    model_name = model_name,
    outcome_type = outcome_type)
  
  if (is.null(final_fit)) return(NULL)
  
  list(
    perf = tibble(
      Outcome = outcome,
      Outcome_Type = outcome_type,
      Predictor_Set = predictor_set,
      Model = model_name,
      N = length(y),
      Num_Predictors = length(unique(ml_dat$predictors)),
      MSE = mse,
      RMSE = ifelse(is.na(mse), NA_real_, sqrt(mse)),
      AUC = auc,
      Selection_Metric = ifelse(outcome_type == "Continuous", "Lowest CV MSE", "Highest CV AUC"),
      Score_For_Selection = score),
    fit = list(
      outcome = outcome,
      outcome_type = outcome_type,
      predictor_set = predictor_set,
      model_name = model_name,
      predictors = ml_dat$predictors,
      fit_obj = final_fit))}

predict_ml_fit <- function(fit_record, newdata) {
  predictors <- fit_record$predictors
  fit_obj <- fit_record$fit_obj
  outcome_type <- fit_record$outcome_type
  
  if (!all(predictors %in% names(newdata))) {
    stop("Missing predictors in newdata: ",
         paste(setdiff(predictors, names(newdata)), collapse = ", "))}
  
  new_x_df <- newdata %>%
    select(all_of(predictors)) %>%
    mutate(across(everything(), ~ suppressWarnings(as.numeric(.))))
  
  if (nrow(new_x_df) != 1) {
    stop("new_x_df does not have exactly 1 row.")}
  
  if (!complete.cases(new_x_df)) {
    stop("Predictors became NA after numeric conversion: ",
         paste(names(new_x_df)[is.na(new_x_df[1, ])], collapse = ", "))}
  
  if (fit_obj$type == "glmnet") {
    new_x <- align_ml_matrix(new_x_df, fit_obj$x_cols)
    return(as.numeric(predict(
      fit_obj$model,
      newx = new_x,
      s = "lambda.min",
      type = ifelse(outcome_type == "Binary", "response", "link")
    )))}
  
  if (fit_obj$type == "ranger") {
    p <- predict(fit_obj$model, data = new_x_df)$predictions
    if (outcome_type == "Binary") return(as.numeric(p[, "1"]))
    return(as.numeric(p))}
  
  if (fit_obj$type == "xgboost") {
    new_x <- align_ml_matrix(new_x_df, fit_obj$x_cols)
    return(as.numeric(predict(fit_obj$model, newdata = new_x)))}
  
  stop("Unknown model type.")}

# ============================================================
# UI
# ============================================================

ui <- fluidPage(
  titlePanel("Kidney Transplant Outcomes Dashboard"),
  
  sidebarLayout(
    sidebarPanel(
      helpText("Main TX models use all available complete transplanted subjects. Non-TX prediction-only models use an 80/20 TX training split."),
      
      radioButtons(
        "model_mode", "eGFR model type",
        choices = c(
          "KDPI-only" = "kdpi_only",
          "Multivariable donor predictors" = "multi"),
        selected = "multi"),
      
      numericInput("seed", "Random seed for Non-TX 80/20 training split", value = 123, min = 1, step = 1),
      
      sliderInput(
        "train_frac",
        "Training fraction for Non-TX prediction model",
        min = 0.5, max = 0.9, value = 0.8, step = 0.05),
      
      hr(),
      
      selectInput(
        "cohort_view",
        "CKD Stage Plot Cohort",
        choices = c("Transplant (Observed)", "Non-Transplant (Prediction only)"),
        selected = "Transplant (Observed)"),
      
      selectInput(
        "outcome_view",
        "Outcome for CKD Stage Plot / Pred vs True",
        choices = c("All outcomes (faceted)", eGFR_outcomes),
        selected = "All outcomes (faceted)"),
      
      hr(),
      h4("DGF logistic regression"),
      sliderInput(
        "top_n_dgf",
        "Top DGF predictors to display",
        min = 3, max = 10, value = 5, step = 1),
      
      hr(),
      h4("Machine learning models"),
      selectInput(
        "ml_selected_outcome",
        "ML outcome to display",
        choices = c(eGFR_outcomes, DGF_col),
        selected = DGF_col),
      selectInput(
        "ml_selected_predictor_set",
        "ML predictor set",
        choices = ml_predictor_set_names,
        selected = "Joint MRMR-selected donor clinical + image features"),
      actionButton("run_ml", "Run / Refresh ML models"),
      
      hr(),
      h4("Kidney recommendation"),
      selectInput(
        "rec_egfr_outcome",
        "eGFR outcome for recommendation",
        choices = eGFR_outcomes,
        selected = "eGFR_CKD_EPI_12M"),
      numericInput(
        "rec_egfr_cutoff",
        "Minimum acceptable predicted eGFR",
        value = 45, min = 0, max = 120, step = 5),
      numericInput(
        "rec_dgf_cutoff",
        "Maximum acceptable predicted DGF probability",
        value = 0.50, min = 0, max = 1, step = 0.05),
      
      hr(),
      h4("Individual recipient (TX cohort)"),
      uiOutput("recipient_ui"),
      
      hr(),
      h4("Individual Non-TX donor (prediction-only cohort)"),
      uiOutput("donor_ui")),
    
    mainPanel(
      tabsetPanel(
        tabPanel(
          "eGFR Performance",
          br(),
          plotOutput("mse_plot", height = "320px"),
          br(),
          tableOutput("perf_table"),
          br(),
          h4("Outcome availability (TX cohort only)"),
          tableOutput("availability_table"),
          br(),
          h4("Predictor coefficients from all-subject TX model for selected outcome"),
          tags$p(
            style = "color:#555;",
            "Select a single eGFR outcome in the sidebar to populate this. KDPI-only models have one coefficient."),
          tableOutput("coef_table")),
        
        tabPanel(
          "Pred vs True",
          br(),
          uiOutput("pred_true_msg"),
          plotOutput("pred_true_plot", height = "520px")),
        
        tabPanel(
          "CKD Stage Distributions",
          br(),
          plotOutput("ckd_plot", height = "550px"),
          br(),
          tags$ul(
            tags$li("Transplant (Observed) = fitted predictions from all available transplanted subjects."),
            tags$li("Non-Transplant cohort = prediction-only rows without observed recipient eGFR outcomes."),
            tags$li("Non-TX predictions are generated from models trained on the 80% TX training split only."))),
        
        tabPanel(
          "Individual Recipient",
          br(),
          h4("Predicted vs true eGFR + CKD stage using all-subject TX model"),
          tableOutput("recipient_compare_table")),
        
        tabPanel(
          "Individual Non-TX Donor Predictions",
          br(),
          tags$p(
            style = "color:#b00;",
            "Prediction-only: Non-TX rows do not have observed recipient eGFR outcomes. Predictions use an 80/20-trained TX model."),
          tableOutput("donor_pred_table")),
        
        tabPanel(
          "DGF Performance",
          br(),
          h4("Delayed graft function model performance using all available TX subjects"),
          tableOutput("dgf_perf_table"),
          br(),
          h4("Top predictors by coefficient magnitude"),
          tableOutput("dgf_coef_table")),
        
        tabPanel(
          "DGF ROC",
          br(),
          plotOutput("dgf_roc_plot", height = "520px")),
        
        tabPanel(
          "DGF Individual Predictions",
          br(),
          h4("Selected recipient predicted probability of DGF using all-subject TX model"),
          tableOutput("recipient_dgf_table"),
          br(),
          h4("Selected Non-TX donor predicted probability of DGF"),
          tags$p(
            style = "color:#b00;",
            "Prediction-only: Non-TX rows do not have observed recipient DGF outcomes. Predictions use an 80/20-trained TX logistic model."),
          tableOutput("donor_dgf_table")),
        
        tabPanel(
          "ML Model Performance",
          br(),
          
          tags$div(
            style = "color:#555; font-size:14px;",
            
            tags$p("Click Run / Refresh ML models after changing MRMR settings."),
            
            tags$ul(
              tags$li("Performance uses internal cross-validation within transplanted subjects."),
              tags$li("Not hold-out 20% test-set metrics."),
              tags$li("Best eGFR model = lowest CV MSE."),
              tags$li("Best DGF model = highest CV AUC."),
              tags$li("Final models are refit on all complete TX rows."))),
          
          plotOutput("ml_perf_plot", height = "420px"),
          br(),
          tableOutput("ml_perf_table")),
        
        tabPanel(
          "Best ML Model",
          br(),
          uiOutput("ml_best_msg"),
          tableOutput("ml_best_table")),
        
        tabPanel(
          "Kidney Recommendation",
          br(),
          tags$p(
            style = "color:#b00;",
            "Exploratory prediction-only recommendation for the selected Non-TX donor. This is not clinical guidance. Cutoffs can be edited on the side bar."),
          tableOutput("kidney_recommendation_table"))))))

# ============================================================
# SERVER
# ============================================================

server <- function(input, output, session) {
  
  observe({
    if (!file.exists(DATA_PATH)) {
      showNotification(
        paste0(
          "ERROR: Data file not found: '", DATA_PATH, "'. ",
          "Make sure the CSV is in the same folder as app.R and deployed with the app."
        ),
        type = "error",
        duration = NULL)}})
  
  df_raw <- reactive({
    validate(
      need(file.exists(DATA_PATH), paste0("Missing file: ", DATA_PATH, ". Upload/deploy it with the app.")))
    
    d <- read.csv(
      DATA_PATH,
      stringsAsFactors = FALSE,
      na.strings = c("", " ", "NA"))
    
    names(d) <- trimws(names(d))
    d})
  
  df <- reactive({
    d <- df_raw()
    
    numeric_predictors <- unique(c(
      predictors_full,
      "Donor_height",
      "Donor_weight_Kg",
      "Donor_cause_of_death_code",
      "Donor_eGFR_CKD_EPI_final_creatinine",
      kdpi_predictors))
    
    for (x in numeric_predictors) {
      if (x %in% names(d)) d[[x]] <- safe_as_numeric(d[[x]])}
    
    for (y in eGFR_outcomes) {
      if (y %in% names(d)) d[[y]] <- safe_as_numeric(d[[y]])}
    
    if (DGF_col %in% names(d)) {
      d[[DGF_col]] <- safe_binary_01(d[[DGF_col]])}
    
    d})
  
  # ------------------------------------------------------------
  # Cohorts
  # ------------------------------------------------------------
  
  cohorts <- reactive({
    d <- df()
    
    tx_df <- d %>% filter(if_any(all_of(eGFR_outcomes), ~ !is.na(.x)))
    non_tx_df <- d %>% filter(if_all(all_of(eGFR_outcomes), ~ is.na(.x)))
    
    list(tx = tx_df, non_tx = non_tx_df)})
  
  predictors_active <- reactive({
    if (input$model_mode == "kdpi_only") kdpi_predictors else predictors_full})
  
  tx_split_for_nontx <- reactive({
    set.seed(input$seed)
    
    tx <- cohorts()$tx
    n <- nrow(tx)
    if (n < 10) return(NULL)
    
    train_n <- floor(input$train_frac * n)
    idx <- sample(seq_len(n), size = train_n)
    
    list(
      train = tx[idx, , drop = FALSE],
      test = tx[-idx, , drop = FALSE])})
  
  # ------------------------------------------------------------
  # eGFR all-subject TX models
  # ------------------------------------------------------------
  
  egfr_models_all_tx <- reactive({
    tx <- cohorts()$tx
    preds <- predictors_active()
    mods <- list()
    
    for (outcome in eGFR_outcomes) {
      needed <- c(outcome, preds)
      sub <- tx %>% select(all_of(needed)) %>% drop_na()
      
      if (nrow(sub) < max(8, length(preds) + 5)) next
      
      f <- as.formula(paste(outcome, "~", paste(preds, collapse = " + ")))
      mods[[outcome]] <- list(
        formula = f,
        fit = lm(f, data = sub),
        data = sub)}
    
    mods})
  
  egfr_models_train_split <- reactive({
    sp <- tx_split_for_nontx()
    if (is.null(sp)) return(list())
    
    train_df <- sp$train
    preds <- predictors_active()
    mods <- list()
    
    for (outcome in eGFR_outcomes) {
      needed <- c(outcome, preds)
      sub <- train_df %>% select(all_of(needed)) %>% drop_na()
      
      if (nrow(sub) < max(8, length(preds) + 5)) next
      
      f <- as.formula(paste(outcome, "~", paste(preds, collapse = " + ")))
      mods[[outcome]] <- list(
        formula = f,
        fit = lm(f, data = sub),
        data = sub)}
    
    mods})
  
  egfr_performance <- reactive({
    mods <- egfr_models_all_tx()
    if (length(mods) == 0) return(NULL)
    
    rows <- lapply(names(mods), function(outcome) {
      obj <- mods[[outcome]]
      dat <- obj$data
      pred <- as.numeric(predict(obj$fit, newdata = dat))
      mse <- mean((dat[[outcome]] - pred)^2)
      
      tibble(
        Outcome = outcome,
        N = nrow(dat),
        MSE = mse,
        RMSE = rmse_from_mse(mse))})
    
    bind_rows(rows) %>% arrange(match(Outcome, eGFR_outcomes))})
  
  availability_table <- reactive({
    tx <- cohorts()$tx
    preds <- predictors_active()
    
    lapply(eGFR_outcomes, function(outcome) {
      n_true <- tx %>% filter(!is.na(.data[[outcome]])) %>% nrow()
      n_complete <- tx %>% select(all_of(c(outcome, preds))) %>% drop_na() %>% nrow()
      
      tibble(
        Outcome = outcome,
        TX_true_available = n_true,
        TX_complete_preds_plus_true = n_complete)}) %>% bind_rows()})
  
  coef_table <- reactive({
    sel <- input$outcome_view
    if (!(sel %in% eGFR_outcomes)) return(NULL)
    
    mods <- egfr_models_all_tx()
    if (!(sel %in% names(mods))) return(NULL)
    
    co <- coef(mods[[sel]]$fit)
    co <- co[names(co) != "(Intercept)"]
    
    tibble(
      Predictor = names(co),
      Coefficient = as.numeric(co),
      Abs = abs(as.numeric(co))) %>%
      arrange(desc(Abs)) %>%
      select(Predictor, Coefficient)})
  
  output$perf_table <- renderTable({
    egfr_performance()}, digits = 3)
  
  output$availability_table <- renderTable({
    availability_table()}, digits = 0)
  
  output$coef_table <- renderTable({
    coef_table()}, digits = 4)
  
  output$mse_plot <- renderPlot({
    perf <- egfr_performance()
    req(perf)
    
    title_text <- if (input$model_mode == "kdpi_only") {
      "MSE from all available TX subjects (KDPI-only OLS)"} 
    else {"MSE from all available TX subjects (Multivariable OLS)"}
    
    ggplot(perf, aes(x = Outcome, y = MSE)) +
      geom_col() +
      labs(
        title = title_text,
        x = "Outcome",
        y = "MSE (eGFR²)") +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))})
  
  # ------------------------------------------------------------
  # CKD stage distributions
  # ------------------------------------------------------------
  
  pred_long <- reactive({
    tx <- cohorts()$tx
    non_tx <- cohorts()$non_tx
    preds <- predictors_active()
    
    if (input$cohort_view == "Transplant (Observed)") {
      mods <- egfr_models_all_tx()
      source_dat <- tx
      label <- "Transplant (Observed)"} else {
      mods <- egfr_models_train_split()
      source_dat <- non_tx
      label <- "Non-Transplant (Prediction only)"}
    
    if (length(mods) == 0) return(NULL)
    
    base <- source_dat %>% select(all_of(preds)) %>% drop_na()
    if (nrow(base) == 0) return(NULL)
    
    out_list <- list()
    
    for (outcome in names(mods)) {
      pred <- suppressWarnings(as.numeric(predict(mods[[outcome]]$fit, newdata = base)))
      
      out_list[[outcome]] <- tibble(
        Cohort = label,
        Outcome = outcome,
        Pred_eGFR = pred,
        Pred_CKD_Stage = factor(ckd_stage(pred), levels = stage_levels))}
    
    bind_rows(out_list) %>% filter(!is.na(Pred_CKD_Stage))})
  
  output$ckd_plot <- renderPlot({
    pdat <- pred_long()
    req(pdat)
    
    if (input$outcome_view == "All outcomes (faceted)") {
      ggplot(pdat, aes(x = Pred_CKD_Stage)) +
        geom_bar() +
        facet_wrap(~ Outcome, scales = "free_y") +
        labs(
          title = paste("Predicted CKD Stage Distribution –", input$cohort_view),
          x = "Predicted CKD Stage",
          y = "Count"
        ) +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
    } else {
      pdat2 <- pdat %>% filter(Outcome == input$outcome_view)
      
      ggplot(pdat2, aes(x = Pred_CKD_Stage)) +
        geom_bar() +
        labs(
          title = paste("Predicted CKD Stage –", input$cohort_view, "-", input$outcome_view),
          x = "Predicted CKD Stage",
          y = "Count") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))}})
  
  # ------------------------------------------------------------
  # Recipient selector + comparison table
  # ------------------------------------------------------------
  
  output$recipient_ui <- renderUI({
    tx <- cohorts()$tx
    if (!("Recipient_code" %in% names(tx))) {
      return(helpText("Recipient_code column not found in dataset."))}
    
    choices <- tx %>%
      mutate(Recipient_code = as.character(Recipient_code)) %>%
      distinct(Recipient_code) %>%
      arrange(Recipient_code) %>%
      pull(Recipient_code) %>%
      head(800)
    
    selectizeInput(
      "recipient_code",
      "Select Recipient_code",
      choices = choices,
      selected = if (length(choices) > 0) choices[1] else NULL,
      options = list(placeholder = "Type to search…"))})
  
  output$recipient_compare_table <- renderTable({
    tx <- cohorts()$tx
    req("Recipient_code" %in% names(tx), !is.null(input$recipient_code))
    
    rec <- tx %>%
      mutate(Recipient_code = as.character(Recipient_code)) %>%
      filter(Recipient_code == as.character(input$recipient_code)) %>%
      dplyr::slice(1)
    
    if (nrow(rec) != 1) return(NULL)
    
    mods <- egfr_models_all_tx()
    preds <- predictors_active()
    
    if (length(mods) == 0) {
      return(tibble(Message = "No eGFR models available with the current settings."))}
    
    rec_pred_base <- rec %>% select(all_of(preds)) %>% drop_na()
    if (nrow(rec_pred_base) != 1) {
      return(tibble(Message = "Selected recipient has missing predictors; cannot generate individualized prediction."))}
    
    rows <- lapply(eGFR_outcomes, function(outcome) {
      obj <- mods[[outcome]]
      if (is.null(obj)) return(NULL)
      
      true <- safe_as_numeric(rec[[outcome]])
      
      if (input$model_mode == "kdpi_only") {
        pred <- as.numeric(predict(obj$fit, newdata = rec_pred_base))
        
        tibble(
          Recipient_code = as.character(rec$Recipient_code),
          Outcome = outcome,
          Predicted_eGFR = pred,
          Predicted_CKD_Stage = ckd_stage(pred),
          True_eGFR = true,
          True_CKD_Stage = ckd_stage(true))} else {
        ints <- predict_with_intervals(obj$fit, rec_pred_base)
        
        tibble(
          Recipient_code = as.character(rec$Recipient_code),
          Outcome = outcome,
          Predicted_eGFR = ints$Predicted_eGFR,
          Pred_CI_L = ints$Pred_CI_L,
          Pred_CI_U = ints$Pred_CI_U,
          Pred_PI_L = ints$Pred_PI_L,
          Pred_PI_U = ints$Pred_PI_U,
          Predicted_CKD_Stage = ckd_stage(ints$Predicted_eGFR),
          True_eGFR = true,
          True_CKD_Stage = ckd_stage(true) )}})
    
    bind_rows(rows)}, digits = 2)
  
  # ------------------------------------------------------------
  # Non-TX donor selector + prediction table
  # ------------------------------------------------------------
  
  output$donor_ui <- renderUI({
    non_tx <- cohorts()$non_tx
    
    if (nrow(non_tx) == 0) {
      return(helpText("No Non-Transplant (prediction-only) rows found in dataset."))}
    
    id_col <- pick_id_col(non_tx)
    
    choices <- non_tx %>%
      distinct(.data[[id_col]]) %>%
      slice_head(n = 800) %>%
      pull()
    
    selectizeInput(
      "donor_id",
      paste0("Select Non-TX Donor ID (", id_col, ")"),
      choices = choices,
      selected = if (length(choices) > 0) choices[1] else NULL,
      options = list(placeholder = "Type to search…"))})
  
  output$donor_pred_table <- renderTable({
    non_tx <- cohorts()$non_tx
    req(nrow(non_tx) > 0, !is.null(input$donor_id))
    
    id_col <- pick_id_col(non_tx)
    rec <- non_tx %>% filter(.data[[id_col]] == input$donor_id) %>% dplyr::slice(1)
    if (nrow(rec) != 1) return(NULL)
    
    mods <- egfr_models_train_split()
    preds <- predictors_active()
    
    if (length(mods) == 0) {
      return(tibble(Message = "No eGFR models available with the current 80/20 training split."))}
    
    rec_pred_base <- rec %>% select(all_of(preds)) %>% drop_na()
    if (nrow(rec_pred_base) != 1) {
      return(tibble(Message = "Selected Non-TX donor has missing predictors; cannot generate individualized prediction."))}
    
    rows <- lapply(eGFR_outcomes, function(outcome) {
      obj <- mods[[outcome]]
      if (is.null(obj)) return(NULL)
      
      if (input$model_mode == "kdpi_only") {
        pred <- as.numeric(predict(obj$fit, newdata = rec_pred_base))
        
        tibble(
          NonTX_ID = as.character(rec[[id_col]]),
          Outcome = outcome,
          Predicted_eGFR = pred,
          Predicted_CKD_Stage = ckd_stage(pred),
          Ground_Truth_Available = "No")} 
      else {ints <- predict_with_intervals(obj$fit, rec_pred_base)
        
        tibble(
          NonTX_ID = as.character(rec[[id_col]]),
          Outcome = outcome,
          Predicted_eGFR = ints$Predicted_eGFR,
          Pred_CI_L = ints$Pred_CI_L,
          Pred_CI_U = ints$Pred_CI_U,
          Pred_PI_L = ints$Pred_PI_L,
          Pred_PI_U = ints$Pred_PI_U,
          Predicted_CKD_Stage = ckd_stage(ints$Predicted_eGFR),
          Ground_Truth_Available = "No")}})
    
    bind_rows(rows)}, digits = 2)
  
  # ------------------------------------------------------------
  # Predicted vs true plot
  # ------------------------------------------------------------
  
  output$pred_true_msg <- renderUI({
    if (input$outcome_view == "All outcomes (faceted)") {
      tags$p(
        style = "color:#b00;",
        "Select a single outcome, such as eGFR_CKD_EPI_6M, to view Predicted vs True.")} else {
      NULL }})
  
  output$pred_true_plot <- renderPlot({
    sel <- input$outcome_view
    req(sel != "All outcomes (faceted)")
    
    mods <- egfr_models_all_tx()
    preds <- predictors_active()
    
    obj <- mods[[sel]]
    req(!is.null(obj))
    
    dat <- obj$data
    pred <- as.numeric(predict(obj$fit, newdata = dat))
    dplt <- tibble(Pred = pred, True = dat[[sel]])
    
    cal_fit <- lm(True ~ Pred, data = dplt)
    r2 <- summary(cal_fit)$r.squared
    
    model_label <- if (input$model_mode == "kdpi_only") "KDPI-only" else "Multivariable OLS"
    
    ggplot(dplt, aes(x = Pred, y = True)) +
      geom_point(alpha = 0.7) +
      geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
      geom_smooth(method = "lm", se = FALSE) +
      labs(
        title = paste("Predicted vs True eGFR —", sel, "(", model_label, ")"),
        subtitle = paste0(
          "All available complete TX subjects. Dashed = y = x line. Solid = calibration. R² = ",
          round(r2, 3)),
        x = "Predicted eGFR",
        y = "True eGFR") +
      theme_bw()})
  
  # ------------------------------------------------------------
  # DGF logistic regression using all available TX subjects
  # ------------------------------------------------------------
  
  dgf_logit_data_all_tx <- reactive({
    tx <- cohorts()$tx
    req(DGF_col %in% names(tx))
    
    tx %>%
      select(all_of(c(DGF_col, dgf_predictors))) %>%
      drop_na()})
  
  dgf_logit_model_all_tx <- reactive({
    dat <- dgf_logit_data_all_tx()
    
    validate(
      need(nrow(dat) >= 15, "Not enough complete DGF cases to fit logistic regression."),
      need(length(unique(dat[[DGF_col]])) == 2, "DGF outcome must contain both 0 and 1."))
    
    f <- as.formula(paste(DGF_col, "~", paste(dgf_predictors, collapse = " + ")))
    glm(f, data = dat, family = binomial())})
  
  dgf_train_split_data <- reactive({
    set.seed(input$seed)
    
    dat <- dgf_logit_data_all_tx()
    if (nrow(dat) < 20) return(NULL)
    
    train_n <- floor(input$train_frac * nrow(dat))
    idx <- sample(seq_len(nrow(dat)), size = train_n)
    dat[idx, , drop = FALSE]})
  
  dgf_logit_model_train_split <- reactive({
    dat <- dgf_train_split_data()
    
    validate(
      need(!is.null(dat), "Not enough DGF cases for the 80/20 Non-TX training model."),
      need(nrow(dat) >= 15, "Not enough complete DGF cases to fit the 80/20 Non-TX training model."),
      need(length(unique(dat[[DGF_col]])) == 2, "80/20 DGF training data must contain both 0 and 1."))
    
    f <- as.formula(paste(DGF_col, "~", paste(dgf_predictors, collapse = " + ")))
    glm(f, data = dat, family = binomial())})
  
  dgf_kdpi_direct_data <- reactive({
    tx <- cohorts()$tx
    req(DGF_col %in% names(tx), "KDPI_2024" %in% names(tx))
    
    tx %>%
      select(all_of(c(DGF_col, "KDPI_2024"))) %>%
      drop_na()})
  
  dgf_perf <- reactive({
    fit_logit <- dgf_logit_model_all_tx()
    logit_dat <- dgf_logit_data_all_tx()
    kdpi_dat <- dgf_kdpi_direct_data()
    
    validate(
      need(length(unique(logit_dat[[DGF_col]])) == 2, "DGF outcome must contain both 0 and 1 for logistic AUC."),
      need(length(unique(kdpi_dat[[DGF_col]])) == 2, "DGF outcome must contain both 0 and 1 for KDPI AUC."))
    
    pred_logit <- as.numeric(predict(fit_logit, newdata = logit_dat, type = "response"))
    kdpi_score <- kdpi_dat$KDPI_2024
    
    tibble(
      Model = c("Logistic regression", "KDPI_2024 direct score"),
      N = c(nrow(logit_dat), nrow(kdpi_dat)),
      AUC = c(
        calc_auc_safe(logit_dat[[DGF_col]], pred_logit),
        calc_auc_safe(kdpi_dat[[DGF_col]], kdpi_score)))})
  
  dgf_coef_table_reactive <- reactive({
    fit <- dgf_logit_model_all_tx()
    top_coef_table(fit, top_n = input$top_n_dgf)})
  
  output$dgf_perf_table <- renderTable({
    dgf_perf()}, digits = 3)
  
  output$dgf_coef_table <- renderTable({
    dgf_coef_table_reactive()}, digits = 4)
  
  output$dgf_roc_plot <- renderPlot({
    fit_logit <- dgf_logit_model_all_tx()
    logit_dat <- dgf_logit_data_all_tx()
    kdpi_dat <- dgf_kdpi_direct_data()
    
    validate(
      need(length(unique(logit_dat[[DGF_col]])) == 2, "DGF outcome must contain both 0 and 1 for logistic ROC."),
      need(length(unique(kdpi_dat[[DGF_col]])) == 2, "DGF outcome must contain both 0 and 1 for KDPI ROC."))
    
    logit_prob <- as.numeric(predict(fit_logit, newdata = logit_dat, type = "response"))
    kdpi_score <- kdpi_dat$KDPI_2024
    
    roc_logit <- roc_df_safe(logit_dat[[DGF_col]], logit_prob)
    roc_kdpi <- roc_df_safe(kdpi_dat[[DGF_col]], kdpi_score)
    
    auc_logit <- calc_auc_safe(logit_dat[[DGF_col]], logit_prob)
    auc_kdpi <- calc_auc_safe(kdpi_dat[[DGF_col]], kdpi_score)
    
    req(!is.null(roc_logit), !is.null(roc_kdpi))
    
    ggplot() +
      geom_line(data = roc_logit, aes(x = fpr, y = sensitivity), linewidth = 1) +
      geom_line(data = roc_kdpi, aes(x = fpr, y = sensitivity), linewidth = 1, linetype = "dashed") +
      geom_abline(intercept = 0, slope = 1, linetype = "dotted") +
      labs(
        title = "ROC Curve for Delayed Graft Function (DGF)",
        subtitle = paste0(
          "All complete TX subjects. Solid = Logistic regression (AUC = ", round(auc_logit, 3),
          "); Dashed = KDPI_2024 direct score (AUC = ", round(auc_kdpi, 3), ")"),
        x = "False Positive Rate (1 - Specificity)",
        y = "Sensitivity") +
      theme_bw()})
  
  output$recipient_dgf_table <- renderTable({
    tx <- cohorts()$tx
    req("Recipient_code" %in% names(tx), !is.null(input$recipient_code))
    
    fit <- dgf_logit_model_all_tx()
    
    rec <- tx %>%
      mutate(Recipient_code = as.character(Recipient_code)) %>%
      filter(Recipient_code == as.character(input$recipient_code)) %>%
      dplyr::slice(1)
    
    if (nrow(rec) != 1) return(NULL)
    
    rec_base <- rec %>% select(all_of(dgf_predictors)) %>% drop_na()
    if (nrow(rec_base) != 1) {
      return(tibble(Message = "Selected recipient has missing DGF predictors; cannot generate prediction."))}
    
    pred_prob <- as.numeric(predict(fit, newdata = rec_base, type = "response"))
    true_dgf <- if (DGF_col %in% names(rec)) safe_binary_01(rec[[DGF_col]]) else NA_real_
    
    tibble(
      Recipient_code = as.character(rec$Recipient_code),
      Predicted_DGF_Probability = pred_prob,
      Observed_DGF = true_dgf)}, digits = 4)
  
  output$donor_dgf_table <- renderTable({
    non_tx <- cohorts()$non_tx
    req(nrow(non_tx) > 0, !is.null(input$donor_id))
    
    fit <- dgf_logit_model_train_split()
    id_col <- pick_id_col(non_tx)
    
    rec <- non_tx %>% filter(.data[[id_col]] == input$donor_id) %>% dplyr::slice(1)
    if (nrow(rec) != 1) return(NULL)
    
    rec_base <- rec %>% select(all_of(dgf_predictors)) %>% drop_na()
    if (nrow(rec_base) != 1) {
      return(tibble(Message = "Selected Non-TX donor has missing DGF predictors; cannot generate prediction."))}
    
    pred_prob <- as.numeric(predict(fit, newdata = rec_base, type = "response"))
    
    tibble(
      NonTX_ID = as.character(rec[[id_col]]),
      Predicted_DGF_Probability = pred_prob,
      Ground_Truth_Available = "No")}, digits = 4)
  
  # ------------------------------------------------------------
  # Machine learning model comparison
  # ------------------------------------------------------------
  
  ml_results <- eventReactive(input$run_ml, {
    tx <- cohorts()$tx
    validate(need(nrow(tx) >= 20, "Not enough transplanted subjects for ML modeling."))
    
    outcomes <- c(eGFR_outcomes, DGF_col)
    perf_rows <- list()
    fit_rows <- list()
    idx <- 1
    
    withProgress(message = "Fitting machine learning models...", value = 0, {
      total <- length(outcomes) * length(ml_predictor_set_names) * length(ml_model_names)
      step <- 0
      
      for (outcome in outcomes) {
        outcome_type <- ifelse(outcome == DGF_col, "Binary", "Continuous")
        pred_sets <- ml_predictor_sets_for_outcome(
          data_names = names(tx),
          outcome_name = outcome,
          top_n = 100)
        
        for (set_name in names(pred_sets)) {
          predictors <- pred_sets[[set_name]]
          
          for (model_name in ml_model_names) {
            step <- step + 1
            incProgress(1 / total, detail = paste(outcome, set_name, model_name))
            
            fitted <- fit_ml_model(
              dat = tx,
              outcome = outcome,
              predictors = predictors,
              predictor_set = set_name,
              model_name = model_name,
              outcome_type = outcome_type)
            
            if (!is.null(fitted)) {
              perf_rows[[idx]] <- fitted$perf
              fit_rows[[idx]] <- tibble(
                Outcome = outcome,
                Outcome_Type = outcome_type,
                Predictor_Set = set_name,
                Model = model_name,
                Fit = list(fitted$fit))
              idx <- idx + 1}}}}})
    
    list(
      perf = bind_rows(perf_rows),
      fits = bind_rows(fit_rows))}, ignoreInit = TRUE)
  
  ml_perf <- reactive({
    req(input$run_ml > 0)
    res <- ml_results()
    validate(need(nrow(res$perf) > 0, "Something went wrong :("))
    res$perf})
  
  ml_fits <- reactive({
    req(input$run_ml > 0)
    res <- ml_results()
    validate(need(nrow(res$fits) > 0, "Something went wrong :("))
    res$fits})
  
  ml_best <- reactive({
    perf <- ml_perf()
    
    bind_rows(lapply(split(perf, paste(perf$Outcome, perf$Predictor_Set, sep = "||")), function(dd) {
      if (unique(dd$Outcome_Type) == "Continuous") {
        dd %>% filter(!is.na(MSE)) %>% slice_min(order_by = MSE, n = 1, with_ties = FALSE)
      } else {
        dd %>% filter(!is.na(AUC)) %>% slice_max(order_by = AUC, n = 1, with_ties = FALSE)}}))})
  
  output$ml_perf_table <- renderTable({
    ml_perf() %>%
      filter(
        Outcome == input$ml_selected_outcome,
        Predictor_Set == input$ml_selected_predictor_set
      ) %>%
      arrange(Outcome, Predictor_Set, Model) %>%
      select(Outcome, Outcome_Type, Predictor_Set, Model, N, Num_Predictors, MSE, RMSE, AUC, Selection_Metric)}, digits = 3)
  
  output$ml_perf_plot <- renderPlot({
    dd <- ml_perf() %>%
      filter(
        Outcome == input$ml_selected_outcome,
        Predictor_Set == input$ml_selected_predictor_set)
    req(nrow(dd) > 0)
    
    if (unique(dd$Outcome_Type) == "Continuous") {
      ggplot(dd, aes(x = Model, y = MSE)) +
        geom_col() +
        labs(
          title = paste("ML model comparison for", input$ml_selected_outcome),
          subtitle = paste(input$ml_selected_predictor_set, "— lower cross-validated MSE is better"),
          x = "Model",
          y = "MSE") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))} else {
      ggplot(dd, aes(x = Model, y = AUC)) +
        geom_col() +
        labs(
          title = paste("ML model comparison for", input$ml_selected_outcome),
          subtitle = paste(input$ml_selected_predictor_set, "— higher cross-validated AUC is better"),
          x = "Model",
          y = "AUC") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))}})
  
  output$ml_best_table <- renderTable({
    ml_best() %>%
      arrange(match(Outcome, c(eGFR_outcomes, DGF_col)), Predictor_Set) %>%
      select(Outcome, Outcome_Type, Predictor_Set, Model, N, Num_Predictors, MSE, RMSE, AUC, Selection_Metric)}, digits = 3)
  output$ml_best_msg <- renderUI({
    best <- ml_best() %>%
      filter(
        Outcome == input$ml_selected_outcome,
        Predictor_Set == input$ml_selected_predictor_set)
    req(nrow(best) == 1)
    
    metric_value <- if (best$Outcome_Type[1] == "Continuous") best$MSE[1] else best$AUC[1]
    metric_name <- if (best$Outcome_Type[1] == "Continuous") "MSE" else "AUC"
    direction <- if (best$Outcome_Type[1] == "Continuous") "lowest" else "highest"
    
    tags$p(
      style = "font-weight:600; color:#1b5e20;",
      paste0(
        "Best model for ", best$Outcome[1], " using ", best$Predictor_Set[1], ": ",
        best$Model[1], " selected by ", direction, " ", metric_name,
        " = ", round(metric_value, 3), ". Final prediction model is refit on all complete TX rows."))})
  
  # ------------------------------------------------------------
  # Point prediction kidney recommendation
  # ------------------------------------------------------------
  
  output$kidney_recommendation_table <- renderTable({
    validate(
      need(input$run_ml > 0,
           "Click 'Run / Refresh ML models' before generating a recommendation."))
    
    non_tx <- cohorts()$non_tx
    req(nrow(non_tx) > 0, !is.null(input$donor_id))
    
    id_col <- pick_id_col(non_tx)
    
    rec <- non_tx %>%
      filter(.data[[id_col]] == input$donor_id) %>%
      dplyr::slice(1)
    
    if (nrow(rec) != 1) return(NULL)
    
    best <- ml_best()
    fits <- ml_fits()
    
    preferred_set <- "Joint MRMR-selected donor clinical + image features"
    
    best_egfr <- best %>%
      filter(
        Outcome == input$rec_egfr_outcome,
        Predictor_Set == preferred_set) %>%dplyr::slice(1)
    
    best_dgf <- best %>%
      filter(
        Outcome == DGF_col,
        Predictor_Set == preferred_set) %>%
      dplyr::slice(1)
    
    validate(
      need(nrow(best_egfr) == 1,
           "No best eGFR ML model available for the recommendation."),
      need(nrow(best_dgf) == 1,
           "No best DGF ML model available for the recommendation."))
    
    egfr_fit <- fits %>%
      filter(
        Outcome == best_egfr$Outcome[1],
        Predictor_Set == best_egfr$Predictor_Set[1],
        Model == best_egfr$Model[1]) %>%
      dplyr::slice(1) %>%
      pull(Fit) %>%
      .[[1]]
    
    dgf_fit <- fits %>%
      filter(
        Outcome == best_dgf$Outcome[1],
        Predictor_Set == best_dgf$Predictor_Set[1],
        Model == best_dgf$Model[1]) %>%
      dplyr::slice(1) %>%
      pull(Fit) %>%
      .[[1]]
    
    pred_egfr <- predict_ml_fit(egfr_fit, rec)
    pred_dgf  <- predict_ml_fit(dgf_fit, rec)
    
    recommendation <- ifelse(
      !is.na(pred_egfr) &&
        !is.na(pred_dgf) &&
        pred_egfr >= input$rec_egfr_cutoff &&
        pred_dgf <= input$rec_dgf_cutoff,
      "Yes - meets selected thresholds",
      "No / Caution - does not meet selected thresholds")
    
    tibble(
      NonTX_ID = as.character(rec[[id_col]]),
      Best_eGFR_Model = best_egfr$Model[1],
      eGFR_Outcome = input$rec_egfr_outcome,
      Predicted_eGFR = pred_egfr,
      Predicted_CKD_Stage = ckd_stage(pred_egfr),
      eGFR_Cutoff = input$rec_egfr_cutoff,
      Best_DGF_Model = best_dgf$Model[1],
      Predicted_DGF_Probability = pred_dgf,
      DGF_Cutoff = input$rec_dgf_cutoff,
      Kidney_Recommendation = recommendation,
      Ground_Truth_Available = "No")}, digits = 3)}


shinyApp(ui = ui, server = server)
