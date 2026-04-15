# app.R
# Kidney Transplant Outcomes Dashboard
# Abby Sun

library(shiny)
library(tidyverse)
library(pROC)

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

stage_levels <- c(
  "G1 (≥90)", "G2 (60–89)", "G3a (45–59)",
  "G3b (30–44)", "G4 (15–29)", "G5 (<15)")

DATA_PATH <- "Reduced_Features_NO_COLOR 2.csv"


ckd_stage <- function(e) {
  dplyr::case_when(
    is.na(e) ~ NA_character_,
    e >= 90  ~ "G1 (≥90)",
    e >= 60  ~ "G2 (60–89)",
    e >= 45  ~ "G3a (45–59)",
    e >= 30  ~ "G3b (30–44)",
    e >= 15  ~ "G4 (15–29)",
    TRUE     ~ "G5 (<15)")}

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
# UI
# ============================================================
ui <- fluidPage(
  titlePanel("Kidney Transplant Outcomes Dashboard"),
  
  sidebarLayout(
    sidebarPanel(
      helpText("eGFR models use an 80/20 hold-out split. DGF uses logistic regression on all available transplant recipients."),
      
      radioButtons(
        "model_mode", "eGFR model type",
        choices = c(
          "KDPI-only (80/20 split)" = "kdpi_split",
          "Multivariable predictors (80/20 split)" = "multi_cv"),
        selected = "multi_cv"),
      
      numericInput("seed", "Random seed", value = 123, min = 1, step = 1),
      
      conditionalPanel(
        condition = "input.model_mode == 'kdpi_split'",
        sliderInput("train_frac", "Training fraction", min = 0.5, max = 0.9, value = 0.8, step = 0.05)),
      
      conditionalPanel(
        condition = "input.model_mode == 'multi_cv'",
        sliderInput("train_frac_multi", "Training fraction", min = 0.5, max = 0.9, value = 0.8, step = 0.05)),
      
      hr(),
      
      selectInput(
        "cohort_view",
        "CKD Stage Plot Cohort",
        choices = c("Training", "Validation", "Non-Transplant (Prediction only)"),
        selected = "Validation"),
      
      selectInput(
        "outcome_view",
        "Outcome for CKD Stage Plot / Pred vs True",
        choices = c("All outcomes (faceted)", eGFR_outcomes),
        selected = "All outcomes (faceted)"),
      
      hr(),
      h4("DGF logistic regression"),
      sliderInput("top_n_dgf", "Top DGF predictors to display", min = 3, max = 10, value = 5, step = 1),
      
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
          h4("Predictor coefficients (lm on TRAINING set) for selected outcome"),
          tags$p(style = "color:#555;", "Select a single outcome in the sidebar to populate this."),
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
            tags$li("Training/Validation cohorts = TX rows (observed outcomes available)."),
            tags$li("Non-Transplant cohort = prediction-only rows (no observed recipient eGFR outcomes)."),
            tags$li("Do not interpret Non-Transplant predictions as truth — they are model projections only."))),
        
        tabPanel(
          "Individual Recipient",
          br(),
          h4("Predicted vs true eGFR + CKD stage (intervals shown for multivariable mode)"),
          tableOutput("recipient_compare_table")),
        
        tabPanel(
          "Individual Non-TX Donor Predictions",
          br(),
          tags$p(
            style = "color:#b00;",
            "Prediction-only: Non-TX rows do not have observed recipient eGFR outcomes."),
          tableOutput("donor_pred_table")),
        
        tabPanel(
          "DGF Performance",
          br(),
          h4("Logistic regression for delayed graft function (DGF)"),
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
          h4("Selected recipient predicted probability of DGF"),
          tableOutput("recipient_dgf_table"),
          br(),
          h4("Selected Non-TX donor predicted probability of DGF"),
          tags$p(
            style = "color:#b00;",
            "Prediction-only: Non-TX rows do not have observed recipient DGF outcomes."),
          tableOutput("donor_dgf_table"))))))

# ============================================================
# SERVER
# ============================================================
server <- function(input, output, session) {
  
  observe({
    if (!file.exists(DATA_PATH)) {
      showNotification(
        paste0(
          "ERROR: Data file not found: '", DATA_PATH, "'. ",
          "Make sure the CSV is in the same folder as app.R and deployed with the app."),
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
    
    for (x in predictors_full) {
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
    if (input$model_mode == "kdpi_split") "KDPI_2024" else predictors_full})
  
  # ============================================================
  # MODE A: KDPI-only 80/20 split
  # ============================================================
  split_tx_kdpi <- reactive({
    req(input$model_mode == "kdpi_split")
    set.seed(input$seed)
    
    tx <- cohorts()$tx
    n <- nrow(tx)
    if (n < 10) return(NULL)
    
    train_n <- floor(input$train_frac * n)
    idx <- sample(seq_len(n), size = train_n)
    
    list(
      train = tx[idx, , drop = FALSE],
      test  = tx[-idx, , drop = FALSE])})
  
  kdpi_models <- reactive({
    req(input$model_mode == "kdpi_split")
    sp <- split_tx_kdpi()
    if (is.null(sp)) return(list())
    
    train_df <- sp$train
    mods <- list()
    
    for (outcome in eGFR_outcomes) {
      sub_train <- train_df %>%
        select(KDPI_2024, all_of(outcome)) %>%
        drop_na()
      
      if (nrow(sub_train) < 8) next
      mods[[outcome]] <- lm(as.formula(paste(outcome, "~ KDPI_2024")), data = sub_train)} 
    mods})
  
  perf_kdpi <- reactive({
    req(input$model_mode == "kdpi_split")
    sp <- split_tx_kdpi()
    mods <- kdpi_models()
    if (is.null(sp) || length(mods) == 0) return(NULL)
    
    rows <- lapply(names(mods), function(outcome) {
      model <- mods[[outcome]]
      
      train_sub <- sp$train %>% select(KDPI_2024, all_of(outcome)) %>% drop_na()
      test_sub  <- sp$test  %>% select(KDPI_2024, all_of(outcome)) %>% drop_na()
      
      if (nrow(train_sub) < 5 || nrow(test_sub) < 5) return(NULL)
      
      train_pred <- as.numeric(predict(model, newdata = train_sub))
      test_pred  <- as.numeric(predict(model, newdata = test_sub))
      
      train_mse <- mean((train_sub[[outcome]] - train_pred)^2)
      val_mse   <- mean((test_sub[[outcome]]  - test_pred)^2)
      
      tibble(
        Outcome = outcome,
        Train_N = nrow(train_sub),
        Validation_N = nrow(test_sub),
        Train_MSE = train_mse,
        Train_RMSE = rmse_from_mse(train_mse),
        Validation_MSE = val_mse,
        Validation_RMSE = rmse_from_mse(val_mse))})
    
    bind_rows(rows) %>% arrange(match(Outcome, eGFR_outcomes))})
  
  # ============================================================
  # MODE B: Multivariable OLS with 80/20 split
  # ============================================================
  split_tx_multi <- reactive({
    req(input$model_mode == "multi_cv")
    set.seed(input$seed)
    
    tx <- cohorts()$tx
    n <- nrow(tx)
    if (n < 10) return(NULL)
    
    train_n <- floor(input$train_frac_multi * n)
    idx <- sample(seq_len(n), size = train_n)
    
    list(
      train = tx[idx, , drop = FALSE],
      test  = tx[-idx, , drop = FALSE])})
  
  multi_models <- reactive({
    req(input$model_mode == "multi_cv")
    sp <- split_tx_multi()
    if (is.null(sp)) return(list())
    
    preds <- predictors_active()
    train_df <- sp$train
    mods <- list()
    
    for (outcome in eGFR_outcomes) {
      needed <- c(outcome, preds)
      sub_train <- train_df %>% select(all_of(needed)) %>% drop_na()
      
      if (nrow(sub_train) < (length(preds) + 8)) next
      
      f <- as.formula(paste(outcome, "~", paste(preds, collapse = " + ")))
      mods[[outcome]] <- list(
        formula = f,
        fit_train = lm(f, data = sub_train),
        train_data = sub_train)}
    mods})
  
  perf_multi <- reactive({
    req(input$model_mode == "multi_cv")
    sp <- split_tx_multi()
    mods <- multi_models()
    if (is.null(sp) || length(mods) == 0) return(NULL)
    
    preds <- predictors_active()
    
    rows <- lapply(names(mods), function(outcome) {
      obj <- mods[[outcome]]
      fit <- obj$fit_train
      
      train_sub <- sp$train %>% select(all_of(c(outcome, preds))) %>% drop_na()
      test_sub  <- sp$test  %>% select(all_of(c(outcome, preds))) %>% drop_na()
      
      if (nrow(train_sub) < 5 || nrow(test_sub) < 5) return(NULL)
      
      train_pred <- as.numeric(predict(fit, newdata = train_sub))
      test_pred  <- as.numeric(predict(fit, newdata = test_sub))
      
      train_mse <- mean((train_sub[[outcome]] - train_pred)^2)
      val_mse   <- mean((test_sub[[outcome]]  - test_pred)^2)
      
      tibble(
        Outcome = outcome,
        Train_N = nrow(train_sub),
        Validation_N = nrow(test_sub),
        Train_MSE = train_mse,
        Train_RMSE = rmse_from_mse(train_mse),
        Validation_MSE = val_mse,
        Validation_RMSE = rmse_from_mse(val_mse))})
    
    bind_rows(rows) %>% arrange(match(Outcome, eGFR_outcomes))})
  
  # ============================================================
  # Availability table
  # ============================================================
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
    
    if (input$model_mode == "kdpi_split") {
      mods <- kdpi_models()
      if (!(sel %in% names(mods))) return(NULL)
      fit <- mods[[sel]]} 
    else {
      mods <- multi_models()
      if (!(sel %in% names(mods))) return(NULL)
      fit <- mods[[sel]]$fit_train}
    
    co <- coef(fit)
    co <- co[names(co) != "(Intercept)"]
    
    tibble(
      Predictor = names(co),
      Coefficient = as.numeric(co),
      Abs = abs(as.numeric(co))) %>%
      arrange(desc(Abs)) %>%
      select(Predictor, Coefficient)})
  
  # ============================================================
  # eGFR outputs
  # ============================================================
  output$perf_table <- renderTable({
    if (input$model_mode == "kdpi_split") perf_kdpi() else perf_multi()}, digits = 3)
  
  output$availability_table <- renderTable({
    availability_table()}, digits = 0)
  
  output$coef_table <- renderTable({
    coef_table()}, digits = 4)
  
  output$mse_plot <- renderPlot({
    perf <- if (input$model_mode == "kdpi_split") perf_kdpi() else perf_multi()
    req(perf)
    
    perf_long <- perf %>%
      select(Outcome, Train_MSE, Validation_MSE) %>%
      pivot_longer(
        cols = c(Train_MSE, Validation_MSE),
        names_to = "Dataset",
        values_to = "MSE") %>%
      mutate(Dataset = recode(Dataset, Train_MSE = "Training", Validation_MSE = "Validation"))
    
    ggplot(perf_long, aes(x = Outcome, y = MSE, fill = Dataset)) +
      geom_col(position = "dodge") +
      labs(
        title = if (input$model_mode == "kdpi_split") {
          "Train vs Validation MSE (KDPI-only OLS)"} 
        else {"Train vs Validation MSE (Multivariable OLS)"},
        x = "Outcome",
        y = "MSE (eGFR^2)") +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))})
  
  # ============================================================
  # CKD stage distributions
  # ============================================================
  pred_long <- reactive({
    tx <- cohorts()$tx
    non_tx <- cohorts()$non_tx
    
    if (input$model_mode == "kdpi_split") {
      sp <- split_tx_kdpi()
      mods <- kdpi_models()
      req(sp)
      if (length(mods) == 0) return(NULL)
      
      make_pred_df <- function(dat, label) {
        out_list <- list()
        for (outcome in names(mods)) {
          model <- mods[[outcome]]
          sub <- dat %>% select(KDPI_2024) %>% drop_na()
          if (nrow(sub) == 0) next
          
          pred <- suppressWarnings(as.numeric(predict(model, newdata = sub)))
          
          out_list[[outcome]] <- tibble(
            Cohort = label,
            Outcome = outcome,
            Pred_eGFR = pred,
            Pred_CKD_Stage = factor(ckd_stage(pred), levels = stage_levels))}
        bind_rows(out_list)}
      
      bind_rows(
        make_pred_df(sp$train, "Training"),
        make_pred_df(sp$test, "Validation"),
        make_pred_df(non_tx, "Non-Transplant (Prediction only)")) %>% filter(!is.na(Pred_CKD_Stage))} 
    else {sp <- split_tx_multi()
      mods <- multi_models()
      preds <- predictors_active()
      req(sp)
      if (length(mods) == 0) return(NULL)
      
      make_pred_df <- function(dat, label) {
        out_list <- list()
        base <- dat %>% select(all_of(preds)) %>% drop_na()
        if (nrow(base) == 0) return(tibble())
        
        for (outcome in names(mods)) {
          fit <- mods[[outcome]]$fit_train
          pred <- suppressWarnings(as.numeric(predict(fit, newdata = base)))
          
          out_list[[outcome]] <- tibble(
            Cohort = label,
            Outcome = outcome,
            Pred_eGFR = pred,
            Pred_CKD_Stage = factor(ckd_stage(pred), levels = stage_levels))}
        bind_rows(out_list)}
      
      bind_rows(
        make_pred_df(sp$train, "Training"),
        make_pred_df(sp$test, "Validation"),
        make_pred_df(non_tx, "Non-Transplant (Prediction only)")) %>% filter(!is.na(Pred_CKD_Stage))}})
  
  output$ckd_plot <- renderPlot({
    pdat <- pred_long()
    req(pdat)
    
    cohort_sel <- input$cohort_view
    pdat <- pdat %>% filter(Cohort == cohort_sel)
    
    if (input$outcome_view == "All outcomes (faceted)") {
      ggplot(pdat, aes(x = Pred_CKD_Stage)) +
        geom_bar() +
        facet_wrap(~ Outcome, scales = "free_y") +
        labs(
          title = paste("Predicted CKD Stage Distribution –", cohort_sel),
          x = "Predicted CKD Stage",
          y = "Count") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))} 
    else {pdat2 <- pdat %>% filter(Outcome == input$outcome_view)
      
      ggplot(pdat2, aes(x = Pred_CKD_Stage)) +
        geom_bar() +
        labs(
          title = paste("Predicted CKD Stage –", cohort_sel, "-", input$outcome_view),
          x = "Predicted CKD Stage",
          y = "Count") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))}})
  
  # ============================================================
  # Recipient selector + comparison table
  # ============================================================
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
      slice(1)
    
    if (nrow(rec) != 1) return(NULL)
    
    if (input$model_mode == "kdpi_split") {
      mods <- kdpi_models()
      if (length(mods) == 0) {
        return(tibble(Message = "No KDPI-only models available (insufficient data)."))}
      
      kdpi_val <- safe_as_numeric(rec$KDPI_2024[1])
      if (is.na(kdpi_val)) {
        return(tibble(Message = "Selected recipient is missing KDPI_2024, cannot predict in KDPI-only mode."))}
      
      rows <- lapply(eGFR_outcomes, function(outcome) {
        model <- mods[[outcome]]
        if (is.null(model)) return(NULL)
        
        pred <- as.numeric(predict(model, newdata = tibble(KDPI_2024 = kdpi_val)))
        true <- safe_as_numeric(rec[[outcome]])
        
        tibble(
          Recipient_code = as.character(rec$Recipient_code),
          Outcome = outcome,
          Predicted_eGFR = pred,
          Predicted_CKD_Stage = ckd_stage(pred),
          True_eGFR = true,
          True_CKD_Stage = ckd_stage(true))})
      bind_rows(rows)} 
    else {
      mods <- multi_models()
      preds <- predictors_active()
      
      if (length(mods) == 0) {
        return(tibble(Message = "No multivariable models available (insufficient complete data)."))}
      
      rec_pred_base <- rec %>% select(all_of(preds)) %>% drop_na()
      if (nrow(rec_pred_base) != 1) {
        return(tibble(Message = "Selected recipient has missing predictors; cannot use multivariable model for individualized prediction."))}
      
      rows <- lapply(eGFR_outcomes, function(outcome) {
        obj <- mods[[outcome]]
        if (is.null(obj)) return(NULL)
        
        ints <- predict_with_intervals(obj$fit_train, rec_pred_base)
        true <- safe_as_numeric(rec[[outcome]])
        
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
          True_CKD_Stage = ckd_stage(true))})
      
      bind_rows(rows)}}, digits = 2)
  
  # ============================================================
  # Non-TX donor selector + prediction table
  # ============================================================
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
    rec <- non_tx %>% filter(.data[[id_col]] == input$donor_id) %>% slice(1)
    if (nrow(rec) != 1) return(NULL)
    
    if (input$model_mode == "kdpi_split") {
      mods <- kdpi_models()
      if (length(mods) == 0) return(tibble(Message = "No KDPI-only models available."))
      
      kdpi_val <- safe_as_numeric(rec$KDPI_2024[1])
      if (is.na(kdpi_val)) {
        return(tibble(Message = "Selected Non-TX row is missing KDPI_2024; cannot predict in KDPI-only mode."))}
      
      rows <- lapply(eGFR_outcomes, function(outcome) {
        model <- mods[[outcome]]
        if (is.null(model)) return(NULL)
        
        pred <- as.numeric(predict(model, newdata = tibble(KDPI_2024 = kdpi_val)))
        
        tibble(
          NonTX_ID = as.character(rec[[id_col]]),
          KDPI_2024 = kdpi_val,
          Outcome = outcome,
          Predicted_eGFR = pred,
          Predicted_CKD_Stage = ckd_stage(pred))})
      
      bind_rows(rows)} 
    else {
      mods <- multi_models()
      preds <- predictors_active()
      if (length(mods) == 0) return(tibble(Message = "No multivariable models available."))
      
      rec_pred_base <- rec %>% select(all_of(preds)) %>% drop_na()
      if (nrow(rec_pred_base) != 1) {
        return(tibble(Message = "Selected Non-TX row has missing predictors; cannot use multivariable model for individualized prediction."))}
      
      rows <- lapply(eGFR_outcomes, function(outcome) {
        obj <- mods[[outcome]]
        if (is.null(obj)) return(NULL)
        
        ints <- predict_with_intervals(obj$fit_train, rec_pred_base)
        
        tibble(
          NonTX_ID = as.character(rec[[id_col]]),
          Outcome = outcome,
          Predicted_eGFR = ints$Predicted_eGFR,
          Pred_CI_L = ints$Pred_CI_L,
          Pred_CI_U = ints$Pred_CI_U,
          Pred_PI_L = ints$Pred_PI_L,
          Pred_PI_U = ints$Pred_PI_U,
          Predicted_CKD_Stage = ckd_stage(ints$Predicted_eGFR))})
      
      bind_rows(rows)}}, digits = 2)
  
  # ============================================================
  # Pred vs True plot
  # ============================================================
  output$pred_true_msg <- renderUI({
    if (input$outcome_view == "All outcomes (faceted)") {
      tags$p(
        style = "color:#b00;",
        "Select a single outcome (e.g. eGFR_CKD_EPI_6M) to view Predicted vs True.")} else {
      NULL}})
  
  output$pred_true_plot <- renderPlot({
    sel <- input$outcome_view
    req(sel != "All outcomes (faceted)")
    
    if (input$model_mode == "kdpi_split") {
      sp <- split_tx_kdpi()
      mods <- kdpi_models()
      req(sp)
      
      model <- mods[[sel]]
      req(!is.null(model))
      
      dat <- sp$test %>% select(KDPI_2024, all_of(sel)) %>% drop_na()
      req(nrow(dat) >= 10)
      
      pred <- as.numeric(predict(model, newdata = dat))
      dplt <- tibble(Pred = pred, True = dat[[sel]])
      
      cal_fit <- lm(True ~ Pred, data = dplt)
      r2 <- summary(cal_fit)$r.squared
      
      ggplot(dplt, aes(x = Pred, y = True)) +
        geom_point(alpha = 0.7) +
        geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
        geom_smooth(method = "lm", se = FALSE) +
        labs(
          title = paste("Predicted vs True eGFR —", sel, "(KDPI-only)"),
          subtitle = paste0(
            "Validation (hold-out). Dashed = y = x line. Solid = calibration. R² = ",
            round(r2, 3)),
          x = "Predicted eGFR",
          y = "True eGFR") +
        theme_bw()} 
    else {
      sp <- split_tx_multi()
      mods <- multi_models()
      preds <- predictors_active()
      req(sp)
      
      obj <- mods[[sel]]
      req(!is.null(obj))
      fit <- obj$fit_train
      
      dat <- sp$test %>% select(all_of(c(sel, preds))) %>% drop_na()
      req(nrow(dat) >= 10)
      
      pred <- as.numeric(predict(fit, newdata = dat))
      dplt <- tibble(Pred = pred, True = dat[[sel]])
      
      cal_fit <- lm(True ~ Pred, data = dplt)
      r2 <- summary(cal_fit)$r.squared
      
      ggplot(dplt, aes(x = Pred, y = True)) +
        geom_point(alpha = 0.7) +
        geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
        geom_smooth(method = "lm", se = FALSE) +
        labs(
          title = paste("Predicted vs True eGFR —", sel, "(Multivariable OLS)"),
          subtitle = paste0(
            "Validation (hold-out). Dashed = y = x line. Solid = calibration. R² = ",
            round(r2, 3)),
          x = "Predicted eGFR",
          y = "True eGFR") +
        theme_bw()}})
  
  # ============================================================
  # DGF logistic regression
  # ============================================================
  dgf_tx_data <- reactive({
    d <- cohorts()$tx
    req(DGF_col %in% names(d))
    
    d %>% filter(!is.na(.data[[DGF_col]]))})
  
  dgf_logit_data <- reactive({
    dgf_tx_data() %>%
      select(all_of(c(DGF_col, dgf_predictors))) %>%
      drop_na()})
  
  dgf_logit_model <- reactive({
    dat <- dgf_logit_data()
    
    validate(
      need(nrow(dat) >= 15, "Not enough complete DGF cases to fit logistic regression."),
      need(length(unique(dat[[DGF_col]])) == 2, "DGF outcome must contain both 0 and 1."))
    
    f <- as.formula(paste(DGF_col, "~", paste(dgf_predictors, collapse = " + ")))
    glm(f, data = dat, family = binomial())})
  
  dgf_kdpi_data <- reactive({
    dgf_tx_data() %>%
      select(all_of(c(DGF_col, "KDPI_2024"))) %>%
      drop_na()})
  
  dgf_perf <- reactive({
    fit <- dgf_logit_model()
    logit_dat <- dgf_logit_data()
    kdpi_dat <- dgf_kdpi_data()
    
    pred_prob <- as.numeric(predict(fit, newdata = logit_dat, type = "response"))
    
    auc_logit <- calc_auc_safe(logit_dat[[DGF_col]], pred_prob)
    auc_kdpi_direct <- calc_auc_safe(kdpi_dat[[DGF_col]], kdpi_dat$KDPI_2024)
    
    tibble(
      Model = c("Logistic regression", "KDPI_2024 direct score"),
      N = c(nrow(logit_dat), nrow(kdpi_dat)),
      AUC = c(auc_logit, auc_kdpi_direct))})
  
  dgf_coef_table_reactive <- reactive({
    fit <- dgf_logit_model()
    top_coef_table(fit, top_n = input$top_n_dgf)})
  
  output$dgf_perf_table <- renderTable({
    dgf_perf()}, digits = 3)
  
  output$dgf_coef_table <- renderTable({
    dgf_coef_table_reactive()}, digits = 4)
  
  output$dgf_roc_plot <- renderPlot({
    fit <- dgf_logit_model()
    
    logit_dat <- dgf_logit_data()
    kdpi_dat <- dgf_kdpi_data()
    
    logit_prob <- as.numeric(predict(fit, newdata = logit_dat, type = "response"))
    
    roc_logit <- roc_df_safe(logit_dat[[DGF_col]], logit_prob)
    roc_kdpi <- roc_df_safe(kdpi_dat[[DGF_col]], kdpi_dat$KDPI_2024)
    
    auc_logit <- calc_auc_safe(logit_dat[[DGF_col]], logit_prob)
    auc_kdpi <- calc_auc_safe(kdpi_dat[[DGF_col]], kdpi_dat$KDPI_2024)
    
    req(!is.null(roc_logit), !is.null(roc_kdpi))
    
    ggplot() +
      geom_line(data = roc_logit, aes(x = fpr, y = sensitivity), linewidth = 1) +
      geom_line(data = roc_kdpi, aes(x = fpr, y = sensitivity), linewidth = 1, linetype = "dashed") +
      geom_abline(intercept = 0, slope = 1, linetype = "dotted") +
      labs(
        title = "ROC Curve for Delayed Graft Function (DGF)",
        subtitle = paste0(
          "Solid = Logistic regression (AUC = ", round(auc_logit, 3),
          "); Dashed = KDPI_2024 direct score (AUC = ", round(auc_kdpi, 3), ")"),
        x = "False Positive Rate (1 - Specificity)",
        y = "Sensitivity") +
      theme_bw()})
  
  output$recipient_dgf_table <- renderTable({
    tx <- cohorts()$tx
    req("Recipient_code" %in% names(tx), !is.null(input$recipient_code))
    
    fit <- dgf_logit_model()
    
    rec <- tx %>%
      mutate(Recipient_code = as.character(Recipient_code)) %>%
      filter(Recipient_code == as.character(input$recipient_code)) %>%
      slice(1)
    
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
    
    fit <- dgf_logit_model()
    id_col <- pick_id_col(non_tx)
    
    rec <- non_tx %>% filter(.data[[id_col]] == input$donor_id) %>% slice(1)
    if (nrow(rec) != 1) return(NULL)
    
    rec_base <- rec %>% select(all_of(dgf_predictors)) %>% drop_na()
    if (nrow(rec_base) != 1) {
      return(tibble(Message = "Selected Non-TX donor has missing DGF predictors; cannot generate prediction."))}
    
    pred_prob <- as.numeric(predict(fit, newdata = rec_base, type = "response"))
    
    tibble(
      NonTX_ID = as.character(rec[[id_col]]),
      Predicted_DGF_Probability = pred_prob)}, digits = 4)}

shinyApp(ui = ui, server = server)