library(readxl)
library(ggplot2)


plot_outcome_eGFR <- function(sheet_name) {
  df <- read_excel("ml_sensitivity.xlsx", sheet = sheet_name)
  month_order <- c(
    "eGFR_CKD_EPI_1M",
    "eGFR_CKD_EPI_3M",
    "eGFR_CKD_EPI_6M",
    "eGFR_CKD_EPI_12M",
    "eGFR_CKD_EPI_3years",
    "eGFR_CKD_EPI_5years"
  )
  df$Month <- factor(df$Month, levels = month_order)
  
ggplot(df, aes(x = Month, y = MSE)) +
  geom_col(fill = "steelblue") +
  geom_text(aes(label = MSE), vjust = -0.5) +
  labs(title = paste("eGFR MSE by", sheet_name), x = "Outcome", y = "MSE") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5)  # 0 = left, 0.5 = center, 1 = right
  )
}

plot_outcome_eGFR("Elastic Net")
plot_outcome_eGFR("Lasso")
plot_outcome_eGFR("RF")
plot_outcome_eGFR("Ridge")
plot_outcome_eGFR("XGBoost")


#I think the dashboard already plots this, but I just did it to learn LOL
plot_outcome_DGF <- function(sheet_name) {
  df <- read_excel("ml_sensitivity.xlsx", sheet = sheet_name)
  model_order <- c(
    "Elastic Net",
    "Lasso",
    "RF",
    "Ridge"
  )
  df$Model <- factor(df$Model, levels = model_order)
  
  ggplot(df, aes(x = Model, y = AUC)) +
    geom_col(fill = "steelblue") +
    geom_text(aes(label = AUC), vjust = -0.5) +
    labs(title = paste("ML Model Comparison for", sheet_name), x = "Model", y = "AUC") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(hjust = 0.5)  # 0 = left, 0.5 = center, 1 = right
    )
}


plot_outcome_DGF("DGF")
