# 04_ml_models.R
# Train and evaluate Random Forest, Logistic Regression, SVM, and Elastic Net models.

suppressPackageStartupMessages({
  library(caret)
  library(readr)
  library(dplyr)
  library(pROC)
})

root_dir <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
prep <- readRDS(file.path(root_dir, "data", "raw", "FRANZOSA_IBD_2019", "franzosa2019_preprocessed.rds"))

train_x <- prep$train_data
test_x  <- prep$test_data
train_y <- prep$train_y
test_y  <- prep$test_y

train_df <- cbind(train_x, Group = train_y)
test_df  <- cbind(test_x, Group = test_y)

message(sprintf("Training: %d Samples, %d Features", nrow(train_df), ncol(train_df) - 1))


#Cross Validation Setup

ctrl <- trainControl(
  method = "cv",
  number = 10,
  classProbs = TRUE,
  savePredictions = "final"
)

#Die 4 Modelle trainieren

#Modell 1: Random Forest

set.seed(42)
message("Training Random Forest...")
model_rf <- train(
  Group ~ .,
  data = train_df,
  method = "rf",
  trControl = ctrl,
  metric = "Accuracy",
  tuneLength = 5
)
message("Random Forest fertig!")


#Modell 2: Logistische Regression

set.seed(42)
message("Training Logistic Regression (glmnet multinomial)...")
model_glm <- train(
  Group ~ .,
  data = train_df,
  method = "glmnet",
  trControl = ctrl,
  metric = "Accuracy",
  tuneGrid = expand.grid(alpha = 0, lambda = 10^seq(-3, 1, length = 10))
)
message("Logistic Regression fertig!")


#Modell 3: Support Vector Machine (SVM) 

set.seed(42)
message("Training SVM...")
model_svm <- train(
  Group ~ .,
  data = train_df,
  method = "svmRadial",
  trControl = ctrl,
  metric = "Accuracy",
  tuneLength = 5
)
message("SVM fertig!")

#Modell 4: Elastic Net

set.seed(42)
message("Training Elastic Net...")
model_enet <- train(
  Group ~ .,
  data = train_df,
  method = "glmnet",
  trControl = ctrl,
  metric = "Accuracy",
  tuneLength = 10
)
message("Elastic Net fertig!")


#Evaluierung auf dem Testset

models <- list(
  RandomForest = model_rf,
  LogisticRegression = model_glm,
  SVM = model_svm,
  ElasticNet = model_enet
)

evaluate_model <- function(model, test_df) {
  pred <- predict(model, newdata = test_df)
  cm <- confusionMatrix(pred, test_df$Group)
  out <- tibble(
    model = model$method,
    accuracy = unname(cm$overall["Accuracy"]),
    kappa = unname(cm$overall["Kappa"])
  )
  list(summary = out, confusion = cm$table)
}

evals <- lapply(models, evaluate_model, test_df = test_df)
perf_table <- bind_rows(lapply(evals, `[[`, "summary"))
perf_table <- perf_table %>% mutate(model_name = names(models)) %>% select(model_name, everything())

print(perf_table)


#Speichern und ROC-Kurven


fig_dir <- file.path(root_dir, "output", "figures")
tab_dir <- file.path(root_dir, "output", "tables")

write_csv(perf_table, file.path(tab_dir, "ml_model_performance.csv"))

for (nm in names(evals)) {
  conf_df <- as.data.frame(evals[[nm]]$confusion)
  write_csv(conf_df, file.path(tab_dir, paste0("ml_confusion_", nm, ".csv")))
}

saveRDS(models, file.path(root_dir, "data", "raw", "FRANZOSA_IBD_2019", "franzosa2019_models.rds"))

for (nm in names(models)) {
  mod <- models[[nm]]
  probs <- predict(mod, newdata = test_df, type = "prob")
  classes <- colnames(probs)
  
  png(file.path(fig_dir, paste0("ml_roc_", nm, ".png")), width = 1400, height = 1000, res = 180)
  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(0, 1))
  axis(1); axis(2); box()
  title(main = paste("One-vs-Rest ROC -", nm), xlab = "False Positive Rate", ylab = "True Positive Rate")
  abline(0, 1, lty = 2, col = "gray50")
  
  cols <- rainbow(length(classes))
  leg_txt <- character(0)
  for (i in seq_along(classes)) {
    cls <- classes[i]
    truth <- as.numeric(test_df$Group == cls)
    roc_obj <- roc(truth, probs[[cls]], quiet = TRUE)
    lines(1 - roc_obj$specificities, roc_obj$sensitivities, col = cols[i], lwd = 2)
    leg_txt <- c(leg_txt, sprintf("%s (AUC=%.3f)", cls, as.numeric(auc(roc_obj))))
  }
  legend("bottomright", legend = leg_txt, col = cols, lwd = 2, cex = 0.8)
  dev.off()
}

message("ML training and evaluation completed!")
