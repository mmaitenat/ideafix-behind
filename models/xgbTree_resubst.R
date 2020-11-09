trainxgbTree <- function(X, Y, id) {
  set.seed(id)
  
  # these are xgbtree's default parameters
  tunegrid <- expand.grid(nrounds = 100,
                          max_depth = 6,
                          eta = 0.3,
                          gamma = 0,
                          colsample_bytree = 1,
                          min_child_weight = 1,
                          subsample = 1)
  
  fitControl.noCV <- trainControl(
    method = "none",
    summaryFunction = twoClassSummary,
    classProbs = TRUE,
    verboseIter = TRUE)
  
  
  resubst.model <- caret::train(x = as.matrix(X), 
                        y = Y, 
                        method = "xgbTree", 
                        tuneGrid = tunegrid,
                        trControl = fitControl.noCV, 
                        preProcess = c("center", "scale"),
                        metric = "ROC")
  
  resubst.pred <- predict(resubst.model, as.matrix(X), type = "prob")
  resubst.out <- pROC::roc(Y, resubst.pred[, 1], plot = FALSE)
  
  # save
  saveRDS(resubst.out, file = sprintf("xgbTree_resubst-error_%s.RDS", id))
}