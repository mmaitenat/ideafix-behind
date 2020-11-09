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
  
  fitControl <- trainControl(
    method = "cv",
    p = 0.8,
    number = 10,
    summaryFunction = twoClassSummary,
    classProbs = TRUE,
    verboseIter = TRUE)

  
  model <- caret::train(x = as.matrix(X), 
                        y = Y, 
                        method = "xgbTree", 
                        tuneGrid = tunegrid,
                        trControl = fitControl, 
                        preProcess = c("center", "scale"),
                        metric = "ROC")
  
  # save
  saveRDS(model$resample, file = sprintf("xgbTree_model_%s.RDS", id))
}