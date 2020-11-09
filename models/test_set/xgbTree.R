trainxgbTree <- function(X.train, Y.train, X.test, Y.test, id) {
  set.seed(id)
  
  fitControl.noCV <- trainControl(
    method = "none",
    summaryFunction = twoClassSummary,
    classProbs = TRUE,
    verboseIter = TRUE)
  
  # these are xgbtree's default parameters
  tunegrid <- expand.grid(nrounds = 100,
                          max_depth = 6,
                          eta = 0.3,
                          gamma = 0,
                          colsample_bytree = 1,
                          min_child_weight = 1,
                          subsample = 1)
  
  model <- caret::train(x = as.matrix(X.train), 
                        y = Y.train, 
                        method = "xgbTree", 
                        tuneGrid = tunegrid,
                        trControl = fitControl.noCV, 
                        preProcess = c("center", "scale"),
                        metric = "ROC")
  
  test.pred <- predict(model, as.matrix(X.test), type= "prob")
  test.out <- pROC::roc(Y.test, test.pred[, 1], plot = FALSE)
  
  # save
  saveRDS(test.out, file = sprintf("xgbTree_test-error_%s.RDS", id))
}
