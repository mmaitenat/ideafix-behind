trainLogReg <- function(X, Y, id) {
  set.seed(id)
  
  folds <- 10
  cvIndex <- createFolds(factor(Y), folds, returnTrain = T)
  fitControl <- trainControl(
    method = "cv",
    index = cvIndex,
    number = folds,
    summaryFunction = twoClassSummary,
    classProbs = TRUE,
    verboseIter = TRUE)
  
  model <- caret::train(x = X, 
                        y = Y, 
                        method = "glm", 
                        family = "binomial",
                        trControl = fitControl, 
                        preProcess = c("center", "scale"),
                        metric = "ROC",
                        na.action = na.pass)
  
  # save
  saveRDS(model$resample, file = sprintf("logReg_model_%s.RDS", id))

}