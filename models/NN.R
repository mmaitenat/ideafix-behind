trainNN <- function(X, Y, id) {
  set.seed(id)
  
  tunegrid <- data.frame(.size = 10, 
                         .decay = 0)
  
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
                        method = "nnet", 
                        tuneGrid = tunegrid,
                        trControl = fitControl, 
                        preProcess = c("center", "scale"),
                        metric = "ROC",
                        na.action = na.pass)
  
  # save
  saveRDS(model$resample, file = sprintf("NN_model_%s.RDS", id))
}