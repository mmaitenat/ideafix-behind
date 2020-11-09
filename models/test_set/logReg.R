trainLogReg <- function(X.train, Y.train, X.test, Y.test, id) {
  set.seed(id)
  
  fitControl.noCV <- trainControl(
    method = "none",
    summaryFunction = twoClassSummary,
    classProbs = TRUE,
    verboseIter = TRUE)
  
  model <- caret::train(x = X.train, 
                        y = Y.train, 
                        method = "glm", 
                        family = "binomial",
                        trControl = fitControl.noCV, 
                        preProcess = c("center", "scale"),
                        metric = "ROC",
                        na.action = na.pass)
  
  test.pred <- predict(model, X.test, type = "prob")
  test.out <- pROC::roc(Y.test, test.pred[, 1], plot = FALSE)
  
  # save
  saveRDS(test.out, file = sprintf("logReg_test-error_%s.RDS", id))
}
