trainNN <- function(X, Y, id) {
  set.seed(id)
  
  tunegrid <- data.frame(.size = 10, 
                         .decay = 0)
  
  fitControl.noCV <- trainControl(
    method = "none",
    summaryFunction = twoClassSummary,
    classProbs = TRUE,
    verboseIter = TRUE)
  
  resubst.model <- caret::train(x = X, 
                                y = Y, 
                                method = "nnet", 
                                tuneGrid = tunegrid,
                                trControl = fitControl.noCV, 
                                preProcess = c("center", "scale"),
                                metric = "ROC",
                                na.action = na.pass)
  
  resubst.pred <- predict(resubst.model, X, type= "prob")
  resubst.out <- pROC::roc(Y, resubst.pred[, 1], plot = FALSE)
  
  # save
  saveRDS(resubst.out, file = sprintf("NN_resubst-error_%s.RDS", id))
}