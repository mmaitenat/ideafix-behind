trainRF <- function(X.train, Y.train, X.test, Y.test, seed.id) {
  library(h2o)
  Sys.unsetenv("http_proxy")
  seed.id <- as.numeric(seed.id)
  
  ## Create an H2O cloud 
  h2o.init(
    port = 54321 + sample(1:27, 1) + (seed.id*2),
    nthreads = -1,            ## -1: use all available threads
    max_mem_size = "20G")    ## specify the memory size for the H2O cloud
  # h2o.removeAll() # Clean slate - just in case the cluster was already running
  
  ## For h2o, we need X and Y matrices to be together in the same matrix
  train <- bind_cols(X.train, Y = Y.train)
  test <- bind_cols(X.test, Y = Y.test)
  
  # Type transform
  train.h2o <- as.h2o(train)
  test.h2o <- as.h2o(test)
  
  rate.per.class.list <- c(0.6320000291, 0.6320000291)
  
  ## run our first predictive model
  model <- h2o.randomForest(
    training_frame = train.h2o,    ## the H2O frame for training
    x = 1:(ncol(train.h2o)-1),     ## the predictor columns, by column index
    y = ncol(train.h2o),           ## the target index (what we are predicting)
    categorical_encoding = "AUTO",
    model_id = "TEST",    ## name the model in H2O
    sample_rate_per_class = rate.per.class.list,  ## setting this to 1 wont do any sampling 
    ntrees = 500,
    mtries = -1, # If set to -1, defaults to sqrtp for classification and p/3 for regression (where p is the # of predictors Defaults to -1.
    max_depth = 20, # same as nodesize
    stopping_rounds = 0,  
    seed = seed.id)   ## Set the random seed so that this can be reproduced
  
  # la prediccion tambien se podria haber hecho con h2o.performance()
  test.pred <- h2o.predict(object = model, 
                           newdata = test.h2o)
  test.out <- pROC::roc(Y.test, as.data.frame(test.pred)$X1, plot = FALSE)
  
  saveRDS(test.out, file = sprintf("RF_test-error_%s.RDS", seed.id))
}
