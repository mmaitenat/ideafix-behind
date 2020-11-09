trainRF <- function(X, Y, seed.id) {
  library(h2o)
  Sys.unsetenv("http_proxy")
  seed.id <- as.numeric(seed.id)
  
  ## Create an H2O cloud 
  h2o.init(
    port = 54321,
    nthreads = -1,            ## -1: use all available threads
    max_mem_size = "8G")    ## specify the memory size for the H2O cloud
  # h2o.removeAll() # Clean slate - just in case the cluster was already running
  
  ## For h2o, we need X and Y matrices to be together in the same matrix
  train <- bind_cols(X, Y = Y)
  
  # Type transform
  train.h2o <- as.h2o(train)
  
  rate.per.class.list <- c(0.6320000291, 0.6320000291)
  
  ## run our first predictive model
  model <- h2o.randomForest(
    training_frame = train.h2o,    ## the H2O frame for training
    x = 1:(ncol(train.h2o)-1),     ## the predictor columns, by column index
    y = ncol(train.h2o),           ## the target index (what we are predicting)
    nfolds = 10, 
    fold_assignment = "Stratified",
    score_each_iteration = TRUE, ## Predict against training and validation for
    keep_cross_validation_predictions = TRUE,
    categorical_encoding = "AUTO",
    model_id = "KCV",    ## name the model in H2O
    sample_rate_per_class = rate.per.class.list,  ## setting this to 1 wont do any sampling 
    ntrees = 500,
    mtries = -1, # If set to -1, defaults to sqrtp for classification and p/3 for regression (where p is the # of predictors Defaults to -1.
    max_depth = 20, # same as nodesize
    stopping_rounds = 0,  
    seed = seed.id)   ## Set the random seed so that this can be reproduced
  
  saveRDS(model, file = sprintf("RF_model_%s.RDS", seed.id))
}
