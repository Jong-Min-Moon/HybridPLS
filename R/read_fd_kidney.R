################# Read functions for kidney data ###################

read_fd_kidney <- function(
    training_ratio,
    n_basis,
    normalize = list("curve" = TRUE, "scalar"= TRUE, "between" = TRUE)
    ){

  kidney_value <- extract_value_kidney()
  kidney_value$y <- response_mean_diagnosis(kidney_value$y)
  kidney_value <- preprocess_reno(kidney_value) #preprocessing, only for Emory kidney data
  kidney_value_split <- train_test_split(kidney_value, 0.3)
  kidney_value_train <- kidney_value_split$train
  kidney_value_test <- kidney_value_split$test

  kidney_predictor_train <- create_hybrid_predictors_kidney(kidney_value_train, 15)
  kidney_predictor_test <- create_hybrid_predictors_kidney(kidney_predictor_test, 15)

  #curve normalization
  if(normalize$"curve" == TRUE){
    kidney_predictor_curvenormalized <- curve_normalize_train_test(kidney_predictor_train, kidney_predictor_test)
    kidney_predictor_train <- kidney_predictor_curvenormalized$train
    kidney_predictor_test <- kidney_predictor_curvenormalized$test
  }

  # scalar normalization
  if(normalize$"scalar" == TRUE){
    kidney_predictor_scalarormalized <- scalar_normalize_train_test(kidney_predictor_train, kidney_predictor_test)
    kidney_predictor_train <- kidney_predictor_scalarormalized$train
    kidney_predictor_test <- kidney_predictor_scalarormalized$test
  }

  if(normalize$"between" == TRUE){
    kidney_predictor_btwnormalized <- btwn_normalize_train_test(kidney_predictor_train, kidney_predictor_test)
    kidney_predictor_train <- kidney_predictor_btwnormalized$train
    kidney_predictor_test <- kidney_predictor_btwnormalized$test
  }
  # btwn normalization


  result <- list()
  cat("Training data:\n")
  result$"y_train" = kidney_value_train$y
  result$"W_train" = kidney_predictor_train

  cat("Test data:\n")
  result$"y_test" = kidney_value_test$y
  result$"W_test" = kidney_predictor_test

  return(result)
}
