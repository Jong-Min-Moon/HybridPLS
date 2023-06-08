################# Read functions for kidney data ###################

read_fd_kidney <- function(
    response_function,
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

  kidney_predictor_curvenormalized <- curve_normalize_train_test(kidney_predictor_train, kidney_predictor_test)

  kidney_variables <- separate_variables_kidney(
    response_function,
    training_ratio,
    list("curve" = normalize$curve, "scalar" = normalize$scalar)
    )


  result <- list()
  cat("Training data:\n")
  result$"y_train" = kidney_variables$training_set$y
  result$"W_train" = turn_into_hybrid_kidney(kidney_variables$training_set, n_basis)

  cat("Training data (mean function):\n")
  result$"W_train_mean" = turn_into_hybrid_kidney(kidney_variables$training_mean, n_basis)

  cat("Test data:\n")
  result$"y_test" = kidney_variables$test_set$y
  result$"W_test" = turn_into_hybrid_kidney(kidney_variables$test_set, n_basis)



    return(result)
}
