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


  if (normalize$between){
    cat("Scale the scalar predictors so that the variability between the functional and scalar predictors are comparable:")
    X_1_star <- result$W_train@predictor_functional_list[[1]]
    X_2_star <- result$W_train@predictor_functional_list[[2]]
    Z_s <- result$W_train@Z
    omega <- (get_norm_sqrd(X_1_star) + get_norm_sqrd(X_2_star))/ sum((Z_s)^2)

    result$W_train@Z <- sqrt(omega) * (result$W_train@Z)
    result$W_test@Z <- sqrt(omega) * (result$W_test@Z)
    }

    return(result)
}
