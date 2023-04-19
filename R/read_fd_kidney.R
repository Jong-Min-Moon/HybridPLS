read_fd_kidney <- function(
    response_function,
    training_ratio,
    n_basis,
    scale = TRUE){
  kidney_variables <- .separate_variables_kidney(
    response_function,
    training_ratio,
    scale)

  return(list(
    #training data
    y_train = kidney_variables$training_set$y,
    W_train = .turn_into_hybrid_kidney(kidney_variables$training_set, n_basis),
    W_train_mean = .turn_into_hybrid_kidney(kidney_variables$training_mean, n_basis),

    #test data
    y_test = kidney_variables$test_set$y,
    W_test = .turn_into_hybrid_kidney(kidney_variables$test_set, n_basis)
  ))
}
