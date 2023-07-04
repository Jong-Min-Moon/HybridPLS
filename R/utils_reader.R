pick_value <- function(value_object, idx){
  #value_object is a matrix, not yet turnt into a functional data object
  value_object_picked <- value_object
  n_x_functional <- length(value_object$x_functional)
  value_object_picked$y <-value_object$y[idx, ]
  value_object_picked$x_scalar <-value_object$x_scalar[idx, ]
  for (i in 1:n_x_functional){
    value_object_picked$x_functional[[i]]$value <- (value_object$x_functional[[i]]$value)[idx, ]
    value_object_picked$x_functional[[i]]$timestamp <- (value_object$x_functional[[i]]$timestamp)[idx, ]
  }


  return(value_object_picked)
}

train_test_split <- function(value_object, test_ratio){
  #value_object is a matrix, not yet turnt into a functional data object
  n_sample <- nrow(value_object$x_scalar)
  n_test <- floor(test_ratio * n_sample)
  test_idx <- sample(n_sample, n_test)
  print(test_idx)
  value_object_train <- pick_value(value_object, -test_idx)
  value_object_test <- pick_value(value_object, test_idx)

  return(
    list(
      train = value_object_train,
      test = value_object_test
    )
  )
}
