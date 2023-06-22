curve_normalize_train_test <- function(predictor_train, predictor_test){
  n_sample <- nrow(predictor_train@Z)
  predictor_train_normalized <- predictor_train
  predictor_test_normalized <- predictor_test

  for (k in 1:(predictor_train@n_predictor_functional)){ # for kth functional predictor
    mean_train <- get_mean(predictor_train@predictor_functional_list[[k]]) # train mean function
    centered_train <-
      subtr_broadcast(predictor_train@predictor_functional_list[[k]], mean_train)
    deno_train <- sqrt(get_sum_of_norm_sqrd(centered_train)/(n_sample-1))

    predictor_train_normalized@predictor_functional_list[[k]] <-
      curve_normalize(
        predictor_train_normalized@predictor_functional_list[[k]],
        mean_train,
        deno_train
      )

    predictor_test_normalized@predictor_functional_list[[k]] <-
      curve_normalize(
        predictor_test_normalized@predictor_functional_list[[k]],
        mean_train,
        deno_train
      )
  }
  return(
    list(
      train = predictor_train_normalized,
      test = predictor_test_normalized
    )
  )
}

curve_normalize <- function(functional_predictor, mean, deno){
  # operation also takes care of the original_X
  functional_predictor_normalized <- subtr_broadcast(functional_predictor, mean) #numerator
  functional_predictor_normalized <- scalar_mul(functional_predictor_normalized, 1/deno)
  return(functional_predictor_normalized)
}

scalar_normalize_train_test <- function(predictor_train, predictor_test){
  predictor_train_normalized <- predictor_train
  predictor_test_normalized <- predictor_test

  mean_train <- get_mean(predictor_train@Z) # train mean function
  sd_train <- get_sd(predictor_train@Z)

  predictor_train_normalized@Z <-
    scalar_normalize(
      predictor_train@Z,
      mean_train,
      sd_train
      )

  predictor_test_normalized@Z <-
    scalar_normalize(
      predictor_test@Z,
      mean_train,
      sd_train
    )

  return(
    list(
      train = predictor_train_normalized,
      test = predictor_test_normalized
    )
  )
}

scalar_normalize <- function(scalar_predictor, mean, deno){
  scalar_predictor_normalized <- subtr_broadcast(scalar_predictor, mean) #numerator
  deno_mat <- matrix(rep(1,nrow(scalar_predictor))) %*% deno
  scalar_predictor_normalized <- scalar_predictor_normalized / deno_mat
  return(scalar_predictor_normalized)
}

btwn_normalize_train_test <- function(predictor_train, predictor_test){
  cat("Scale the scalar predictors so that the variability between the functional and scalar predictors are comparable:")
  Z_s <- predictor_train@Z
  predictor_train_normalized <- predictor_train
  predictor_test_normalized <- predictor_test

  omega <- 0
  for (k in 1:(predictor_train@n_predictor_functional)){ # for kth functional predictor
    omega <- omega + get_sum_of_norm_sqrd(predictor_train@predictor_functional_list[[k]])
  }
    omega <- omega / sum((Z_s)^2)

    predictor_train_normalized@Z <- sqrt(omega) * (predictor_train@Z)
    predictor_test_normalized@Z <- sqrt(omega) * (predictor_test@Z)

    return(
      list(
        train = predictor_train_normalized,
        test = predictor_test_normalized
      )
    )
}
