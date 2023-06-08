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
  functional_predictor_normalized <- subtr_broadcast(functional_predictor, mean) #numerator
  functional_predictor_normalized <- scalar_mul(functional_predictor_normalized, 1/deno)
  return(functional_predictor_normalized)
}


