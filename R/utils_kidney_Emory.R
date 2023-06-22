# exlusive for Emory kidney data
#output:
## 1. reno_base, 2. reno_post : matrix
### stores functional data
### - nrow = number of samples
### - ncol = number of observed function value points
## 3. scalar_predictors: matrix
### stores scalar data
### - nrow = number of samples
### - ncol = number of scalar predictors
## 4. y : matrix (vertical vector)
### stores response data
### - nrow = number of samples
### - ncol = dimension of the response variable
## 5. time_interval_stamp
extract_value_kidney <- function(){
  list_ID <- unique(kidney$ID) #unique identifier of the kidneys
  n_sample <- length(list_ID)

  reno_base_value <- reno_base_timestamp <- array(numeric(), c(n_sample, 59)) #baseline renogram curves
  reno_post_value <- reno_post_timestamp <- array(numeric(), c(n_sample, 40)) #post-furosemide renogram curves

  scalar_predictors <- array(numeric(), c(n_sample, 15))
  response <- array(numeric(), c(n_sample, 3))

  for (i in 1:n_sample){
    ID <- list_ID[i] # one object
    data_observation <- kidney[kidney$ID == ID, ] #data of one object

    data_base <- data_observation[data_observation$Study == "Baseline", ]
    reno_base_value[i, ] <- (data_base$renogram_value)
    reno_base_timestamp[i,] <- (data_base$Time_Interval_Stamp)/max(data_base$Time_Interval_Stamp) #timestamp rescaled to [0,1]

    data_post <- data_observation[data_observation$Study != "Baseline", ]
    reno_post_value[i, ] <- (data_post$renogram_value)
    reno_post_timestamp[i,] <- (data_post$Time_Interval_Stamp)/max(data_post$Time_Interval_Stamp) #timestamp rescaled to [0,1]

    scalar_predictors[i, ] <- as.numeric(data_base[1, c(12, 16:29)]) #scalar predictor variables. 12 = age. 16-29
    response[i,] <- as.numeric(data_base[1, c(13, 14, 15)])
    }
  return(
    list(
      "x_functional" = list(
        "first" = list("value" = reno_base_value, "timestamp" = reno_base_timestamp),
        "second" = list("value" = reno_post_value, "timestamp" = reno_post_timestamp)
        ),
      "x_scalar" = scalar_predictors,
      "y" = response
    )
  )
}

response_mean_diagnosis <- function(y){
  return(
    matrix(apply(y, 1, mean))
    )
  }

preprocess_reno <- function(kidney_value){
  # Pre-process (only for our data) the renogram curves by dividing:
  ### a) each baseline renogram curve by its maximum; and
  ### b) each post-furosemide (diuretic) renogram curve by the maximum of the baseline renogram curve
  new_kidney_value <- kidney_value
  reno_base <- kidney_value$x_functional$first$value
  reno_post <- kidney_value$x_functional$second$value
  n_sample <- nrow(reno_base)


  for (i in 1:n_sample){
    max_base <- max(reno_base[i, ])
    new_kidney_value$x_functional$first$value[i, ] <- reno_base[i, ] / max_base # a)
    new_kidney_value$x_functional$second$value[i, ] <- reno_post[i, ] / max_base # b)
  }
  return(new_kidney_value)
}

################# Read functions for kidney data ###################

read_fd_kidney <- function(
    test_ratio,
    n_basis,
    normalize = list("curve" = TRUE, "scalar"= TRUE, "between" = TRUE)
){

  kidney_value <- extract_value_kidney()
  kidney_value$y <- response_mean_diagnosis(kidney_value$y)
  kidney_value <- preprocess_reno(kidney_value) #preprocessing, only for Emory kidney data
  kidney_value_split <- train_test_split(kidney_value, test_ratio)
  kidney_value_train <- kidney_value_split$train
  kidney_value_test <- kidney_value_split$test

  kidney_predictor_train <- create_hybrid_predictors_kidney(kidney_value_train, n_basis)
  kidney_predictor_test <- create_hybrid_predictors_kidney(kidney_value_test, n_basis)

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

create_hybrid_predictors_kidney <- function(value_object, n_basis){
  Z <- value_object$x_scalar

  cat(paste("\t For base curves, "))
  split_fit_base <- fit_spine_2d(
    argvals = value_object$x_functional$first$timestamp[1,],
    evals = value_object$x_functional$first$value,
    n_basis = n_basis)

  cat(paste("\t For post curves, "))
  predictor_functional_1 <- create_predictor_functional(
    split_fit_base$C,
    split_fit_base$J,
    split_fit_base$J_dotdot,
    value_object$x_functional$first$timestamp,
    value_object$x_functional$first$value
  )

  split_fit_post <- fit_spine_2d(
    argvals = value_object$x_functional$second$timestamp[1,],
    evals = value_object$x_functional$second$value,
    n_basis = n_basis)

  predictor_functional_2 <- create_predictor_functional(
    split_fit_post$C,
    split_fit_post$J,
    split_fit_post$J_dotdot,
    value_object$x_functional$second$timestamp,
    value_object$x_functional$second$value
  )

  predictor_object <- new("predictor_hybrid",
                          Z = Z,
                          predictor_functional_list = list(
                            predictor_functional_1,
                            predictor_functional_2
                          ),
                          n_predictor_functional = 2,
                          n_sample = nrow(Z)
  )
  return(predictor_object)
}


