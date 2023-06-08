








separate_variables_kidney <- function(
    response_function,
    training_ratio,
    normalize = list("curve" = TRUE, "scalar"= TRUE)
    ){

  # 1. Prelim
  ## pre-specified values. These values are dataset-specific. DO NOT CHANGE THESE VALUES.
  K = 2 # number of functional predictor variables






  ## initialize data-saving objects






  # 2. separate variables







  #preserve scalar pridictor variable names
  colnames(scalar_predictors) <- colnames(data_scalar)



  # 4. Data processing
  ## procedure follows the pdf file "Data Processing.pdf"

  ## 4.1. Pre-process (only for our data) the renogram curves by dividing:
  ### a) each baseline renogram curve by its maximum; and
  ### b) each post-furosemide (diuretic) renogram curve by the maximum of the baseline renogram curve
  reno_preprocessed <- preprocess_reno(reno_base, reno_post)
  reno_base <- reno_preprocessed$base
  reno_post <- reno_preprocessed$post



  #Normalize the renogram curves, for both of training and test dataset
  ## 4.2. Center the normalized renogram curves:
  ### training dataset:
  reno_base_train <- reno_base[training_idx,]
  reno_base_test <- reno_base[-training_idx,]

  - ( matrix(rep(1, n_samples_test), ncol =1) %*% reno_base_mean )
  normalize_curve_2d <- function(curve_train, curve_test){
    reno_base_mean <- apply(curve_train, 2, mean)
    reno_base_sd <- apply(curve_train, 2, sd)
  }




  numerator <-
  reno_base_centered <- scale(reno_base_train, scale = FALSE)


  reno_post_train <- reno_post[training_idx,]
  reno_post_centered <- scale(reno_post_train, scale = FALSE)
  reno_post_mean <- matrix(attr(reno_post_centered, "scaled:center"), nrow=1)

  ### test dataset, using the training mean

  reno_post_test <- reno_post[-training_idx,]- ( matrix(rep(1, n_samples_test), ncol =1) %*% reno_post_mean )


  ## 4.3. Standardize the scalar predictors
  ### training dataset:
  scalar_predictors_train <- scalar_predictors[training_idx,]

  if(normalize$scalar){
    scalar_predictors_centered <- scale(scalar_predictors_train, scale = TRUE)
    scalar_predictors_sd <- matrix(attr(scalar_predictors_centered, "scaled:scale"), nrow=1)
  }else{
    scalar_predictors_centered <- scale(scalar_predictors_train, scale = FALSE)
    scalar_predictors_sd <- matrix(rep(1,15), nrow = 1)
  }

  scalar_predictors_mean <- matrix(attr(scalar_predictors_centered, "scaled:center"), nrow=1)

  ### test dataset:
  scalar_predictors_test <- scalar_predictors[-training_idx,]
  scalar_predictors_test <- scalar_predictors_test - ( matrix(rep(1, n_samples_test), ncol =1) %*% scalar_predictors_mean )
  scalar_predictors_test <- scalar_predictors_test / ( matrix(rep(1, n_samples_test), ncol =1) %*% scalar_predictors_sd )


  ## 4.4. Scale the scalar predictors so that the variability between the functional and scalar predictors are comparable:
  ### need to calculate the functional norm, so we'll do it later



  # 5. save results




  # mean function for training
  training_mean <- list(
    reno_base_value = reno_base_mean,
    reno_post_value = reno_post_mean,
    reno_base_x = argvals_base,
    reno_post_x = argvals_post,
    scalar_predictors = scalar_predictors_mean,
    scalar_predictors_sd = scalar_predictors_sd
  )
  # training set
  training_set <- list(
    reno_base_value = reno_base_centered,
    reno_post_value = reno_post_centered,
    reno_base_x = argvals_base,
    reno_post_x = argvals_post,
    scalar_predictors = scalar_predictors_centered,
    y = y[training_idx]
  )

  # testing set
  test_set <- list(
    reno_base_value = reno_base_test,
    reno_post_value = reno_post_test,
    reno_base_x = argvals_base,
    reno_post_x = argvals_post,
    scalar_predictors = scalar_predictors_test,
    y = y_test
  )

  return(list(training_set = training_set, test_set = test_set, training_mean = training_mean))
}
