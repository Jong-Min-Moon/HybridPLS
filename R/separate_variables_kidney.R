separate_variables_kidney <- function(
    response_function,
    training_ratio,
    normalize = list("curve" = TRUE, "scalar"= TRUE)
    ){

  # 1. Prelim
  ## pre-specified values
  K = 2 # number of functional predictor variables
  p = 15 #number of scalar predictor variables
  d_base <- 59 # number of observed points of one baseline renogram curve
  d_post <- 40 # number of observed points of one post-furosemide renogram curve

  list_ID <- unique(kidney$ID) #unique identifier of the kidneys.
  n_sample <- length(list_ID)

  ## initialize data-saving objects
  reno_base <- array(numeric(), c(n_sample, d_base)) #baseline renogram curves
  reno_post <- array(numeric(), c(n_sample, d_post)) #post-furosemide  renogram curves
  scalar_predictors <- matrix(nrow = n_sample, ncol = p)
  y <- rep(NA, n_sample)



  # 2. separate variables


  for (i in 1:n_sample){
    ID <- list_ID[i] # one object
    data_observation <- kidney[kidney$ID == ID,] #data of one object
    data_base <- data_observation[data_observation$Study == "Baseline", ]
    data_post <- data_observation[data_observation$Study != "Baseline", ]
    data_scalar <- data_observation[1, c(12, 16:29)] #scalar predictor variables. 12 = age. 16-29

    # save values in the matrix
    reno_base[i, ] <- (data_base$renogram_value)
    reno_post[i, ] <- (data_post$renogram_value)
    scalar_predictors[i, ] <- as.numeric(data_scalar)
    y[i] <- response_function(data_observation) #response function is an input variable, so we can change this
  }

  training_idx <- sample(n_sample, floor(training_ratio * n_sample))
  y_test <-  y[-training_idx]
  print(y_test)
  n_samples_test <- length(y_test)


  # 3. save some info
  ## observed function evaluation points rescaled into [0,1], same for all i
  argvals_base <- data_base$Time_Interval_Stamp / max(data_base$Time_Interval_Stamp)
  argvals_post <- data_post$Time_Interval_Stamp / max(data_post$Time_Interval_Stamp)

  #preserve scalar pridictor variable names
  colnames(scalar_predictors) <- colnames(data_scalar)



  # 4. Data processing
  ## procedure follows the pdf file "Data Processing.pdf"

  ## 4.1. Normalize the renogram curves, for both of training and test dataset
  if (normalize$curve){
    for (i in 1:n_sample){
      max_base <- max(reno_base[i, ])
      reno_base[i, ] <- reno_base[i, ] / max_base
      reno_post[i, ] <- reno_post[i, ] / max_base
      }
    }


  ## 4.2. Center the normalized renogram curves:
  ### training dataset:
  reno_base_train <- reno_base[training_idx,]
  reno_base_centered <- scale(reno_base_train, scale = FALSE)
  reno_base_mean <- matrix(attr(reno_base_centered, "scaled:center"), nrow=1)

  reno_post_train <- reno_post[training_idx,]
  reno_post_centered <- scale(reno_post_train, scale = FALSE)
  reno_post_mean <- matrix(attr(reno_post_centered, "scaled:center"), nrow=1)

  ### test dataset, using the training mean
  reno_base_test <- reno_base[-training_idx,]- ( matrix(rep(1, n_samples_test), ncol =1) %*% reno_base_mean )
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
