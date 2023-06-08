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
    reno_base_timestamp[i,] <- (data_base$Time_Interval_Stamp)/max(data_base$Time_Interval_Stamp)

    data_post <- data_observation[data_observation$Study != "Baseline", ]
    reno_post_value[i, ] <- (data_post$renogram_value)
    reno_post_timestamp[i,] <- (data_post$Time_Interval_Stamp)/max(data_post$Time_Interval_Stamp)

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
