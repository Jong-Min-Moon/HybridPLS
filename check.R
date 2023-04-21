

par(mfrow=c(2,2))
y_limit = c(0.04, 0.15)
min_mse <- rep(NA, 4)
# scale=True, seed = 1 ----------------------------------------------------



# 1. No scaling

seed = 1
set.seed(seed) # for training/test split

#set the number of max iterations
L_max = 20
lambda = c(1e-6,5*1e-6)

kidney_obj <- read_fd_kidney(
  response_function = response_mean_diagnosis,
  training_ratio = 0.7,
  n_basis = 20,
  normalize = list("curve" = FALSE, "scalar"= FALSE, "between" = FALSE)
  )


# W is already centered
W_train_centered <- kidney_obj$W_train
y_train <- kidney_obj$y_train

# minmax scale and logit transform and centering
y_train_min <- min(y_train)
y_train_max <- max(y_train)

y_train_logit <- response_transform_min_max_logit(y_train,
                                                  y_train_min,
                                                  y_train_max)
y_train_logit_mean <- mean(y_train_logit) #centering again, for PLS algorithm
y_train_logit_centered <- y_train_logit - y_train_logit_mean


# test data set
W_test_centered <- kidney_obj$W_test
y_test <- kidney_obj$y_test


#MSE trend

L_trend_0_logit <- rep(NA, L_max)
for(L in 1:L_max){
  pls_model_transform <- nipals_pen_hybrid(
    W_train_centered, y_train_logit_centered, L, lambda, 0)

  #inverse transform
  y_pred_pls <- predict_test(pls_model_transform, W_test_centered)
  y_pred_pls <- y_pred_pls + y_train_logit_mean
  y_pred_pls <- exp(y_pred_pls)/( 1 + exp(y_pred_pls))
  y_pred_pls <- y_pred_pls * (y_train_max- y_train_min)
  y_pred_pls <- y_pred_pls + y_train_min

  L_trend_0_logit[L] <- mse(y_pred_pls, y_test)
}

plot(L_trend_0_logit[1:L_max], main ="no scaling", ylab = "MSE", xlab = "number of iteration", pch=".", ylim=y_limit)
lines(L_trend_0_logit)
min_mse[1] <- min(L_trend_0_logit, omit.NA = TRUE)



######
# 2. curve scaling

seed = 1
set.seed(seed) # for training/test split

#set the number of max iterations
L_max = 20
lambda = c(1e-6,5*1e-6)

kidney_obj <- read_fd_kidney(
  response_function = response_mean_diagnosis,
  training_ratio = 0.7,
  n_basis = 20,
  normalize = list("curve" = TRUE, "scalar"= FALSE, "between" = FALSE)
)


# W is already centered
W_train_centered <- kidney_obj$W_train
y_train <- kidney_obj$y_train

# minmax scale and logit transform and centering
y_train_min <- min(y_train)
y_train_max <- max(y_train)

y_train_logit <- response_transform_min_max_logit(y_train,
                                                  y_train_min,
                                                  y_train_max)
y_train_logit_mean <- mean(y_train_logit) #centering again, for PLS algorithm
y_train_logit_centered <- y_train_logit - y_train_logit_mean


# test data set
W_test_centered <- kidney_obj$W_test
y_test <- kidney_obj$y_test


#MSE trend

L_trend_0_logit <- rep(NA, L_max)
for(L in 1:L_max){
  pls_model_transform <- nipals_pen_hybrid(
    W_train_centered, y_train_logit_centered, L, lambda, 0)

  #inverse transform
  y_pred_pls <- predict_test(pls_model_transform, W_test_centered)
  y_pred_pls <- y_pred_pls + y_train_logit_mean
  y_pred_pls <- exp(y_pred_pls)/( 1 + exp(y_pred_pls))
  y_pred_pls <- y_pred_pls * (y_train_max- y_train_min)
  y_pred_pls <- y_pred_pls + y_train_min

  L_trend_0_logit[L] <- mse(y_pred_pls, y_test)
}

plot(L_trend_0_logit[1:L_max], main ="curve scaling", ylab = "MSE", xlab = "number of iteration", pch=".", ylim=y_limit)
lines(L_trend_0_logit)
min_mse[2] <- min(L_trend_0_logit, omit.NA = TRUE)



######
# 3. curve and scalar scaling

seed = 1
set.seed(seed) # for training/test split

#set the number of max iterations
L_max = 20
lambda = c(1e-6,5*1e-6)

kidney_obj <- read_fd_kidney(
  response_function = response_mean_diagnosis,
  training_ratio = 0.7,
  n_basis = 20,
  normalize = list("curve" = TRUE, "scalar"= TRUE, "between" = FALSE)
)


# W is already centered
W_train_centered <- kidney_obj$W_train
y_train <- kidney_obj$y_train

# minmax scale and logit transform and centering
y_train_min <- min(y_train)
y_train_max <- max(y_train)

y_train_logit <- response_transform_min_max_logit(y_train,
                                                  y_train_min,
                                                  y_train_max)
y_train_logit_mean <- mean(y_train_logit) #centering again, for PLS algorithm
y_train_logit_centered <- y_train_logit - y_train_logit_mean


# test data set
W_test_centered <- kidney_obj$W_test
y_test <- kidney_obj$y_test


#MSE trend

L_trend_0_logit <- rep(NA, L_max)
for(L in 1:L_max){
  pls_model_transform <- nipals_pen_hybrid(
    W_train_centered, y_train_logit_centered, L, lambda, 0)

  #inverse transform
  y_pred_pls <- predict_test(pls_model_transform, W_test_centered)
  y_pred_pls <- y_pred_pls + y_train_logit_mean
  y_pred_pls <- exp(y_pred_pls)/( 1 + exp(y_pred_pls))
  y_pred_pls <- y_pred_pls * (y_train_max- y_train_min)
  y_pred_pls <- y_pred_pls + y_train_min

  L_trend_0_logit[L] <- mse(y_pred_pls, y_test)
}

plot(L_trend_0_logit[1:L_max], main ="curve + scalar scaling", ylab = "MSE", xlab = "number of iteration", pch=".", ylim=y_limit)
lines(L_trend_0_logit)
min_mse[3] <- min(L_trend_0_logit, omit.NA = TRUE)



######
# 4. curve and scalar and between scaling

seed = 1
set.seed(seed) # for training/test split

#set the number of max iterations
L_max = 20
lambda = c(1e-6,5*1e-6)

kidney_obj <- read_fd_kidney(
  response_function = response_mean_diagnosis,
  training_ratio = 0.7,
  n_basis = 20,
  normalize = list("curve" = TRUE, "scalar"= TRUE, "between" = TRUE)
)


# W is already centered
W_train_centered <- kidney_obj$W_train
y_train <- kidney_obj$y_train

# minmax scale and logit transform and centering
y_train_min <- min(y_train)
y_train_max <- max(y_train)

y_train_logit <- response_transform_min_max_logit(y_train,
                                                  y_train_min,
                                                  y_train_max)
y_train_logit_mean <- mean(y_train_logit) #centering again, for PLS algorithm
y_train_logit_centered <- y_train_logit - y_train_logit_mean


# test data set
W_test_centered <- kidney_obj$W_test
y_test <- kidney_obj$y_test


#MSE trend

L_trend_0_logit <- rep(NA, L_max)
for(L in 1:L_max){
  pls_model_transform <- nipals_pen_hybrid(
    W_train_centered, y_train_logit_centered, L, lambda, 0)

  #inverse transform
  y_pred_pls <- predict_test(pls_model_transform, W_test_centered)
  y_pred_pls <- y_pred_pls + y_train_logit_mean
  y_pred_pls <- exp(y_pred_pls)/( 1 + exp(y_pred_pls))
  y_pred_pls <- y_pred_pls * (y_train_max- y_train_min)
  y_pred_pls <- y_pred_pls + y_train_min

  L_trend_0_logit[L] <- mse(y_pred_pls, y_test)
}

plot(L_trend_0_logit[1:L_max], main ="curve + scalar + between scaling", ylab = "MSE", xlab = "number of iteration", pch=".", ylim=y_limit)
lines(L_trend_0_logit)
min_mse[4] <- min(L_trend_0_logit, omit.NA = TRUE)

print(min_mse)
