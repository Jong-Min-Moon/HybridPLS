

par(mfrow=c(1,1))
# scale=True, seed = 1 ----------------------------------------------------

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

plot(L_trend_0_logit[1:L_max], main =
       paste("kappa=1e-6 and 5*1e-6, tau = 0, L_max = 20 ",
             "\n MSE for PLS with logit transform\n", "seed = ", seed, " curves scaled"
       ), ylab = "MSE", xlab = "number of iteration", pch=".", ylim=c(0.04, 0.11))
lines(L_trend_0_logit)
print(min(L_trend_0_logit, omit.NA = TRUE))



# scale=True, seed = 10 ----------------------------------------------------

seed = 10
set.seed(seed) # for training/test split

#set the number of max iterations
L_max = 20


kidney_obj <- read_fd_kidney_obj(
  dir = "~/GitHub/fingernail/data/renogram_data.csv",
  response_function = response_mean_diagnosis,
  training_ratio = 0.7,
  n_basis = 20,
  scale = TRUE
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
    W_train_centered, y_train_logit_centered, L, c(1e-6,5*1e-6), 0)

  #inverse transform
  y_pred_pls <- predict_test(pls_model_transform, W_test_centered)
  y_pred_pls <- y_pred_pls + y_train_logit_mean
  y_pred_pls <- exp(y_pred_pls)/( 1 + exp(y_pred_pls))
  y_pred_pls <- y_pred_pls * (y_train_max- y_train_min)
  y_pred_pls <- y_pred_pls + y_train_min

  L_trend_0_logit[L] <- mse(y_pred_pls, y_test)
}

plot(L_trend_0_logit[1:L_max], main =
       paste("kappa=1e-6 and 5*1e-6, tau = 0, L_max = 20 ",
             "\n MSE for PLS with logit transform\n", "seed = ", seed, " curves scaled"
       ), ylab = "MSE", xlab = "number of iteration", pch=".", ylim=c(0.04, 0.16))
lines(L_trend_0_logit)
print(min(L_trend_0_logit, omit.NA = TRUE))



# scale= True, seed = 100 -------------------------------------------------
seed = 100
set.seed(seed) # for training/test split

#set the number of max iterations
L_max = 20


kidney_obj <- read_fd_kidney(
  dir = "~/GitHub/fingernail/data/renogram_data.csv",
  response_function = response_mean_diagnosis,
  training_ratio = 0.7,
  n_basis = 20,
  scale = TRUE
)

# W is already centered
W_train_centered <- kidney$W_train
y_train <- kidney$y_train

# minmax scale and logit transform and centering
y_train_min <- min(y_train)
y_train_max <- max(y_train)

y_train_logit <- response_transform_min_max_logit(y_train,
                                                  y_train_min,
                                                  y_train_max)
y_train_logit_mean <- mean(y_train_logit) #centering again, for PLS algorithm
y_train_logit_centered <- y_train_logit - y_train_logit_mean


# test data set
W_test_centered <- kidney$W_test
y_test <- kidney$y_test


#MSE trend

L_trend_0_logit <- rep(NA, L_max)
for(L in 1:L_max){
  pls_model_transform <- nipals_pen_hybrid(
    W_train_centered, y_train_logit_centered, L, c(1e-6,5*1e-6), 0)

  #inverse transform
  y_pred_pls <- predict_test(pls_model_transform, W_test_centered)
  y_pred_pls <- y_pred_pls + y_train_logit_mean
  y_pred_pls <- exp(y_pred_pls)/( 1 + exp(y_pred_pls))
  y_pred_pls <- y_pred_pls * (y_train_max- y_train_min)
  y_pred_pls <- y_pred_pls + y_train_min

  L_trend_0_logit[L] <- mse(y_pred_pls, y_test)
}

plot(L_trend_0_logit[1:L_max], main =
       paste("kappa=1e-6 and 5*1e-6, tau = 0, L_max = 20 ",
             "\n MSE for PLS with logit transform\n", "seed = ", seed, " curves scaled"
       ), ylab = "MSE", xlab = "number of iteration", pch=".", ylim=c(0.04, 0.16))
lines(L_trend_0_logit)
print(min(L_trend_0_logit, omit.NA = TRUE))

source("class_hybrid_predictors.R")
source("functions.R")



# scale= True, seed = 1000 -------------------------------------------------
seed = 1000
set.seed(seed) # for training/test split

#set the number of max iterations
L_max = 20


kidney <- read_fd_kidney(
  dir = "~/GitHub/fingernail/data/renogram_data.csv",
  response_function = response_mean_diagnosis,
  training_ratio = 0.7,
  n_basis = 20,
  scale = TRUE
)

# W is already centered
W_train_centered <- kidney$W_train
y_train <- kidney$y_train

# minmax scale and logit transform and centering
y_train_min <- min(y_train)
y_train_max <- max(y_train)

y_train_logit <- response_transform_min_max_logit(y_train,
                                                  y_train_min,
                                                  y_train_max)
y_train_logit_mean <- mean(y_train_logit) #centering again, for PLS algorithm
y_train_logit_centered <- y_train_logit - y_train_logit_mean


# test data set
W_test_centered <- kidney$W_test
y_test <- kidney$y_test


#MSE trend

L_trend_0_logit <- rep(NA, L_max)
for(L in 1:L_max){
  pls_model_transform <- nipals_pen_hybrid(
    W_train_centered, y_train_logit_centered, L, c(1e-6,5*1e-6), 0)

  #inverse transform
  y_pred_pls <- predict_test(pls_model_transform, W_test_centered)
  y_pred_pls <- y_pred_pls + y_train_logit_mean
  y_pred_pls <- exp(y_pred_pls)/( 1 + exp(y_pred_pls))
  y_pred_pls <- y_pred_pls * (y_train_max- y_train_min)
  y_pred_pls <- y_pred_pls + y_train_min

  L_trend_0_logit[L] <- mse(y_pred_pls, y_test)
}

plot(L_trend_0_logit[1:L_max], main =
       paste("kappa=1e-6 and 5*1e-6, tau = 0, L_max = 20 ",
             "\n MSE for PLS with logit transform\n", "seed = ", seed, " curves scaled"
       ), ylab = "MSE", xlab = "number of iteration", pch=".", ylim=c(0.04, 0.16))
lines(L_trend_0_logit)
print(min(L_trend_0_logit, omit.NA = TRUE))

source("class_hybrid_predictors.R")
source("functions.R")

# scale=FALSE, seed = 1 ----------------------------------------------------

seed = 1
set.seed(seed) # for training/test split

#set the number of max iterations
L_max = 20


kidney <- read_fd_kidney(
  dir = "~/GitHub/fingernail/data/renogram_data.csv",
  response_function = response_mean_diagnosis,
  training_ratio = 0.7,
  n_basis = 20,
  scale = FALSE)


# W is already centered
W_train_centered <- kidney$W_train
y_train <- kidney$y_train

# minmax scale and logit transform and centering
y_train_min <- min(y_train)
y_train_max <- max(y_train)

y_train_logit <- response_transform_min_max_logit(y_train,
                                                  y_train_min,
                                                  y_train_max)
y_train_logit_mean <- mean(y_train_logit) #centering again, for PLS algorithm
y_train_logit_centered <- y_train_logit - y_train_logit_mean


# test data set
W_test_centered <- kidney$W_test
y_test <- kidney$y_test


#MSE trend

L_trend_0_logit <- rep(NA, L_max)
for(L in 1:L_max){
  pls_model_transform <- nipals_pen_hybrid(
    W_train_centered, y_train_logit_centered, L, c(1e-6,5*1e-6), 0)

  #inverse transform
  y_pred_pls <- predict_test(pls_model_transform, W_test_centered)
  y_pred_pls <- y_pred_pls + y_train_logit_mean
  y_pred_pls <- exp(y_pred_pls)/( 1 + exp(y_pred_pls))
  y_pred_pls <- y_pred_pls * (y_train_max- y_train_min)
  y_pred_pls <- y_pred_pls + y_train_min

  L_trend_0_logit[L] <- mse(y_pred_pls, y_test)
}

plot(L_trend_0_logit[1:L_max], main =
       paste("kappa=1e-6 and 5*1e-6, tau = 0, L_max = 20 ",
             "\n MSE for PLS with logit transform\n", "seed = ", seed, " curves NOT scaled"
       ), ylab = "MSE", xlab = "number of iteration", pch=".", ylim=c(0.04, 0.16))
lines(L_trend_0_logit)
print(min(L_trend_0_logit, omit.NA = TRUE))



# scale=FALSE, seed = 10 ----------------------------------------------------

seed = 10
set.seed(seed) # for training/test split

#set the number of max iterations
L_max = 20


kidney <- read_fd_kidney(
  dir = "~/GitHub/fingernail/data/renogram_data.csv",
  response_function = response_mean_diagnosis,
  training_ratio = 0.7,
  n_basis = 20,
  scale = FALSE
)

# W is already centered
W_train_centered <- kidney$W_train
y_train <- kidney$y_train

# minmax scale and logit transform and centering
y_train_min <- min(y_train)
y_train_max <- max(y_train)

y_train_logit <- response_transform_min_max_logit(y_train,
                                                  y_train_min,
                                                  y_train_max)
y_train_logit_mean <- mean(y_train_logit) #centering again, for PLS algorithm
y_train_logit_centered <- y_train_logit - y_train_logit_mean


# test data set
W_test_centered <- kidney$W_test
y_test <- kidney$y_test


#MSE trend

L_trend_0_logit <- rep(NA, L_max)
for(L in 1:L_max){
  pls_model_transform <- nipals_pen_hybrid(
    W_train_centered, y_train_logit_centered, L, c(1e-6,5*1e-6), 0)

  #inverse transform
  y_pred_pls <- predict_test(pls_model_transform, W_test_centered)
  y_pred_pls <- y_pred_pls + y_train_logit_mean
  y_pred_pls <- exp(y_pred_pls)/( 1 + exp(y_pred_pls))
  y_pred_pls <- y_pred_pls * (y_train_max- y_train_min)
  y_pred_pls <- y_pred_pls + y_train_min

  L_trend_0_logit[L] <- mse(y_pred_pls, y_test)
}

plot(L_trend_0_logit[1:L_max], main =
       paste("kappa=1e-6 and 5*1e-6, tau = 0, L_max = 20 ",
             "\n MSE for PLS with logit transform\n", "seed = ", seed, " curves NOT scaled"
       ), ylab = "MSE", xlab = "number of iteration", pch=".", ylim=c(0.04, 0.16))
lines(L_trend_0_logit)
print(min(L_trend_0_logit, omit.NA = TRUE))



# scale= True, seed = 100 -------------------------------------------------
seed = 100
set.seed(seed) # for training/test split

#set the number of max iterations
L_max = 20


kidney <- read_fd_kidney(
  dir = "~/GitHub/fingernail/data/renogram_data.csv",
  response_function = response_mean_diagnosis,
  training_ratio = 0.7,
  n_basis = 20,
  scale = FALSE
)

# W is already centered
W_train_centered <- kidney$W_train
y_train <- kidney$y_train

# minmax scale and logit transform and centering
y_train_min <- min(y_train)
y_train_max <- max(y_train)

y_train_logit <- response_transform_min_max_logit(y_train,
                                                  y_train_min,
                                                  y_train_max)
y_train_logit_mean <- mean(y_train_logit) #centering again, for PLS algorithm
y_train_logit_centered <- y_train_logit - y_train_logit_mean


# test data set
W_test_centered <- kidney$W_test
y_test <- kidney$y_test


#MSE trend

L_trend_0_logit <- rep(NA, L_max)
for(L in 1:L_max){
  pls_model_transform <- nipals_pen_hybrid(
    W_train_centered, y_train_logit_centered, L, c(1e-6,5*1e-6), 0)

  #inverse transform
  y_pred_pls <- predict_test(pls_model_transform, W_test_centered)
  y_pred_pls <- y_pred_pls + y_train_logit_mean
  y_pred_pls <- exp(y_pred_pls)/( 1 + exp(y_pred_pls))
  y_pred_pls <- y_pred_pls * (y_train_max- y_train_min)
  y_pred_pls <- y_pred_pls + y_train_min

  L_trend_0_logit[L] <- mse(y_pred_pls, y_test)
}

plot(L_trend_0_logit[1:L_max], main =
       paste("kappa=1e-6 and 5*1e-6, tau = 0, L_max = 20 ",
             "\n MSE for PLS with logit transform\n", "seed = ", seed, " curves NOT scaled"
       ), ylab = "MSE", xlab = "number of iteration", pch=".", ylim=c(0.04, 0.16))
lines(L_trend_0_logit)
print(min(L_trend_0_logit, omit.NA = TRUE))

# scale= True, seed = 1000 -------------------------------------------------
seed = 1000
set.seed(seed) # for training/test split

#set the number of max iterations
L_max = 20


kidney <- read_fd_kidney(
  dir = "~/GitHub/fingernail/data/renogram_data.csv",
  response_function = response_mean_diagnosis,
  training_ratio = 0.7,
  n_basis = 20,
  scale = FALSE
)

# W is already centered
W_train_centered <- kidney$W_train
y_train <- kidney$y_train

# minmax scale and logit transform and centering
y_train_min <- min(y_train)
y_train_max <- max(y_train)

y_train_logit <- response_transform_min_max_logit(y_train,
                                                  y_train_min,
                                                  y_train_max)
y_train_logit_mean <- mean(y_train_logit) #centering again, for PLS algorithm
y_train_logit_centered <- y_train_logit - y_train_logit_mean


# test data set
W_test_centered <- kidney$W_test
y_test <- kidney$y_test


#MSE trend

L_trend_0_logit <- rep(NA, L_max)
for(L in 1:L_max){
  pls_model_transform <- nipals_pen_hybrid(
    W_train_centered, y_train_logit_centered, L, c(1e-6,5*1e-6), 0)

  #inverse transform
  y_pred_pls <- predict_test(pls_model_transform, W_test_centered)
  y_pred_pls <- y_pred_pls + y_train_logit_mean
  y_pred_pls <- exp(y_pred_pls)/( 1 + exp(y_pred_pls))
  y_pred_pls <- y_pred_pls * (y_train_max- y_train_min)
  y_pred_pls <- y_pred_pls + y_train_min

  L_trend_0_logit[L] <- mse(y_pred_pls, y_test)
}

plot(L_trend_0_logit[1:L_max], main =
       paste("kappa=1e-6 and 5*1e-6, tau = 0, L_max = 20 ",
             "\n MSE for PLS with logit transform\n", "seed = ", seed, " curves NOT scaled"
       ), ylab = "MSE", xlab = "number of iteration", pch=".", ylim=c(0.04, 0.16))
lines(L_trend_0_logit)
print(min(L_trend_0_logit, omit.NA = TRUE))
