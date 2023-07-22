#
n_sample <- 170
n_eval <- 49
n_basis <- 10
M <- matrixNormal::J(n_sample, n_basis) # mean = all 1, n << N
U <- matrixNormal::I(n_sample) # n x n matrix. assume independence between observations
V_1 <- matrix(rnorm(n_basis*n_basis), nrow = n_basis) # N x N matrix. assume dependence across time
V_1 <- V_1 %*% t(V_1)
V_2 <- matrix(rnorm(n_basis*n_basis), nrow = n_basis) # N x N matrix. assume dependence across time
V_2 <- V_2 %*% t(V_2)
between_predictor_correlation <- matrix(rnorm(n_sample*n_basis), nrow = n_sample)

basis_coef_1 <- matrixNormal::rmatnorm(M = M, U = U, V=V_1) + between_predictor_correlation
basis_coef_2 <- matrixNormal::rmatnorm(M = M, U = U, V=V_2) + between_predictor_correlation

my_basis <- create.bspline.basis(rangeval = c(0,1), nbasis = n_basis)
predictor_func_1 <- fda::fd(coef = t(basis_coef_1), basisobj = my_basis)
predictor_func_2 <- fda::fd(coef = t(basis_coef_2), basisobj = my_basis)
par(mfrow=c(1,1))
plot(predictor_func_1[1])

t <- seq(0,1, length = n_eval)
eval.fd(predictor_func_1, t)
#scalar predictors
p = 15
V_Z <- rnorm(p,0,1)
V_Z <- V_Z %*% t(V_Z)
Z <- mvtnorm::rmvnorm(n = n_sample, mean = rep (0, p), sigma = V_Z)


#y
my_basis <- create.bspline.basis(rangeval = c(0,1), nbasis = n_basis)
argvals <- t
PhiB <- predict(my_basis, argvals)
PhiB <- PhiB[nrow(PhiB):1, ncol(PhiB):1]
J <- get_gram_2d(argvals, PhiB)

b_1 <- rnorm(n_basis, 0, 1)
b_2 <- rnorm(n_basis, 0, 1)
alpha <- rnorm(p, 0, 1)
response <- basis_coef_1 %*% J %*% b_1 +
  basis_coef_2 %*% J %*% b_2 +
  Z %*% alpha +
  rnorm(n_sample, 0, 1)
response <- matrix(response, ncol = 1)

value_object <- list(
  "x_functional" = list(
    "first" = list("value" = t(eval.fd(predictor_func_1, argvals)), "timestamp" = t(t %*% t(rep(1, n_sample)))),
    "second" = list("value" = t(eval.fd(predictor_func_2, argvals)), "timestamp" = t(t %*% t(rep(1, n_sample))))
  ),
  "x_scalar" = Z,
  "y" = response
)


kidney_value <- preprocess_reno(value_object) #preprocessing, only for Emory kidney data
kidney_value_split <- train_test_split(kidney_value, 0.3)
kidney_value_train <- kidney_value_split$train
kidney_value_test <- kidney_value_split$test

simul_data_functional <- read_fd_simul(
  test_ratio = 0.3,
  n_basis = n_basis,
  normalize = list("curve" = TRUE, "scalar"= TRUE, "between" = TRUE),
  value_object
)



# W is already centered
W_train_centered <- simul_data_functional$W_train
y_train <- simul_data_functional$y_train



# test data set
W_test_centered <- simul_data_functional$W_test
y_test <- simul_data_functional$y_test


##############################################################################################
# 1. scalar-on-function regression + scalar covariate
# table to summarize the result
regression_result <- matrix(NA, nrow = 2, ncol = 3)
rownames(regression_result) <- c("train", "test")
colnames(regression_result) <- c("first", "second", "all")

# load data - matrix format to use refund package
# since the pfr function does the smoothing for us
timepoints_pre <- W_train_centered@predictor_functional_list[[1]]@original_t[1,]
timepoints_post <- W_train_centered@predictor_functional_list[[2]]@original_t[1,]
dataset_for_regression <- dataset_for_regression_refund(y_train, W_train_centered)
dataset_for_prediction <- dataset_for_prediction_refund(W_test_centered)


soft.fit <- refund::pfr(
  response ~
    refund::lf( # first functional predictor
      precurve,
      argvals = timepoints_pre,
      presmooth="bspline",
      presmooth.opts=list(nbasis=n_basis)
    ),
  data = dataset_for_regression
)


# first functional predictors
soft.fit <- refund::pfr(
  response ~
    refund::lf( # first functional predictor
      precurve,
      argvals = timepoints_pre,
      presmooth="bspline",
      presmooth.opts=list(nbasis=N)
    )
  +
    X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + X11 + X12 + X13 + X14 + X15,
  data = dataset_for_regression
)


regression_result[1,1] <- mse(predict(soft.fit), y_train)
regression_result[2,1] <- mse(predict(soft.fit, dataset_for_prediction), y_test)


# second functional predictors
soft.fit <- refund::pfr(
  y_train ~
    refund::lf( # second functional predictor
      postcurve,
      argvals = timepoints_post,
      presmooth="bspline",
      presmooth.opts=list(nbasis=n_basis)
    ) +
    X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + X11 + X12 + X13 + X14 + X15,
  data = dataset_for_regression
)

regression_result[1,2] <- mse(predict(soft.fit), y_train)
regression_result[2,2] <- mse(predict(soft.fit, dataset_for_prediction), y_test)


# all functional predictors
soft.fit <- refund::pfr(
  y_train ~
    refund::lf( # first functional predictor
      precurve,
      argvals = timepoints_pre,
      presmooth="bspline",
      presmooth.opts=list(nbasis=n_basis)
    )
  +
    lf( # second functional predictor
      postcurve,
      argvals = timepoints_post,
      presmooth="bspline",
      presmooth.opts=list(nbasis=n_basis)
    ) +
    X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + X11 + X12 + X13 + X14 + X15,
  data = dataset_for_regression
)


regression_result[1,3] <- mse(predict(soft.fit), y_train)
regression_result[2,3] <- mse(predict(soft.fit, dataset_for_prediction), y_test)
mse[1] <- mse(predict(soft.fit, dataset_for_prediction), y_test)
mse[1]
regression_result




#2. PCA regression: MFPCA score of functional predictors + PCA score of scalar variable
############### CANNOT USE TWO FUNCTIONAL PREDICTORS#########







# 5.Our method: HybridPLS


#MSE trend
L_max = 40
lambda = c(1e-6,5*1e-6) #hyperparameters
L_trend_0_logit <- rep(NA, L_max)
for(L in 1:L_max){
  #learn the model
  pls_model_transform <- nipals_pen_hybrid(
    W_train_centered, y_train, L, lambda, 0)

  # predict
  y_pred_pls <- predict_test(pls_model_transform, W_test_centered)

  #inverse transform

  # calculate the mse and save
  L_trend_0_logit[L] <- mse(y_pred_pls, y_test)
}

plot(L_trend_0_logit[1:L_max], main ="hybridPLS", ylab = "MSE", xlab = "number of iteration", pch=".")
lines(L_trend_0_logit)
mse[5] <- min(L_trend_0_logit)
mse
