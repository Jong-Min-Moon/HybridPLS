#
n_simul = 10
n_method = 2
rmse <- matrix(NA, nrow = n_simul, ncol = n_method)
colnames(rmse) <- c("lin", "hybpls")
n_sample <- 170
test.ratio <- 0.3
n_eval <- 49
n_basis <- 10
eval_point <- seq(0,1, length = n_eval)
n_scalar_predictor <- 5

n_basis_predictor <- 15
n_basis_coef <- 7
basis_reg_coef <- create.bspline.basis(rangeval = c(0,1), nbasis = n_basis)
my_basis <- create.bspline.basis(rangeval = c(0,1), nbasis = n_basis)

##### regression coefficients, fixed over simulations #####
## first functional regression coefficient
set.seed(4)
reg_coef_1 <- fda::Data2fd(argvals = eval_point, y = as.vector(arima.sim(model = list(ar = -0.5), n = n_eval)), basisobj = basis_reg_coef)
plot(reg_coef_1)

## second functional regression coefficient
set.seed(5)
reg_coef_2 <- fda::Data2fd(argvals = eval_point, y = as.vector(arima.sim(model = list(ar = 0.5), n = n_eval)), basisobj = basis_reg_coef)
plot(reg_coef_2)

# scalar regression coefficient
alpha <- rnorm(n_scalar_predictor, 0, 1)
alpha


##### predictor generation #####
for (simul_number in 1:n_simul){
set.seed(simul_number)
ar_slope_1 <- 0.9
ar_slope_2 <- 0.5
functional_data_eval_1 <- functional_data_eval_2 <-  matrix(NA, nrow = n_sample, ncol = n_eval)
for (i in 1:n_sample){
  functional_data_eval_1[i,] <- arima.sim(model = list(ar = ar_slope_1), n = n_eval)
  functional_data_eval_2[i,] <- arima.sim(model = list(ar = ar_slope_2), n = n_eval)
}
functional_data_eval_1 <- functional_data_eval_1

functional_data_eval_2 <- functional_data_eval_2
#between_predictor_correlation <- matrix(rnorm(n_sample*n_basis), nrow = n_sample)


#scalar predictors
#V_Z <- rnorm(n_scalar_predictor,0,1/5)
#V_Z <- V_Z %*% t(V_Z)
V_Z <- eye(5)/30
Z <- mvtnorm::rmvnorm(n = n_sample, mean = rep (0, n_scalar_predictor), sigma = V_Z)
plot(Z)
# noise
noise <- rnorm(n_sample, 0, 1/10)

##### response generation #####
predictor_func_1 <- fda::Data2fd(argvals = eval_point, y = t(functional_data_eval_1), basisobj = my_basis)
predictor_func_2 <- fda::Data2fd(argvals = eval_point, y = t(functional_data_eval_2), basisobj = my_basis)

signal <- inprod(reg_coef_1, predictor_func_1) + inprod(reg_coef_2, predictor_func_2) + t(Z %*% alpha) #+ 100
range(inprod(reg_coef_1, predictor_func_1))
range(inprod(reg_coef_2, predictor_func_2))
range(t(Z %*% alpha))
range(noise)
range(signal)
response <- matrix(signal+ noise, ncol = 1)

range(response)

#####

value_object <- list(
  "x_functional" = list(
    "first" = list("value" = functional_data_eval_1, "timestamp" = t(eval_point %*% t(rep(1, n_sample)))),
    "second" = list("value" = functional_data_eval_2, "timestamp" = t(eval_point %*% t(rep(1, n_sample))))
  ),
  "x_scalar" = Z,
  "y" = response
)

simul_data_functional <- read_fd_simul(
  test_ratio = 0.3,
  n_basis = n_basis,
  normalize = list("curve" = TRUE, "scalar"= TRUE, "between" = FALSE),
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
colnames(regression_result) <- c("f1", "f1+f2", "all")

dataset_for_regression <- dataset_for_regression_refund(y_train, W_train_centered)
dataset_for_prediction <- dataset_for_prediction_refund(W_test_centered)

library(refund)
# first functional predictors
soft_f1.fit <- pfr(
  y ~
    lf( # first functional predictor
      F1,
      argvals = eval_point,
      presmooth="bspline",
      presmooth.opts=list(nbasis=n_basis_predictor)
    ),
  data = dataset_for_regression
)


regression_result[1,1] <- mse(predict(soft_f1.fit), dataset_for_regression$y)
regression_result[2,1] <- mse(predict(soft_f1.fit, dataset_for_prediction), y_test)
regression_result

# first and second functional predictors
soft_f_all.fit <- pfr(
  y ~
    lf( # first functional predictor
      F1,
      argvals = eval_point,
      presmooth="bspline",
      presmooth.opts=list(nbasis=n_basis_predictor)
    ) +
    lf( # first functional predictor
      F2,
      argvals = eval_point,
      presmooth="bspline",
      presmooth.opts=list(nbasis=n_basis_predictor)
    ),
  data = dataset_for_regression
)

regression_result[1,2] <- mse(predict(soft_f_all.fit), dataset_for_regression$y)
regression_result[2,2] <- mse(predict(soft_f_all.fit, dataset_for_prediction), y_test)
regression_result

# all predictors
soft_all.fit <- pfr(
  y ~
    lf( # first functional predictor
      F1,
      argvals = eval_point,
      presmooth="bspline",
      presmooth.opts=list(nbasis=n_basis_predictor)
    ) +
    lf( # first functional predictor
      F2,
      argvals = eval_point,
      presmooth="bspline",
      presmooth.opts=list(nbasis=n_basis_predictor)
    ) + X1 + X2 + X3 + X4 + X5,
  data = dataset_for_regression
)


regression_result[1,3] <- mse(predict(soft_all.fit), dataset_for_regression$y)
regression_result[2,3] <- mse(predict(soft_all.fit, dataset_for_prediction), y_test)
regression_result

rmse[simul_number,1] <- regression_result[2,3]
rmse





# 5.Our method: HybridPLS
lambda <- c( #hyperparameters
  max(abs(W_train_centered@predictor_functional_list[[1]]@J))/max(abs((W_train_centered@predictor_functional_list[[1]]@J_dotdot))),
  max(abs(W_train_centered@predictor_functional_list[[2]]@J))/max(abs((W_train_centered@predictor_functional_list[[2]]@J_dotdot)))
)



#MSE trend
n_iter = 40
#lambda = c(1,1)
iter_trend <- rep(NA, n_iter)

#learn the model
hybpls.fit<- nipals_pen_hybrid(W_train_centered, y_train, n_iter, lambda, 0)

# predict
for (iter in 1 : n_iter){
  y_pred_pls <- predict(hybpls.fit, W_test_centered, iter)
  iter_trend[iter] <- mse(y_pred_pls, y_test)
}
plot(iter_trend, main ="hybridPLS", ylab = "MSE", xlab = "number of iteration", pch=".")
lines(iter_trend)
rmse[simul_number,2] <- min(iter_trend)
}

write.csv(rmse, "result.csv")
