#
##### save RMSE #####
n_simul = 100
n_method = 3
rmse <- matrix(NA, nrow = n_simul, ncol = n_method)
colnames(rmse) <- c("lin", "pc", "hybpls")
##### save setting of interest #####
# n_sample
param_list <- c(50)
rmse_of_interest_mean <-  matrix(NA, nrow = length(param_list), ncol = n_method)
rownames(rmse_of_interest_mean) <- param_list
rmse_of_interest_sd <- rmse_of_interest_mean
##################
test.ratio <- 0.3
n_eval <- 201
eval_point <- seq(0,1, length = n_eval)
n_predictor_scalar <- 4
n_predictor_functional <- 2
n_pc_scalar <- floor(n_predictor_scalar/2)
n_pc_functional <- 2

n_basis_predictor <- 7
n_basis_coef <- 7
basis_reg_coef <- create.bspline.basis(rangeval = c(0,1), nbasis = n_basis_coef)
my_basis <- create.bspline.basis(rangeval = c(0,1), nbasis = n_basis_predictor)

##### regression coefficients, fixed over simulations #####
## first functional regression coefficient
set.seed(4)
reg_coef_1 <- fda::Data2fd(argvals = eval_point, y = as.vector(arima.sim(model = list(ar = -0.5), n = n_eval)), basisobj = basis_reg_coef)
plot(reg_coef_1, main = "functional regression coef 1")

## second functional regression coefficient
set.seed(5)
reg_coef_2 <- fda::Data2fd(argvals = eval_point, y = as.vector(arima.sim(model = list(ar = 0.5), n = n_eval)), basisobj = basis_reg_coef)
plot(reg_coef_2, main = "functional regression coef 2")
plot(soft_all.fit)


# scalar regression coefficient
alpha <- rnorm(n_predictor_scalar, 0, 1)


##### predictor generation #####
for (param_number in 1 : length(param_list)){
  n_sample <- param_list[param_number]
  for (simul_number in 1:n_simul){
    set.seed(simul_number)
    ar_slope <- 0.9


    predictor_func_1 <- fda::Data2fd(argvals = eval_point, y = t(functional_data_eval_1), basisobj = my_basis)
    predictor_func_2 <- fda::Data2fd(argvals = eval_point, y = t(functional_data_eval_2), basisobj = my_basis)

    # oversmoothed, so replace the original evals with smoothed eval
    functional_data_eval_1 = t(eval.fd(eval_point, predictor_func_1))
    functional_data_eval_2 = t(eval.fd(eval_point, predictor_func_2))

    predictor_func_1 <- fda::Data2fd(argvals = eval_point, y = t(functional_data_eval_1), basisobj = my_basis)
    predictor_func_2 <- fda::Data2fd(argvals = eval_point, y = t(functional_data_eval_2), basisobj = my_basis)

    plot(predictor_func_1)
    plot(predictor_func_2)


    #between_predictor_correlation <- matrix(rnorm(n_sample*n_basis), nrow = n_sample)


    #scalar predictors
    frequency_scalar_1 = 50
    frequency_scalar_2 = 60
    frequency_scalar_3 = 70
    frequency_scalar_4 = 80

    X1 = apply(
      filter_signal(signal.original, eval_point, frequency_scalar_1, gaussian_bandwidth),
      1,
      mean)

    X2 = apply(
      filter_signal(signal.original, eval_point, frequency_scalar_2, gaussian_bandwidth),
      1,
      mean)

    X3 = apply(
      filter_signal(signal.original, eval_point, frequency_scalar_3, gaussian_bandwidth),
      1,
      mean)

    X4 = apply(
      filter_signal(signal.original, eval_point, frequency_scalar_4, gaussian_bandwidth),
      1,
      mean)

    Z <- cbind(X1, X2, X3, X4)/30

    # noise
    noise <- rnorm(n_sample, 0, 0.1)

    ##### response generation #####


    range(abs(inprod(reg_coef_1, predictor_func_1)))
    range(abs(inprod(reg_coef_2, predictor_func_2)))
    range(t(Z %*% alpha))

    signal <- inprod(reg_coef_1, predictor_func_1) + inprod(reg_coef_2, predictor_func_2) + t(Z %*% alpha)

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


    dataset_for_regression <- dataset_for_regression_refund(y_train, W_train_centered)
    dataset_for_prediction <- dataset_for_prediction_refund(W_test_centered)

    library(refund)

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
        ) + X1 + X2 + X3 + X4,
      data = dataset_for_regression
    )


    rmse[simul_number,1] <- mse(predict(soft_all.fit, dataset_for_prediction), y_test)
    rmse
    ##############################################################################################
    # # 2. FPCA + PCA
    #
    # dataset_for_PCA_train <- dataset_for_PCA(dataset_for_regression, eval_point)
    # dataset_for_PCA_test <- dataset_for_PCA(dataset_for_prediction, eval_point)
    #
    # # The first n samples are for training and the rest testing
    # pca.fit.F1 <- FPCA(dataset_for_PCA_train$F1, dataset_for_PCA_train$eval_point)
    # pca.fit.F2 <- FPCA(dataset_for_PCA_train$F2, dataset_for_PCA_train$eval_point)
    # pca.fit.scalar <- prcomp(W_train_centered@Z) #learn PCA with *train* data
    #
    # pca.data.train <- data.frame(
    #   y_train,
    #   fpca.fit.F1$xiEst[,1:n_pc_functional],
    #   fpca.fit.F2$xiEst[,1:n_pc_functional],
    #   (pca.fit.scalar$x)[,1:n_pc_scalar]
    #   )
    # colnames(pca.data.train)[1] <- "y"
    # pca.data.test <- data.frame(
    #   y_test,
    #   predict(fpca.fit.F1, dataset_for_PCA_test$F1, dataset_for_PCA_test$eval_point, K = n_pc_functional)$scores,
    #   predict(fpca.fit.F2, dataset_for_PCA_test$F2, dataset_for_PCA_test$eval_point, K = n_pc_functional)$scores,
    #   (predict(pca.fit.scalar, W_test_centered@Z))[,1:n_pc_scalar]
    # )
    # colnames(pca.data.test)[1] <- "y"
    #
    # pca.reg.fit <-  lm(y ~. , data = pca.data.train)
    # rmse[simul_number,2] <- mse(predict(pca.reg.fit, pca.data.test), y_test)

    # 5.Our method: HybridPLS
    lambda <- c( #hyperparameters
      max(abs(W_train_centered@predictor_functional_list[[1]]@J))/max(abs((W_train_centered@predictor_functional_list[[1]]@J_dotdot))),
      max(abs(W_train_centered@predictor_functional_list[[2]]@J))/max(abs((W_train_centered@predictor_functional_list[[2]]@J_dotdot)))
    )/10



    #MSE trend
    #n_iter = ceil((n_predictor_scalar + n_predictor_functional)/2)
    n_iter = 10
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
    rmse[simul_number,3] <- min(iter_trend)
  }
  rmse <- sqrt(rmse)
  rmse_of_interest_mean[param_number, ] <- apply(rmse, 2, mean)
  rmse_of_interest_sd[param_number, ] <- apply(rmse, 2, sd)
}
rmse_of_interest_mean
rmse_of_interest_sd
