#
##### save RMSE #####
n_simul = 50
n_iter = 10
n_method = 3
mse <- matrix(NA, nrow = n_simul, ncol = n_method)
colnames(mse) <- c("lin", "pc", "hybpls")
##### save setting of interest #####
# n_sample
param_list <- c(30)
rmse_of_interest_mean <-  matrix(NA, nrow = length(param_list), ncol = n_method)
rownames(rmse_of_interest_mean) <- param_list
rmse_of_interest_sd <- rmse_of_interest_mean
##################
test.ratio <- 0.3
n_eval <- 401
eval_point <- seq(0,1, length = n_eval)
n_predictor_scalar <- 6
n_pc_scalar <- floor(n_predictor_scalar/2)
n_pc_functional <- 2

n_basis_predictor <- 30
n_basis_coef <- 30
basis_reg_coef <- create.fourier.basis(rangeval = c(0,1), nbasis = n_basis_coef)
basis_predictor <- create.fourier.basis(rangeval = c(0,1), nbasis = n_basis_predictor)

##### regression coefficients, fixed over simulations #####
## first functional regression coefficient
set.seed(4)
reg_coef_1 <- fda::Data2fd(argvals = eval_point, y = as.vector(arima.sim(model = list(ar = -0.5), n = n_eval)), basisobj = basis_reg_coef)
plot(reg_coef_1, main = "functional regression coef 1")

## second functional regression coefficient
set.seed(5)
reg_coef_2 <- fda::Data2fd(argvals = eval_point, y = as.vector(arima.sim(model = list(ar = 0.5), n = n_eval)), basisobj = basis_reg_coef)
plot(reg_coef_2, main = "functional regression coef 2")



# scalar regression coefficient
alpha <- rnorm(n_predictor_scalar, 0, 1)/(6*n_predictor_scalar^(1/8))


##### predictor generation #####
for (param_number in 1 : length(param_list)){
  n_sample <- param_list[param_number]
  for (simul_number in 1:n_simul){
    set.seed(simul_number)

    n_sample <- 50
    ar_slope <- 0.9
    gaussian_bandwidth <- 3

    freq_1 <- 10
    freq_2 <- 15
    diff <- freq_2 - freq_1
    scalar_freq <- seq(freq_1 - diff/2, freq_2 + diff/2,length = n_predictor_scalar+2)
    predictor.value <- generate_ar_wavelet_conv(
      n_sample, n_eval, ar_slope, gaussian_bandwidth,
      c(freq_1,freq_2), scalar_freq[2:(n_predictor_scalar+1)]
    )

    predictor_func_1 <- fda::Data2fd(argvals = eval_point,
                                     y = t(predictor.value$x_functional$first$value),
                                     basisobj = basis_predictor)
    predictor_func_2 <- fda::Data2fd(argvals = eval_point,
                                     y = t(predictor.value$x_functional$second$value),
                                     basisobj = basis_predictor)
    # noise


    signal_functional <- inprod(reg_coef_1, predictor_func_1) + inprod(reg_coef_2, predictor_func_2)
    signal_scalar <- t(predictor.value$x_scalar %*% alpha)
    signal <- signal_functional + signal_scalar
    sd(signal_functional)
    sd(signal_scalar)
    noise <- rnorm(n_sample, 0, 2)
    #range(noise)
    #range(signal)
    response <- matrix(signal+ noise, ncol = 1)

    range(response)

    #####

    #data packaging pipeline

    # step 1. for hybridPLS
    predictor.value$y <- response
    simul_data_functional <- read_fd_simul(
      test_ratio = 0.3,
      n_basis = n_basis_predictor,
      normalize = list("curve" = TRUE, "scalar"= TRUE, "between" = FALSE),
      predictor.value
    )

    # W is already centered
    W_train_centered <- simul_data_functional$W_train
    y_train <- simul_data_functional$y_train

    # test data set
    W_test_centered <- simul_data_functional$W_test
    y_test <- simul_data_functional$y_test

    # step 2. for linear regression
    dataset_for_regression <- dataset_for_regression_refund(y_train, W_train_centered)
    dataset_for_prediction <- dataset_for_prediction_refund(W_test_centered)

    ##############################################################################################
    # 1. scalar-on-function regression + scalar covariate
    # table to summarize the result



    library(refund)

    # all predictors
    reg.string <- paste0(
      "soft_all.fit <- pfr(y ~lf(F1,argvals = eval_point,presmooth='bspline',presmooth.opts=list(nbasis=n_basis_predictor))+lf(F2,argvals = eval_point,presmooth='bspline',presmooth.opts=list(nbasis=n_basis_predictor))+",
      paste(paste0("X", 1:n_predictor_scalar), collapse="+"),
      ", data = dataset_for_regression)"
    )
    eval(parse(text = reg.string))


    mse[simul_number,1] <- mse(predict(soft_all.fit, dataset_for_prediction), y_test)
    mse
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
    )*100000



    #MSE trend
    #n_iter = ceil((n_predictor_scalar + n_predictor_functional)/2)
    #lambda = c(1,1)
    iter_trend <- rep(NA, n_iter)

    #learn the model
    hybpls.fit<- nipals_pen_hybrid(W_train_centered, y_train, n_iter, lambda, 0)

    # predict
    for (iter in 1 : n_iter){
      y_pred_pls <- predict(hybpls.fit, W_test_centered, iter)
      iter_trend[iter] <- mse(y_pred_pls, y_test)
    }
    par(mfrow=c(2,1))
    plot(iter_trend, main ="hybridPLS", ylab = "MSE", xlab = "number of iteration", pch=".")
    lines(iter_trend)


    mse[simul_number,3] <- iter_trend[n_iter]
  }


  rmse_scaled <- sqrt(mse)/sd(signal)
  rmse_of_interest_mean[param_number, ] <- apply(rmse_scaled, 2, mean)
  rmse_of_interest_sd[param_number, ] <- apply(rmse_scaled, 2, sd)
}
rmse_of_interest_mean
rmse_of_interest_sd
sd(signal_functional)
sd(signal_scalar)


cor_func <- cor(predictor.value$x_functional$first$value, predictor.value$x_functional$second$value)
hist(diag(cor_func))
cor_between_1 <- cor(predictor.value$x_functional$first$value, predictor.value$x_scalar[,2])
hist(cor_between_1)
cor_between_2 <- cor(predictor.value$x_functional$first$value, predictor.value$x_scalar)
hist(cor_between_2)

cor(predictor.value$x_scalar)

