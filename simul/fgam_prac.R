test.ratio <- 0.3
n_eval <- 401
eval_point <- seq(0,1, length = n_eval)
n_predictor_scalar <- 2
n_pc_scalar <- floor(n_predictor_scalar/2)
n_pc_functional <- 2

n_basis_predictor <- 30
n_basis_coef <- 15
basis_reg_coef <- create.bspline.basis(rangeval = c(0,1), nbasis = n_basis_coef)
basis_predictor <- create.bspline.basis(rangeval = c(0,1), nbasis = n_basis_predictor)

##### regression coefficients, fixed over simulations #####
## first functional regression coefficient
set.seed(10)
reg_coef_1 <- fda::Data2fd(argvals = eval_point, y = as.vector(arima.sim(model = list(ar = -0.5), n = n_eval)), basisobj = basis_reg_coef)
plot(reg_coef_1)
## second functional regression coefficient
set.seed(5)
reg_coef_2 <- fda::Data2fd(argvals = eval_point, y = as.vector(arima.sim(model = list(ar = 0.5), n = n_eval)), basisobj = basis_reg_coef)
plot(reg_coef_2)

alpha <- rnorm(n_predictor_scalar, 0, 1)/(6*n_predictor_scalar^(1/8))


n_sample <- 500
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
plot(predictor.value$x_functional$first$value, pch=".")
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
noise <- rnorm(n_sample, 0, 1)
#range(noise)
#range(signal)
response <- matrix(signal + noise, ncol = 1)

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

library(refund)
soft_all.fit <- pfr(
  y_train~
    fpc(F1,argvals = eval_point, bs = "bs", k = n_basis_predictor, ncomp=20)+
  fpc(F2,argvals = eval_point, bs = "bs", k = n_basis_predictor, ncomp=20)+
  X1+X2,
  data = dataset_for_regression)

predict(soft_all.fit, newdata = list(dataset_for_prediction))
soft_all.fit$coefficients

par(mfrow=c(2,1))
plot(reg_coef_2)
plot(value~X.argvals, data = coef(soft_all.fit), pch=".")
par(mfrow=c(1,1))
plot(reg_coef_1)
str(soft_all.fit)
summary(soft_all.fit)
predictor_func_1_pca <- fda::Data2fd(argvals = eval_point,
                                 y = t(W_train_centered@predictor_functional_list[[1]]@original_X),
                                 basisobj = basis_predictor)
predictor_func_2_pca <- fda::Data2fd(argvals = eval_point,
                                 y = t(W_train_centered@predictor_functional_list[[2]]@original_X),
                                 basisobj = basis_predictor)

fpca_fit_1 <- pca.fd(predictor_func_1_pca, nharm = 5)
fpca_fit_2 <- pca.fd(predictor_func_1_pca, nharm = 5)
pca.fit.scalar <- prcomp(W_train_centered@Z)
fpca_score_train_1 <- fpca_fit_1$scores
fpca_score_train_2 <- fpca_fit_1$scores
pca_score_train <- (pca.fit.scalar$x)[,1:n_pc_scalar]
training_data_pc <- cbind(fpca_score_train_1, fpca_score_train_2, pca_score_train)

J_dotdot <- getbasispenalty(basis_predictor)
penalty_matrix <- bdiag(J_dotdot, J_dotdot, eye(n_predictor_scalar))

