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


# FPCA predictor 1
W_train_fd <- extract_fd(W_train_centered)
W_test_fd <- extract_fd(W_test_centered)

base.pcalist = pca.fd(W_train_fd[[1]], n_basis)
print(base.pcalist$values)
## train data
discrete_data_base_train <- t( eval.fd(W_train_fd[[1]], timepoints_pre) ) #evalauation of the normalized train data at the observed timepoints
PhiB_FPCA_base <- eval.fd(base.pcalist$harmonics, timepoints_pre) # evalauation of the FPCA components(eigenfunctions) at the observed timepoints
dim(PhiB_FPCA_base)
FPCA_score_base_train <- t( MASS::ginv(PhiB_FPCA_base) %*% t(discrete_data_base_train) ) #FPCA score, base, train data
FPCA_score_base_train <- FPCA_score_base_train[,1:2]
## test data
discrete_data_base_test <- t( eval.fd(W_test_fd[[1]], timepoints_pre) )
FPCA_score_base_test <- t( MASS::ginv(PhiB_FPCA_base) %*% t(discrete_data_base_test) )
FPCA_score_base_test <- FPCA_score_base_test[,1:2]

# FPCA predictor 2
t_point_post <- (W_train_centered@predictor_functional_list[[2]]@original_t)[1,]
post.pcalist = pca.fd(W_train_fd[[2]], n_basis)
print(post.pcalist$values)
## train data
discrete_data_post_train <- t( eval.fd(W_train_fd[[2]], timepoints_post) ) #evalauation of the normalized train data at the observed timepoints
PhiB_FPCA_post <- eval.fd(post.pcalist$harmonics, timepoints_post) # evalauation of the FPCA components(eigenfunctions) at the observed timepoints
dim(PhiB_FPCA_post)
FPCA_score_post_train <- t( MASS::ginv(PhiB_FPCA_post) %*% t(discrete_data_post_train) ) #FPCA score, base, train data
FPCA_score_post_train <- FPCA_score_post_train[,1:2]
## test data
discrete_data_post_test <- t( eval.fd(W_test_fd[[1]], timepoints_post) )
FPCA_score_post_test <- t( MASS::ginv(PhiB_FPCA_post) %*% t(discrete_data_post_test) )
FPCA_score_post_test <- FPCA_score_post_test[,1:2]

PCA_scalar <- prcomp(W_train_centered@Z) #learn PCA with train data
summary(PCA_scalar)
PCA_score_train <- (PCA_scalar$x)[,1:7]
PCA_score_test <- (predict(PCA_scalar, W_test_centered@Z))[,1:7]

#combine
PCA_data_train <- data.frame(cbind(y_train, FPCA_score_base_train, FPCA_score_post_train, PCA_score_train))
PCA_data_test <- data.frame(cbind(FPCA_score_base_test, FPCA_score_post_test, PCA_score_test))
colnames(PCA_data_test) <- colnames(PCA_data_train)[-1]
#regress
pcareg <- lm(y_train~., data = PCA_data_train)
#inverse transform

y_pred_pls <- reponse_inverse_transform_min_max_logit(
  predict(pcareg), y_train_logit_mean, y_train_max, y_train_min) # training fit
mse(y_pred_pls, y_train)


y_pred_pls <- reponse_inverse_transform_min_max_logit(
  predict(pcareg, PCA_data_test), y_train_logit_mean, y_train_max, y_train_min) # test fit
mse[2] <- mse(y_pred_pls, y_test)
mse[2]






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
