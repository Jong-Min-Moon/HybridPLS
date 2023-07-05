#
n_sample <- 170
test.ratio <- 0.3
n_eval <- 49
n_basis <- 10
eval_point <- seq(0,1, length = n_eval)
n_scalar_predictor <- 15

n_basis_predictor <- 15
n_basis_coef <- 7
basis_reg_coef <- create.bspline.basis(rangeval = c(0,1), nbasis = n_basis)
my_basis <- create.bspline.basis(rangeval = c(0,1), nbasis = n_basis)

##### regression coefficients, fixed over simulations #####
## first functional regression coefficient
set.seed(4)
reg_coef_1 <- fda::Data2fd(argvals = eval_points, y = as.vector(arima.sim(model = list(ar = -0.5), n = n_eval)), basisobj = basis_reg_coef)
plot(reg_coef_1)

## second functional regression coefficient
set.seed(5)
reg_coef_2 <- fda::Data2fd(argvals = eval_points, y = as.vector(arima.sim(model = list(ar = 0.5), n = n_eval)), basisobj = basis_reg_coef)
plot(reg_coef_2)

# scalar regression coefficient
alpha <- rnorm(n_scalar_predictor, 0, 1)
alpha


##### predictor generation #####
set.seed(1)

ar_slope_1 <- 0.9
ar_slope_2 <- 0.5
functional_data_eval_1 <- functional_data_eval_2 <-  matrix(NA, nrow = n_sample, ncol = n_eval)
for (i in 1:n_sample){
  functional_data_eval_1[i,] <- arima.sim(model = list(ar = ar_slope_1), n = n_eval)
  functional_data_eval_2[i,] <- arima.sim(model = list(ar = ar_slope_2), n = n_eval)
}

# ar process into functional objects
predictor_func_1 <- fda::Data2fd(argvals = eval_points, y = t(functional_data_eval_1), basisobj = my_basis)
predictor_func_2 <- fda::Data2fd(argvals = eval_points, y = t(functional_data_eval_2), basisobj = my_basis)

#between_predictor_correlation <- matrix(rnorm(n_sample*n_basis), nrow = n_sample)


#scalar predictors
V_Z <- rnorm(n_scalar_predictor,0,1/10)
V_Z <- V_Z %*% t(V_Z)
Z <- mvtnorm::rmvnorm(n = n_sample, mean = rep (0, n_scalar_predictor), sigma = V_Z)

# noise
noise <- rnorm(n_sample, 0, 1/10)

##### response generation #####

range(inprod(reg_coef_1, predictor_func_1))
range(inprod(reg_coef_2, predictor_func_2))
range(t(Z %*% alpha))
range(noise)
response <- matrix(
  inprod(reg_coef_1, predictor_func_1) + inprod(reg_coef_2, predictor_func_2) + t(Z %*% alpha) + noise,
  ncol = 1
)

range(response)

#####

##############################################################################################
# 1. scalar-on-function regression + scalar covariate
# table to summarize the result
regression_result <- matrix(NA, nrow = 2, ncol = 3)
rownames(regression_result) <- c("train", "test")
colnames(regression_result) <- c("f1", "f1+f2", "all")

test.idx <- sample(n_sample, floor(test.ratio* n_sample))
dataset_for_regression <- data.frame(y = response[-test.idx])
dataset_for_regression$X1 = functional_data_eval_1[-test.idx,]
dataset_for_regression$X2 = functional_data_eval_2[-test.idx,]
dataset_for_regression$Z = Z[-test.idx,]

dataset_for_prediction <- data.frame(y = response[test.idx])
dataset_for_prediction$X1 = functional_data_eval_1[test.idx,]
dataset_for_prediction$X2 = functional_data_eval_2[test.idx,]
dataset_for_prediction$Z = Z[test.idx,]

library(refund)
# first functional predictors
soft.fit <- pfr(
  y ~
    lf( # first functional predictor
      X1,
      argvals = eval_points,
      presmooth="bspline",
      presmooth.opts=list(nbasis=n_basis_predictor)
    ),
  data = dataset_for_regression
)


regression_result[1,1] <- mse(predict(soft.fit), dataset_for_regression$y)
regression_result[2,1] <- mse(predict(soft.fit, dataset_for_prediction), dataset_for_prediction$y)
regression_result

# first and second functional predictors
soft.fit <- pfr(
  y ~
    lf( # first functional predictor
      X1,
      argvals = eval_points,
      presmooth="bspline",
      presmooth.opts=list(nbasis=n_basis_predictor)
    ) +
    lf( # first functional predictor
      X2,
      argvals = eval_points,
      presmooth="bspline",
      presmooth.opts=list(nbasis=n_basis_predictor)
    ),
  data = dataset_for_regression
)

regression_result[1,2] <- mse(predict(soft.fit), dataset_for_regression$y)
regression_result[2,2] <- mse(predict(soft.fit, dataset_for_prediction), dataset_for_prediction$y)
regression_result

# all predictors
soft.fit <- pfr(
  y ~
    lf( # first functional predictor
      X1,
      argvals = eval_points,
      presmooth="bspline",
      presmooth.opts=list(nbasis=n_basis_predictor)
    ) +
    lf( # first functional predictor
      X2,
      argvals = eval_points,
      presmooth="bspline",
      presmooth.opts=list(nbasis=n_basis_predictor)
    ),
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
