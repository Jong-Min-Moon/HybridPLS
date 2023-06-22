#1. curve normalization 구현
#2

y_limit = c(0.04, 0.15)
mse <- rep(NA, 5)
names(mse) <- c("hybridPLS", "sofreg", "FPCA+PCAreg", "PLS", "PCA+PLS")


# 1. HybridPLS

seed = 1
set.seed(seed) # for training/test split
n_basis = 10
#set the number of max iterations
L_max = 10
lambda = c(1e-6,5*1e-6)

kidney_obj <- read_fd_kidney(
  test_ratio = 0.3,
  n_basis = n_basis,
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
  #learn the model
  pls_model_transform <- nipals_pen_hybrid(
    W_train_centered, y_train_logit_centered, L, lambda, 0)

  # predict
  y_pred_pls <- predict_test(pls_model_transform, W_test_centered)

  #inverse transform
  y_pred_pls <- reponse_inverse_transform_min_max_logit(
    y_pred_pls, y_train_logit_mean, y_train_max, y_train_min)

  # calculate the mse and save
  L_trend_0_logit[L] <- mse(y_pred_pls, y_test)
}

plot(L_trend_0_logit[1:L_max], main ="hybridPLS", ylab = "MSE", xlab = "number of iteration", pch=".", ylim=y_limit)
lines(L_trend_0_logit)
mse[1] <- min(L_trend_0_logit, omit.NA = TRUE)
mse[1]


##############################################################################################
# 2. scalar-on-function regression with scalar covariate


timepoints <- W_train_centered@predictor_functional_list[[1]]@original_t
dataset_for_regression <- dataset_for_regression_Emory(y_train_logit_centered, W_train_centered)
dataset_for_prediction <- dataset_for_prediction_Emory(W_test_centered)
regression_result <- matrix(NA, nrow = 2, ncol = 3)
rownames(regression_result) <- c("train", "test")
colnames(regression_result) <- c("first", "second", "all")




# first functional predictors
soft.fit <- pfr(
  y_train_logit_centered ~
    lf(precurve, k=10, argvals = W_train_centered@predictor_functional_list[[1]]@original_t, bs="bs") +
    X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + X11 + X12 + X13 + X14 + X15,
  data = dataset_for_regression
)

y_pred_pls <- reponse_inverse_transform_min_max_logit(
  predict(soft.fit), y_train_logit_mean, y_train_max, y_train_min) # training fit
regression_result[1,1] <- mse(y_pred_pls, y_train)

y_pred_pls <- reponse_inverse_transform_min_max_logit(
  predict(soft.fit, dataset_for_prediction), y_train_logit_mean, y_train_max, y_train_min) # test fit
regression_result[2,1] <- mse(y_pred_pls, y_test)


# second functional predictors
soft.fit <- pfr(
  y_train_logit_centered ~
    lf(postcurve, k=10, argvals = W_train_centered@predictor_functional_list[[2]]@original_t, bs="bs") +
    X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + X11 + X12 + X13 + X14 + X15,
  data = dataset_for_regression
)

y_pred_pls <- reponse_inverse_transform_min_max_logit(
  predict(soft.fit), y_train_logit_mean, y_train_max, y_train_min) # training fit
regression_result[1,2] <- mse(y_pred_pls, y_train)

y_pred_pls <- reponse_inverse_transform_min_max_logit(
  predict(soft.fit, dataset_for_prediction), y_train_logit_mean, y_train_max, y_train_min) # test fit
regression_result[2,2] <- mse(y_pred_pls, y_test)

# all functional predictors
soft.fit <- pfr(
  y_train_logit_centered ~
    lf(precurve, k=10, argvals = W_train_centered@predictor_functional_list[[1]]@original_t, bs="bs") +
    lf(postcurve, k=10, argvals = W_train_centered@predictor_functional_list[[2]]@original_t, bs="bs") +
    X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + X11 + X12 + X13 + X14 + X15,
  data = dataset_for_regression
)

y_pred_pls <- reponse_inverse_transform_min_max_logit(
  predict(soft.fit), y_train_logit_mean, y_train_max, y_train_min) # training fit
regression_result[1,3] <- mse(y_pred_pls, y_train)

y_pred_pls <- reponse_inverse_transform_min_max_logit(
  predict(soft.fit, dataset_for_prediction), y_train_logit_mean, y_train_max, y_train_min) # test fit
regression_result[2,3] <- mse(y_pred_pls, y_test)

regression_result

# FPCA + PCA




# FPCA predictor 1
t_point_base <- (W_train_centered@predictor_functional_list[[1]]@original_t)[1,]
base.pcalist = pca.fd(W_train_fd[[1]], n_basis)
print(base.pcalist$values)
## train data
discrete_data_base_train <- t( eval.fd(W_train_fd[[1]], t_point_base) ) #evalauation of the normalized train data at the observed timepoints
PhiB_FPCA_base <- eval.fd(base.pcalist$harmonics, t_point_base) # evalauation of the FPCA components(eigenfunctions) at the observed timepoints
dim(PhiB_FPCA_base)
FPCA_score_base_train <- t( MASS::ginv(PhiB_FPCA_base) %*% t(discrete_data_base_train) ) #FPCA score, base, train data
FPCA_score_base_train <- FPCA_score_base_train[,1:2]
## test data
discrete_data_base_test <- t( eval.fd(W_test_fd[[1]], t_point_base) )
FPCA_score_base_test <- t( MASS::ginv(PhiB_FPCA_base) %*% t(discrete_data_base_test) )
FPCA_score_base_test <- FPCA_score_base_test[,1:2]

# FPCA predictor 2
t_point_post <- (W_train_centered@predictor_functional_list[[2]]@original_t)[1,]
post.pcalist = pca.fd(W_train_fd[[2]], n_basis)
print(post.pcalist$values)
## train data
discrete_data_post_train <- t( eval.fd(W_train_fd[[2]], t_point_post) ) #evalauation of the normalized train data at the observed timepoints
PhiB_FPCA_post <- eval.fd(post.pcalist$harmonics, t_point_post) # evalauation of the FPCA components(eigenfunctions) at the observed timepoints
dim(PhiB_FPCA_post)
FPCA_score_post_train <- t( MASS::ginv(PhiB_FPCA_post) %*% t(discrete_data_post_train) ) #FPCA score, base, train data
FPCA_score_post_train <- FPCA_score_post_train[,1:2]
## test data
discrete_data_post_test <- t( eval.fd(W_test_fd[[1]], t_point_post) )
FPCA_score_post_test <- t( MASS::ginv(PhiB_FPCA_post) %*% t(discrete_data_post_test) )
FPCA_score_post_test <- FPCA_score_post_test[,1:2]

PCA_scalar <- prcomp(W_train_centered@Z) #learn PCA with train data
summary(PCA_scalar)
PCA_score_train <- (PCA_scalar$x)[,1:7]
PCA_score_test <- (predict(PCA_scalar, W_test_centered@Z))[,1:7]

#combine
PCA_data_train <- data.frame(cbind(y_train_logit_centered, FPCA_score_base_train, FPCA_score_post_train, PCA_score_train))
PCA_data_test <- data.frame(cbind(FPCA_score_base_test, FPCA_score_post_test, PCA_score_test))
colnames(PCA_data_test) <- colnames(PCA_data_train)[-1]
#regress
pcareg <- lm(y_train_logit_centered~., data = PCA_data_train)
#inverse transform
y_pred_pls <- predict(pcareg, PCA_data_test)
y_pred_pls <- y_pred_pls + y_train_logit_mean
y_pred_pls <- exp(y_pred_pls)/( 1 + exp(y_pred_pls))
y_pred_pls <- y_pred_pls * (y_train_max- y_train_min)
y_pred_pls <- y_pred_pls + y_train_min
mse[3] <- mse(y_pred_pls, y_test)




#############
# 4. PLS
PLS_data_train <- data.frame(cbind(y_train_logit_centered, discrete_data_base_train, discrete_data_post_train, W_train_centered@Z))
PLS_data_test <- data.frame(cbind(discrete_data_base_test, discrete_data_post_test, W_test_centered@Z))
colnames(PLS_data_test) <- colnames(PLS_data_train)[-1]
plsreg <- plsr(y_train_logit_centered ~ .,ncomp=1, data = PLS_data_train)

#inverse transform
y_pred_pls <- predict(plsreg, PLS_data_test)
y_pred_pls <- y_pred_pls + y_train_logit_mean
y_pred_pls <- exp(y_pred_pls)/( 1 + exp(y_pred_pls))
y_pred_pls <- y_pred_pls * (y_train_max- y_train_min)
y_pred_pls <- y_pred_pls + y_train_min
mse[4] <- mse(y_pred_pls, y_test)

#############
# 5. PCA+PLS
plsreg <- plsr(y_train_logit_centered ~ .,ncomp=1, data = PCA_data_train)

#inverse transform
y_pred_pls <- predict(plsreg, PCA_data_test)
y_pred_pls <- y_pred_pls + y_train_logit_mean
y_pred_pls <- exp(y_pred_pls)/( 1 + exp(y_pred_pls))
y_pred_pls <- y_pred_pls * (y_train_max- y_train_min)
y_pred_pls <- y_pred_pls + y_train_min
mse[5] <- mse(y_pred_pls, y_test)
mse






#4. functional PLS using only the functional predictors
W_train_fdata <- list()
W_train_fdata[[1]] <- fda.usc::fdata(W_train_fd[[1]])
W_train_fdata[[2]] <- fda.usc::fdata(W_train_fd[[2]])

W_test_fdata <- list()
W_test_fdata[[1]] <- fda.usc::fdata(W_test_fd[[1]])
W_test_fdata[[2]] <- fda.usc::fdata(W_test_fd[[2]])

res <- fda.usc::fregre.pls(W_train_fdata[[2]], y_train_logit_centered)
