
# 1. HybridPLS

seed = 1
set.seed(seed) # for training/test split
n_basis = 10
#set the number of max iterations


kidney_obj <- read_fd_kidney(
  test_ratio = 0.3,
  n_basis = n_basis,
  normalize = list("curve" = TRUE, "scalar"= TRUE, "between" = TRUE)
)


# W is already centered
W_train_centered <- kidney_obj$W_train
y_train <- kidney_obj$y_train


W_train_fd <- extract_fd(W_train_centered)


timepoints <- W_train_centered@predictor_functional_list[[1]]@original_t


dataset_for_regression <- dataset_for_regression_Emory(y_train_logit_centered, W_train_centered)
soft.fit <- pfr(
  y_train_logit_centered ~ lf(precurve,
                              k=10,
                              argvals = W_train_centered@predictor_functional_list[[1]]@original_t),
  data = dataset_for_regression)
plot(soft.fit, select = 1)

dataset_for_prediction <- dataset_for_prediction_Emory(W_test_centered)


#inverse transform
y_pred_pls <- predict(soft.fit, dataset_for_prediction)
y_pred_pls <- reponse_inverse_transform_min_max_logit(
  y_pred_pls, y_train_logit_mean, y_train_max, y_train_min)
mse[2] <- mse(y_pred_pls, y_test)
mse[2]

