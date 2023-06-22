
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


my_basis <- create.bspline.basis(rangeval = c(0,1), nbasis = n_basis)

W_train_from_matrix <-fda::Data2fd(
  argvals = W_train_centered@predictor_functional_list[[1]]@original_t[1,],
  y = t(W_train_centered@predictor_functional_list[[1]]@original_X),
  basisobj = my_basis)
par(mfrow=c(1,2))
plot(W_train_fd[[1]])
plot(W_train_from_matrix)
par(mrfow=c(1,1))
