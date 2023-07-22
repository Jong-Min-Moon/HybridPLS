get_traning_idx <- function(n_sample, training_ratio){
  return( sample(n_sample, floor(training_ratio * n_sample)) )
}

get_gram_2d <- function(argval, basis_evals){
  n_basis <- dim(basis_evals)[2]
  J <- matrix(nrow= n_basis, ncol = n_basis)
  for (i in 1:n_basis) {
    for (j in 1:n_basis){
      J[i,j] <- pracma::trapz(argval, basis_evals[,i]*basis_evals[,j])
    }}
  return(J)
}

get_gram_half <- function(J) {
  # gram matrix is positive definite,
  # so eigenvalues can be square-rooted, and eigenvectors are orthonormal.
  J_eigen <- eigen(J)
  S_half <- diag(sqrt(J_eigen$values))
  V <- J_eigen$vectors
  J_half <- V %*% S_half %*% t(V)
  return(J_half)
}

fit_spine_2d <- function(argvals, evals, n_basis){
  cat(paste("use", n_basis, "basis functions\n"))
  # create 1d b-spline basis functions
  #PhiB <- splines2::bSpline(argvals, df=n_basis, degree = 3, intercept = TRUE) #Jt(n.eval) X M B spline Basis #old way
  my_basis <- create.bspline.basis(rangeval = c(0,1), nbasis = n_basis)


  # calculate basis coefficients for observed functional data # idea from Dr. Jang
  #C <- t( MASS::ginv(PhiB) %*% t(evals) ) #old way
  fda_obj <- Data2fd(argvals, t(evals), my_basis)
  C <- t(fda_obj$coefs)
  # calculate gram matrices for
  # - basis functions (J)
  # - 2nd order derivative of basis functions (J'')
  J <- inprod(my_basis, my_basis)
  J_dotdot <- getbasispenalty(my_basis)

  return(
    list(
      C = C,
      J = J,
      J_dotdot = J_dotdot,
      fd = fda_obj
    )
  )
}

dataset_for_regression_refund <- function(response, hybrid_predictor){
  dataset_for_regression <- data.frame(hybrid_predictor@Z)
  dataset_for_regression$y <- response
  dataset_for_regression$F1 <-hybrid_predictor@predictor_functional_list[[1]]@original_X
  dataset_for_regression$F2 <-hybrid_predictor@predictor_functional_list[[2]]@original_X
  return(dataset_for_regression)
}

dataset_for_prediction_refund <- function(hybrid_predictor){
  dataset_for_regression <- data.frame(hybrid_predictor@Z)
  dataset_for_regression$F1 <-hybrid_predictor@predictor_functional_list[[1]]@original_X
  dataset_for_regression$F2 <-hybrid_predictor@predictor_functional_list[[2]]@original_X
  return(dataset_for_regression)
}

dataset_for_PCA <- function(dataset_for_refund, eval_point){
  n = nrow(dataset_for_refund)
  F1 <- F2 <- eval_point_PC <- list()
  for (i in 1:n){
    F1[[i]] <- (dataset_for_refund$F1)[i,]
    F2[[i]] <- (dataset_for_refund$F2)[i,]
    eval_point_PC[[i]] <- eval_point
  }
  n_eval <- length(F1[[1]])
  total_dim <- ncol(dataset_for_refund)
  scalar <- dataset_for_refund[1: (total_dim-3)]

  return(
    list(
      "F1" = F1, "F2" = F2, "scalar" = scalar, "eval_point" = eval_point_PC
    )
  )
}
