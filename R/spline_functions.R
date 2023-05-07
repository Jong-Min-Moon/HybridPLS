fit_spine_2d <- function(argvals, evals, n_basis){

  # create 1d b-spline basis functions
  PhiB <- splines2::bSpline(argvals, df=n_basis, degree = 3, intercept = TRUE) #Jt X M B spline Basis
  cat(paste("use", n_basis, "basis functions\n"))

  # calculate basis coefficients for observed functional data
  # idea from Dr. Jang
  C <- t(MASS::ginv(PhiB) %*% t(evals))

  # calculate gram matrices for
  # - basis functions (J)
  # - 2nd order derivative of basis functions (J'')
  J <- get_gram(argvals, PhiB)
  J_dotdot <- get_gram(argvals, deriv(PhiB, 2))

  return(
    list(
      C = C,
      J = J,
      J_dotdot = J_dotdot
    )
  )
}
