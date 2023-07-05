my_basis <- create.bspline.basis(rangeval = c(0,1), nbasis = n_basis)
PhiB <- predict(my_basis, eval_point)
PhiB <- PhiB[nrow(PhiB):1, ncol(PhiB):1]

# calculate basis coefficients for observed functional data # idea from Dr. Jang
#C <- t( MASS::ginv(PhiB) %*% t(evals) ) #old way
fda_obj <- Data2fd(argvals, t(evals), my_basis)
C <- t(fda_obj$coefs)
# calculate gram matrices for
# - basis functions (J)
# - 2nd order derivative of basis functions (J'')
J <- get_gram_2d(eval_point, PhiB)

#J_dotdot <- get_gram_2d(argvals, deriv(PhiB, 2))
deriv_PhiB <- predict(my_basis, eval_point, deriv=2)
deriv_PhiB <- deriv_PhiB[nrow(deriv_PhiB):1, ncol(deriv_PhiB):1]
J_dotdot <- get_gram_2d(eval_point, deriv_PhiB)


deriv.fd(my_basis,2)
#  set up a B-spline basis of order 4 with 13 basis functions
#  and knots at 0.0, 0.1,..., 0.9, 1.0.
basisobj <- create.bspline.basis(c(0,1),13)
#  compute the 13 by 13 matrix of inner products of second derivatives
penmat <- getbasispenalty(basisobj)




M=10
PhiB <- bSpline(eval_point, df=M, degree = 3, intercept = TRUE)
PhiB_d2 <-  deriv(PhiB, 2) # Jt X M second derivatibe B-splines


for (i in 1:M) {
  for (j in 1:M)
    J[i,j] <- trapz(eval_point,PhiB[,i]*PhiB[,j])
}
head(J)
head(inprod(my_basis,my_basis))

for (i in 1:M) {
  for (j in 1:M)
    J_dotdot[i,j] <- trapz(eval_point,PhiB_d2[,i]*PhiB_d2[,j])
}
head(J_dotdot)
head(getbasispenalty(my_basis))
J_dotdot - getbasispenalty(my_basis)
