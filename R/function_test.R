# matrix addition
A.mat <- matrix(rnorm(4), nrow = 2)
B.mat <- matrix(rnorm(4), nrow = 2)
print(A.mat)
print(B.mat)
print( add(A.mat, B.mat, 2) )


# create two bases
base_1 <- create_basis_hybridPLS(coef = A.mat, J = B.mat %*% t(B.mat), J_dotdot = A.mat %*% t(A.mat))
base_2 <- create_basis_hybridPLS(coef = A.mat + 1, J = B.mat %*% t(B.mat) + 1, J_dotdot = A.mat %*% t(A.mat) + 1)

# basis is_same
is_same(base_1, base_2)
is_same(base_1, base_1)




# create hybrid predictors
predictor_1 <- create_hybrid_predictors_kidney(Z = A.mat, basis_1 = base_1, basis_2 = base_2)
predictor_2 <- create_hybrid_predictors_kidney(Z = B.mat, basis_1 = base_2, basis_2 = base_1)

# add
predictor_1@Z
predictor_2@Z
add(predictor_1, predictor_2)@Z #error occurs as expected

#this should work
predictor_1 <- create_hybrid_predictors_kidney(Z = A.mat, basis_1 = base_1, basis_2 = base_2)
predictor_2 <- create_hybrid_predictors_kidney(Z = B.mat, basis_1 = base_1, basis_2 = base_2)
predictor_1@Z
predictor_2@Z
add(predictor_1, predictor_2)@Z
predictor_1@basis_1
predictor_2@basis_1
add(predictor_1, predictor_2)@basis_1
