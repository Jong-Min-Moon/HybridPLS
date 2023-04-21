# matrix addition
A.mat <- matrix(rnorm(4), nrow = 2)
B.mat <- matrix(rnorm(4), nrow = 2)
print(A.mat)
print(B.mat)
print( add(A.mat, B.mat, 2) )


# create two functional predictors
functional_predictor_1 <- create_predictor_functional(coef = A.mat, J = B.mat %*% t(B.mat), J_dotdot = A.mat %*% t(A.mat))
functional_predictor_2 <- create_predictor_functional(coef = A.mat + 1, J = B.mat %*% t(B.mat) + 1, J_dotdot = A.mat %*% t(A.mat) + 1)

# create hybrid predictors with different bases
predictor_1 <- create_hybrid_predictors_kidney(A.mat, functional_predictor_1, functional_predictor_2)
predictor_2 <- create_hybrid_predictors_kidney(B.mat, functional_predictor_1, functional_predictor_1)
is_same_basis(predictor_1, predictor_2)
predictor_1@Z
predictor_2@Z
add(predictor_1, predictor_2)@Z #error occurs as expected

#this should work
predictor_1 <- create_hybrid_predictors_kidney(A.mat, functional_predictor_1, functional_predictor_2)
predictor_2 <- create_hybrid_predictors_kidney(B.mat, functional_predictor_1, functional_predictor_2)
is_same_basis(predictor_1, predictor_2)

predictor_1@Z
predictor_2@Z
add(predictor_1, predictor_2)@Z
A.mat + B.mat #verify

predictor_1@predictor_functional_list[[1]]
predictor_2@predictor_functional_list[[1]]
add(predictor_1, predictor_2)@predictor_functional_list[[1]]

predictor_1@predictor_functional_list[[1]]
predictor_2@predictor_functional_list[[1]]
sub(predictor_1, predictor_2)@predictor_functional_list[[1]]


