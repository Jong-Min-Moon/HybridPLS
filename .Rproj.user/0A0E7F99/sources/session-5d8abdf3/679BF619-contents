setClass(
  "hybrid_predictors_kidney",
  representation(
    #for basis expansion representation
    Z = "matrix", # n x p matrix of real-valued covariates. each row represents an observed finite-dimensional vector.
    basis_1 = "basis_hybridPLS",
    basis_2 = "basis_hybridPLS"
  ))

create_hybrid_predictors_kidney <- function(Z, basis_1, basis_2){
  predictor_object <- new("hybrid_predictors_kidney",
                          Z = Z,
                          basis_1 = basis_1,
                          basis_2 = basis_2
  )
  return(predictor_object)
}
