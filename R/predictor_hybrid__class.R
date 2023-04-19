# n: number of observations
# p: number of scalar predictors
# M: number of basis functions for each functional predictors

setClass(
  "predictor_hybrid",
  representation(
    # scalar predictors:
    ## n x p matrix of real-valued covariates.
    ## each row represents an observed finite-dimensional vector.
    Z = "matrix",

    # functional predictors
    ## basis expansion representation
    ## number of predictors can be arbitrary
    predictor_functional_list = "list",
    n_predictor_functional = "numeric",

    # basic information
    n_sample = "numeric"
  ))

create_hybrid_predictors_kidney <- function(Z, predictor_functional_1, predictor_functional_2){
if(
  !(
    class(predictor_functional_1)[1] == "predictor_functional"
    ) & (
      class(predictor_functional_2)[1] == "predictor_functional"
      )
){
  stop("Inputs should be functional predictor objects")
}else{
  predictor_object <- new("predictor_hybrid",
                          Z = Z,
                          predictor_functional_list = list(
                            predictor_functional_1,
                            predictor_functional_2
                            ),
                          n_predictor_functional = 2,
                          n_sample = dim(Z)[1]
                          )
  return(predictor_object)
  }
}
