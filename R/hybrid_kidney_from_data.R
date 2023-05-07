library(Matrix)
library(MASS)
library(pracma)
library(jpeg)
library(fda)
library(splines2)

################# Read functions for kidney data ###################

hybrid_kidney_from_data <- function(
    argvals_base,
    argvals_post,
    reno_base,
    reno_post,
    scalar_predictors,
    n_basis
){
  cat(paste("\t For base curves, "))
  split_fit_base <- fit_spine_2d(argvals = argvals_base, evals = reno_base, n_basis = n_basis)

  cat(paste("\t For post curves, "))
  split_fit_post <- fit_spine_2d(argvals = argvals_post, evals = reno_post, n_basis = n_basis)

  predictor_functional_1 <- create_predictor_functional(
    split_fit_base$C,
    split_fit_base$J,
    split_fit_base$J_dotdot
    )

  predictor_functional_2 <- create_predictor_functional(
    split_fit_post$C,
    split_fit_post$J,
    split_fit_post$J_dotdot
  )

  hybrid_predictor_dataset <- create_hybrid_predictors_kidney(
    scalar_predictors,
    predictor_functional_1,
    predictor_functional_2
    )

  return(hybrid_predictor_dataset)
}
