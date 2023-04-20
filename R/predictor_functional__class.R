setClass(
  "predictor_functional",
  representation(
    #for basis expansion representation
    coef = "matrix", # n x M matrix of basis coefficients. each row represents an observed function.

    #gram matrices. Bases are only used in the form of Gram matrix.
    J = "matrix", #original gram matrix (M x M)
    J_half = "matrix", # matrix square root of the gram matrix (M x M)
    J_dotdot = "matrix" # gram matrix for second derivative (M x M)
  ))

create_predictor_functional <- function(coef, J, J_dotdot){
  predictor_object <- new("predictor_functional", coef = coef, J = J, J_half = .get_J_half(J), J_dotdot = J_dotdot)
  return(predictor_object)
}
