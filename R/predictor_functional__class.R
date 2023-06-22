
setClass(
  "predictor_functional",
  representation(
    #for basis expansion representation
    coef = "matrix", # n x M matrix of basis coefficients. each row represents an observed function.

    #gram matrices. Bases are only used in the form of Gram matrix.
    J = "matrix", #original gram matrix (M x M)
    J_half = "matrix", # matrix square root of the gram matrix (M x M)
    J_dotdot = "matrix", # gram matrix for second derivative (M x M)
    original_t = "matrix"
  ))

create_predictor_functional <- function(coef, J, J_dotdot, original_t){
  predictor_object <- new("predictor_functional",
                          coef = coef,
                          J = J,
                          J_half = get_gram_half(J),
                          J_dotdot = J_dotdot,
                          original_t = original_t
                          )
  return(predictor_object)
}

setMethod("scalar_mul", signature("predictor_functional", "numeric"),
###########################################
function(input, scalar){
  result <- input
  result@coef <- scalar * (input@coef)
  return(result)
}
###########################################
)
setMethod("add_broadcast", "predictor_functional",
###########################################
function(input, other, alpha = 1) {
  if(prod(input@J != other@J)>0){
  stop("In order to add two predictor_functional objects, they should use the same basis")
    }else{
      result <- input
      result@coef <-add_broadcast(input@coef, other@coef, alpha)
      return(result)
    }
  }
###########################################
)

setMethod("subtr_broadcast", "predictor_functional",
###########################################
function(input, other, alpha = 1) {
###########################################
  return(add_broadcast(input, other, (-1 * alpha)))
  }
###########################################
)
setMethod("get_sum_of_norm_sqrd", "predictor_functional",
###########################################
function(input){
  n <- nrow(input@coef) #sample size
  norm_sqrd <- 0
  for (i in 1:n){
    coef_now <- matrix(input@coef[i,], nrow = 1)
    norm_sqrd <- norm_sqrd + sum(
      (t(coef_now) %*% coef_now) * input@J
    )
    }
  return(norm_sqrd)
  }
###########################################
)

setMethod("get_mean", "predictor_functional",
###########################################
function(input){
  mean_object <- input
  mean_object@coef <- matrix(apply(mean_object@coef, 2, mean), nrow=1)
  return(mean_object)
  }
###########################################
)



# get_Lambda
#
# input:
# W
# lambda: K-dimensional vector
#
# transition from scalar kappa to vector lambda (created October 27th 2022)
#

