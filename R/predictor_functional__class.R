
setClass(
  "predictor_functional",
  representation(
    #for basis expansion representation
    coef = "matrix", # n x M matrix of basis coefficients. each row represents an observed function.

    #gram matrices. Bases are only used in the form of Gram matrix.
    J = "matrix", #original gram matrix (M x M)
    J_half = "matrix", # matrix square root of the gram matrix (M x M)
    J_dotdot = "matrix", # gram matrix for second derivative (M x M)
    original_t = "matrix",
    original_X = "matrix"
  ))

create_predictor_functional <- function(coef, J, J_dotdot, original_t, original_X){
  predictor_object <- new("predictor_functional",
                          coef = coef,
                          J = J,
                          J_half = get_gram_half(J),
                          J_dotdot = J_dotdot,
                          original_t = original_t,
                          original_X = original_X
                          )
  return(predictor_object)
}

setMethod("scalar_mul", signature("predictor_functional", "numeric"),
###########################################
function(input, scalar){
  result <- input
  result@coef <- scalar * (input@coef)
  result@original_X <- scalar * (input@original_X)
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
      result@original_X <-add_broadcast(input@original_X, other@original_X, alpha)
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

setMethod("get_mean", "predictor_functional",
          ###########################################
          function(input){
            mean_object <- input
            mean_object@coef <- matrix(apply(mean_object@coef, 2, mean), nrow=1)
            mean_object@original_X <- matrix(apply(mean_object@original_X, 2, mean), nrow=1)
            return(mean_object)
          }
          ###########################################
)

setMethod("get_sum_of_norm_sqrd", "predictor_functional",
          #output: a scalar value
###########################################
function(input){
  n <- nrow(input@coef) #sample size
  norm_sqrd <- 0
  for (i in 1:n){
    coef_now <- matrix(input@coef[i,], nrow = 1)
    print(coef_now)
    print(dim(input@J))
    norm_sqrd <- norm_sqrd + sum(
      (t(coef_now) %*% coef_now) * input@J
    )
    }
  return(norm_sqrd)
  }
###########################################
)


setMethod("plot", signature("predictor_functional"),
          function(x){

            # Plot using base graphic engine
            n_basis <- dim(x@predictor_functional_list[[1]]@J)[2]
            splines <- create.bspline.basis(rangeval = c(0,1), nbasis = n_basis)
            first <-fd(coef = t(x@predictor_functional_list[[1]]@coef), basisobj = splines)
            second <-fd(coef = t(x@predictor_functional_list[[2]]@coef), basisobj = splines)

            par(mfrow=c(2,1))
            plot(first)
            plot(second)
            par(mfrow=c(1,1))
            #y_1 <- eval.fd(x, psi_1)
            #print(y_1)
            #plot(x, y_1,  pch=".")
            #plot(x, eval.fd(x, psi_2))
          })


# get_Lambda
#
# input:
# W
# lambda: K-dimensional vector
#
# transition from scalar kappa to vector lambda (created October 27th 2022)
#

