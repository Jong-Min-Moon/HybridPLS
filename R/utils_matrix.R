setMethod("add", "matrix",
###########################################
function(input, other, alpha = 1) {
  # dimension check
  if(
    prod( dim(input) != dim(other) ) > 0
  ){
    stop("Cannot add")
  # operation
  }else{
    return(input + alpha * other)
  }
  }
###########################################
)

setMethod("add_broadcast", "matrix",
###########################################
function(input, other, alpha = 1){
  if(is.vector(other)){
    other <- matrix(other, nrow = 1)
  }else if(dim(other)[1] >1){
      (stop("RHS must have single observation"))
    }
  return(input + rep(alpha, nrow(input) ) %*%  other)
}
###########################################
)

setMethod("subtr_broadcast", "matrix",
###########################################
function(input, other, alpha = 1){
  return(add_broadcast(input, other, (-1*alpha)))
  }
###########################################
)

setMethod("get_mean", "matrix",
###########################################
function(input){
  return(matrix(apply(input, 2, mean), nrow=1))
  }
###########################################
)

setMethod("get_sd", "matrix",
###########################################
function(input){
  return(
    matrix(apply(input, 2, sd), nrow=1)
    )
  }
###########################################
)
