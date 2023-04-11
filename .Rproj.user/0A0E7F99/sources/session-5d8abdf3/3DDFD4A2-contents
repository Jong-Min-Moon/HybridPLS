setMethod("add", "hybrid_predictors_kidney",
###########################################

function(input, other, alpha = 1) {
# basis check
  if(
    is_same(input@basis_1, input@basis_2)
  ){
    stop("Cannot add")
  }else{
    (input@basis_1)@coeff <- add((input@basis_1)@coeff, (other@basis_1)@coeff, alpha)
    (input@basis_2)@coeff <- add((input@basis_2)@coeff, (other@basis_2)@coeff, alpha)
    input@Z <- add(input@Z, other@Z, alpha)
    return(input)
  }
  }

###########################################
)
