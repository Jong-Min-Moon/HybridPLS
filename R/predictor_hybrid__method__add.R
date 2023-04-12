setMethod("add", "predictor_hybrid",
###########################################

function(input, other, alpha = 1) {
# can add only when the bases are the same
  if(!(is_same_basis(input, other))){
    stop("Cannot add predictors with different bases")
  }else{
    # since the bases are the same, only update the coef of the bases
    for (i in 1:input@n_predictor_functional){
      input@predictor_functional_list[[i]]@coef <- add(
        input@predictor_functional_list[[i]]@coef,
        other@predictor_functional_list[[i]]@coef,
        alpha
      )
    }

    #update coef of the data
    input@Z <- add(input@Z, other@Z, alpha)
    return(input)
  }
  }

###########################################
)

