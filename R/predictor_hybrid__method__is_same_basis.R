setMethod("is_same_basis", "predictor_hybrid",
###########################################

function(input, other){
  if(
    input@n_predictor_functional != other@n_predictor_functional
  ){
    stop("Number of the functional predictors must be the same")
  }

  is_all_elem_same = 1
  for (i in 1:input@n_predictor_functional){
    is_all_elem_same = is_all_elem_same * (
      prod(
        (input@predictor_functional_list[[i]]@J        == other@predictor_functional_list[[i]]@J) *
        (input@predictor_functional_list[[i]]@J_half   == other@predictor_functional_list[[i]]@J_half) *
        (input@predictor_functional_list[[i]]@J_dotdot == other@predictor_functional_list[[i]]@J_dotdot)
      )
    )
  }
if(is_all_elem_same){
  return(TRUE)
}else{
  return(FALSE)
  warning("Different bases")
  }
}

###########################################
)
