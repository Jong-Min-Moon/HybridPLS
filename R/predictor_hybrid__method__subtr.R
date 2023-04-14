setMethod("subtr", "predictor_hybrid",
#######################################

function(input, other, alpha = 1){
  return(add(input, other, -alpha))
}

#######################################
)
