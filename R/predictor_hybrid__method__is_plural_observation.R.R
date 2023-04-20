setMethod("is_one_observation", "predictor_hybrid",
#####################################################

function(input){
  if(nrow(input@Z) == 1){
    return(TRUE)
  }else{
    cat("more than one observations")
    return(FALSE)
  }
}

###########################################
)
