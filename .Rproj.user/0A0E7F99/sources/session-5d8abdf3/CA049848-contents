setMethod("is_same", "basis_hybridPLS",
###########################################

function(input, other) {
            # basis check
if (
  prod( (input@basis_1)@J != (other@basis_1)@J ) * prod(input@J_2 != other@J_2) > 0
){
  return(FALSE)
  warning("Different bases")
}else{
  return(TRUE)
  }
}
###########################################
)
