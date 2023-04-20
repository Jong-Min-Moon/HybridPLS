# basis_hybridPLS object
setGeneric("is_same_basis", function(input, other)
  standardGeneric("is_same_basis"))

setGeneric("get_norm_sqrd", function(input)
  standardGeneric("get_norm_sqrd"))


setGeneric("add", function(input, other, alpha=1)
  standardGeneric("add"))

setGeneric(
  "subtr",
  function(input, other, alpha=1) standardGeneric("subtr")
  )


setGeneric(
  "pls_pen",
  function(W, y, L) standardGeneric("pls_pen")
  )

setGeneric(
  "get_V_star",
  function(W, y) standardGeneric("get_V_star")
  )


