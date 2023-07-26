


# basis_hybridPLS object
setGeneric("is_same_basis", function(input, other)
  standardGeneric("is_same_basis"))

setGeneric("is_one_observation", function(input)
  standardGeneric("is_one_observation"))

setGeneric("get_sum_of_norm_sqrd", function(input)
  standardGeneric("get_sum_of_norm_sqrd"))


setGeneric("add", function(input, other, alpha=1)
  standardGeneric("add"))

setGeneric("add_broadcast", function(input, other, alpha=1)
  standardGeneric("add_broadcast"))

setGeneric(
  "subtr",
  function(input, other, alpha=1) standardGeneric("subtr")
  )

setGeneric(
  "subtr_broadcast",
  function(input, other, alpha=1)
    standardGeneric("subtr_broadcast")
  )

setGeneric(
  "scalar_mul",
  function(input, scalar) standardGeneric("scalar_mul")
  )

setGeneric(
  "get_mean",
  function(input) standardGeneric("get_mean")
)

setGeneric(
  "get_sd",
  function(input) standardGeneric("get_sd")
)



setGeneric(
  "hybrid_inner_prod",
  valueClass = "numeric",
  function(xi_1, xi_2) standardGeneric("hybrid_inner_prod")
  )

setGeneric(
  "LSE_ptws",
  function(input, rho) standardGeneric("LSE_ptws")
  )



setGeneric(
  "pls_pen",
  function(W, y, L) standardGeneric("pls_pen")
  )

setGeneric(
  "get_V_star",
  function(W, y) standardGeneric("get_V_star")
  )

setGeneric(
  "get_J_star",
  function(W) standardGeneric("get_J_star")
  )

setGeneric(
  "get_J_dotdot_star",
  function(W) standardGeneric("get_J_dotdot_star")
)

setGeneric(
  "get_Lambda",
  function(W, lambda) standardGeneric("get_Lambda")
)

setGeneric(
  "hybrid_from_coef",
  function(format, xi_star) standardGeneric("hybrid_from_coef")
  )

setGeneric(
  "fitted_value",
  function(data, coeff, pls_score) standardGeneric("fitted_value")
  )

setGeneric(
  "mean_norm",
  function(W) standardGeneric("mean_norm")
  )

setGeneric(
  "hybrid_norm",
  valueClass = "numeric",
  function(xi) standardGeneric("hybrid_norm")
  )

setGeneric(
  "index_sample",
  function(W, i) standardGeneric("index_sample")
  )


setGeneric(
  "extract_fd",
  function(input) standardGeneric("extract_fd")
)
