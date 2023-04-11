setClass(
  "basis_hybridPLS",
  representation(
    #for basis expansion representation
    coef = "matrix", # n x M matrix of basis coefficients. each row represents an observed function.

    #gram matrices. Bases are only used in the form of Gram matrix.
    J = "matrix", #original gram matrix
    J_half = "matrix", # matrix square root of the gram matrix
    J_dotdot = "matrix" # gram matrix for second derivative
  ))
