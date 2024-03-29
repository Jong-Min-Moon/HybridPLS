#
# inputs:
#   W = centered predictor
#   y = response
#   L = number of iterations
#   lambda = hyperparameter
#   tau = hyperparmeter
#
# outputs:
# rho = pls score (=random variable created by inner producting)
# delta = predictor basis
# nu = response basis
# xi = pls component
nipals_pen_hybrid <- function(W, y, n_iter, lambda, tau) {

  # 1. initialize the storage
  rho <- delta <- nu <- xi <- sigma <- eta <- list()
  E <- V_star <- eigen_val <-list()
  fitted_value_W <- fitted_value_y <-list()
  resid_y <- W_now <- list() #data for iteration
  first_eigen_val <- mse_W <- mse_y <-rep(NA, n_iter)




  # 2. extract necessary numbers and matrices
  ## 2.1.  gram matrices are the same throughout the loop
  p <- dim(W@Z)[2] #number of scalar predictors
  J_star <- get_J_star(W)
  J_dotdot_star <- get_J_dotdot_star(W)

  ## 2.2. Create L_mat which imposes smoothness regularization to the objective function
  Lambda <- get_Lambda(W, lambda)
  J_Lambda_Jpp <-  J_star + (Lambda %*% J_dotdot_star) # denominator in the equation (6)
  L_mat <- t(chol(J_Lambda_Jpp))

  stopLoop = 0
  W_now[[1]] <- W
  y_now <- y
  for (l in 1:n_iter) {
    cat(paste("#############################################", "\n"))
    cat(paste(l, "th iteration", "\n"))

    # STEP 1. calculate pls score

    pls_pen_result <- pls_pen(W_now[[l]], y_now, L_mat)
    xi[[l]] <- pls_pen_result$xi # pls component. Section 3.2
    rho[[l]] <- hybrid_inner_prod(W_now[[l]], xi[[l]]) # n * 1 pls score. Step 1 in page 4.

    # STEP 2. residulized outcomes and predictors for the next step
    # update y
    nu[[l]] <- LSE_ptws(y_now, rho[[l]]) #hybrid residualization. Section 3.3
    fitted_value_y[[l]] <- fitted_value( y_now, nu[[l]] , rho[[l]])
    y_now <- y_now - fitted_value_y[[l]]
    resid_y[[l]] <- y_now
    # update Ws
    delta[[l]] <- LSE_hybrid(W_now[[l]], rho[[l]], tau) #regression coef (W on rho)
    fitted_value_W[[l]] <- fitted_value(W_now[[l]], delta[[l]], rho[[l]])
    W_now[[l+1]] <- subtr(W_now[[l]], fitted_value_W[[l]])

    #check mean norm
    mse_W[l] <- mean_norm(W_now[[l]])
    mse_y[l] <- Matrix::norm(y_now,"2")

    cat("W:", mse_W[l], "\n")
    cat("y:",  mse_y[l], "\n\n")

    # STEP 3.
    sigma[[l]] <- xi[[l]]
    if (l == 1) {
      eta[[l]] <- scalar_mul(sigma[[l]], nu[[l]])
    }else{ # if l > 1
      for (u in 1:(l - 1)){
        sigma[[l]] <- subtr( sigma[[l]], sigma[[u]], hybrid_inner_prod(delta[[u]], xi[[l]]) )
      } #end of u-loop
      eta[[l]] <- add(eta[[l-1]], sigma[[l]], nu[[l]])
    }#if statement

    # just for records.
    E[[l]] <- pls_pen_result$E # for monitoring symetricity
    V_star[[l]] <- pls_pen_result$V_star
    eigen_val[[l]] <- pls_pen_result$eigen_val
    first_eigen_val[l] <- eigen_val[[l]][1]

  }# for loop: l



pls_object <- new(
  "hybrid_pls_result",
                    eta = eta,
                    xi = xi,
                    nu = nu,
                    rho = rho,
                    delta = delta,
                    L_mat = L_mat,
                    V_star = V_star,
                    E = E,
                    eigen_val = eigen_val,
                    first_eigen_val = first_eigen_val,
                    J_Lambda_Jpp = J_Lambda_Jpp,
                    mse_W = mse_W,
                    mse_y = mse_y,
                    resid_y = resid_y,
                    fitted_value_W = fitted_value_W,
                    fitted_value_y = fitted_value_y,
                    W_now = W_now
)
  return(pls_object)
  cat(paste("#############################################", "\n"))
}#end of the function

setClass(
  "hybrid_pls_result",
  representation(
    eta = "list", #estimated regression coefficient function
    xi = "list",
    nu = "list",
    rho = "list",
    delta = "list",
    E = "list",
    V_star = "list",
    eigen_val = "list",
    first_eigen_val = "vector",
    J_Lambda_Jpp = "dgCMatrix",
    L_mat= "matrix",
    resid_y = "list",
    mse_W = 'vector',
    mse_y = 'vector',
    fitted_value_W = "list",
    fitted_value_y = "list",
    W_now = "list"
  ))

setMethod("predict", signature(object = "hybrid_pls_result"),
          function(object, ...) {
            max_iter = length(eta)
            dots <- list(...)
            if (length(dots) == 0){
              return(
                hybrid_inner_prod(
                  (object@W_now)[[1]], (object@eta)[[max_iter]])
              )
            }else{
              new_data <- dots[[1]]
              iter <- dots[[2]]
              return(
                hybrid_inner_prod(new_data, (object@eta)[[iter]])
              )
            }
          }
)
