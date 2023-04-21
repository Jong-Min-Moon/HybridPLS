library(Matrix)
library(MASS)
library(pracma)
library(jpeg)
library(fda)
library(splines2)


# n: number of observations
# p: number of scalar predictors
# M: number of basis functions for each functional predictors

setClass(
  "predictor_hybrid",
  representation(
    # scalar predictors:
    ## n x p matrix of real-valued covariates.
    ## each row represents an observed finite-dimensional vector.
    Z = "matrix",

    # functional predictors
    ## basis expansion representation
    ## number of predictors can be arbitrary
    predictor_functional_list = "list",
    n_predictor_functional = "numeric",

    # basic information
    n_sample = "numeric"
  ))

create_hybrid_predictors_kidney <- function(Z, predictor_functional_1, predictor_functional_2){
if(
  !(
    class(predictor_functional_1)[1] == "predictor_functional"
    ) & (
      class(predictor_functional_2)[1] == "predictor_functional"
      )
){
  stop("Inputs should be functional predictor objects")
}else{
  predictor_object <- new("predictor_hybrid",
                          Z = Z,
                          predictor_functional_list = list(
                            predictor_functional_1,
                            predictor_functional_2
                            ),
                          n_predictor_functional = 2,
                          n_sample = dim(Z)[1]
                          )
  return(predictor_object)
  }
}



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

setMethod("subtr", "predictor_hybrid",
          #######################################

          function(input, other, alpha = 1){
            return(add(input, other, -alpha))
          }

          #######################################
)


setMethod("scalar_mul", signature("predictor_hybrid", "numeric"),
          ################################################################
          function(W, scalar) {
            for (i in 1:W@n_predictor_functional){
              W@predictor_functional_list[[i]]@coef <- scalar * W@predictor_functional_list[[i]]@coef
            }
            W@Z <- scalar*(W@Z)
            return(W)
          }
          ################################################################
)


setMethod("hybrid_inner_prod", signature("predictor_hybrid", "predictor_hybrid"),
          #####################################################################

          function(xi_1, xi_2) {
            if (is_one_observation(xi_2) &
                is_same_basis(xi_1, xi_2)
            ){
              ip_functional <- 0
              for (i in 1:(xi_1@n_predictor_functional)){
                coeff_1 <- xi_1@predictor_functional_list[[i]]@coef
                coeff_2 <- xi_2@predictor_functional_list[[i]]@coef
                J <- xi_1@predictor_functional_list[[i]]@J
                ip_functional <- ip_functional + as.numeric(coeff_1  %*% J  %*% t(coeff_2))
              }
              ip_scalar <- as.numeric(xi_1@Z %*% t(xi_2@Z))
              return(ip_functional + ip_scalar)
            }else{
              stop("No observation or different bases")
            }

          }
          ###################################################################################################
)



setMethod("hybrid_norm",
          signature("predictor_hybrid"),
          function(xi) {
            return( sqrt(hybrid_inner_prod(xi, xi)) )}
)


setMethod("mean_norm", signature("predictor_hybrid"),
          ######################################################
          function(W){
            norm_W = 0
            n = dim(W@Z)[1]
            print(n)
            for (i in 1:n){
              w_now <- index_sample(W, i)
              norm_W <- norm_W + hybrid_norm(w_now)
            }
            return(norm_W/n)
          }

          ######################################################
)


setGeneric(
  "LSE_hybrid",
  function(W, rho, tau) standardGeneric("LSE_hybrid")
)

setMethod("LSE_hybrid", "predictor_hybrid",
          ###########################################

          function(W, rho, tau) {
            p <- dim(W@Z)[2]

            rho_t_C_star_J_star <- (W@predictor_functional_list[[1]]@coef) %*% (W@predictor_functional_list[[1]]@J)
            for (i in (2 : W@n_predictor_functional)){
              rho_t_C_star_J_star <- cbind(
                rho_t_C_star_J_star,
                (W@predictor_functional_list[[i]]@coef) %*% (W@predictor_functional_list[[i]]@J)
              )
            }
            rho_t_C_star_J_star <- t(rho) %*% cbind(rho_t_C_star_J_star, W@Z)


            # matrices in the statement of Proposition 4
            J_star <-get_J_star(W)
            J_dotdot_star <- get_J_dotdot_star(W)
            convex_J <- (norm(rho, type="2"))^2 * J_star + tau * J_dotdot_star
            d_hat_star <- t(
              solve( t(convex_J), t(rho_t_C_star_J_star) )
            )
            delta_hat <- hybrid_from_coef(W, d_hat_star)
            return(delta_hat)
          }

          #########################################
)



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



setMethod("index_sample", signature("predictor_hybrid"),
          ################################################################
          function(W, i){
            n <- dim(W@Z)[1]
            if(i > n){stop("index out of range")}
            W@Z <- matrix(W@Z[i,], nrow = 1)
            for (j in 1 : W@n_predictor_functional){
              W@predictor_functional_list[[j]]@coef <- matrix(W@predictor_functional_list[[j]]@coef[i,], nrow = 1)
            }
            return(W)
          }

          ################################################################
)



setMethod("fitted_value",
          signature(
            "predictor_hybrid",
            "predictor_hybrid",
            "vector"
          ),
          ####################################################
          function(data, coeff, pls_score) {
            #               n * 1           1 * M
            for (i in 1 : (data@n_predictor_functional) ){
              data@predictor_functional_list[[i]]@coef <- pls_score %*% coeff@predictor_functional_list[[i]]@coef
            }
            data@Z <-  pls_score %*% coeff@Z
            return(data)
          }

          ####################################################
)




setMethod("hybrid_from_coef", signature("predictor_hybrid"),
          #######################################

          function(format, xi_star){
            M <- dim(format@predictor_functional_list[[1]]@coef)[2] #number of basis
            format@predictor_functional_list[[1]]@coef <- t(as.matrix(xi_star[1:M]))
            format@predictor_functional_list[[2]]@coef <- t(as.matrix(xi_star[(M + 1):(2 * M)]))
            format@Z <- t(as.matrix(xi_star[(2 * M + 1):length(xi_star)]))
            return(format)
          }

          #######################################
)


# get_J_star
#
#


setMethod("get_J_star", signature("predictor_hybrid"),
          #######################################

          function(W){
            p <- dim(W@Z)[2]
            J_star <- W@predictor_functional_list[[1]]@J
            for (i in 2:(W@n_predictor_functional)){
              J_star <- Matrix::bdiag(J_star, W@predictor_functional_list[[i]]@J)
            }
            J_star <- Matrix::bdiag(J_star, diag(p))
            return(J_star)
          }

          #######################################
)



setMethod("get_J_dotdot_star", signature("predictor_hybrid"),
          #######################################

          function(W){
            p <- dim(W@Z)[2]
            J_dotdot_star <- W@predictor_functional_list[[1]]@J_dotdot
            for (i in 2:(W@n_predictor_functional)){
              J_dotdot_star <- Matrix::bdiag(J_dotdot_star, W@predictor_functional_list[[i]]@J_dotdot)
            }
            J_dotdot_star <- Matrix::bdiag(J_dotdot_star, matrix(0, nrow = p, ncol = p))
            return(J_dotdot_star)
          }

          #######################################
)




setMethod("get_Lambda", signature("predictor_hybrid", "numeric"),
          #################################################################

          function(W, lambda){
            p <- dim(W@Z)[2]
            n_basis_1 <- dim(W@predictor_functional_list[[1]]@J)[2]
            Lambda <- lambda[1] * diag(n_basis_1)
            for (i in 2:(W@n_predictor_functional)){
              n_basis_now <- dim(W@predictor_functional_list[[i]]@J)[2]
              Lambda <- Matrix::bdiag(Lambda, lambda[i] * diag(n_basis_now))
            }
            Lambda <- Matrix::bdiag(
              Lambda,
              matrix(0, nrow = p, ncol = p)
            )
          }

          #######################################

)






# get_V_star
#
# computes V*,
# which is the core part of the matrix representation of the covariance,
# which is the objective function of the optimization problem
#
# input:
# W: hybrid predictor object
# y: response vector
#
# output:
# matrix V* (from Proposition 1)
# which is a (MK+p) x (MK+p) square matrix,
# where
# M = number of basis
# K = dimension of the functional vector. In our case, K = 3(R, G, B)
# p = dimension of covariate.
#


setMethod("get_V_star", signature("predictor_hybrid", "numeric"),
          #######################################

          function(W, y){
            n <- W@n_sample
            K <- W@n_predictor_functional
            # building blocks
            ## 1. C multiplied by J. Scalable to arbitrary K
            predictor_first <- W@predictor_functional_list[[1]]
            C_J <- predictor_first@coef %*% predictor_first@J
            for (i in 2:K){
              predictor_now <- W@predictor_functional_list[[i]]
              C_J <- cbind(
                C_J,
                predictor_now@coef %*% predictor_now@J
              )
            }

            ## 2. easy ones
            J_Ct_y <- t(C_J) %*% y
            Zt_y <- t(W@Z) %*% y

            # V matrix
            upper_left <- J_Ct_y %*% t(J_Ct_y)
            upper_right <- J_Ct_y %*% t(Zt_y)
            lower_left <- t(upper_right)
            lower_right <- Zt_y %*% t(Zt_y)

            V_star <- cbind(
              rbind(upper_left, lower_left),
              rbind(upper_right, lower_right)
            )/(n^2)

            return(V_star)
          }

          #######################################
)



# pls_pen
# computes hybrid pls component (loading)
# by maximizing the covariance with regularization
# this function is called at each l-th iteration of PLS fitting

#input:
# W: hyprid predictor object
# y: response vector
# L: regularization matrix

setMethod("pls_pen", "predictor_hybrid",
          #######################################

          function(W, y, L){
            V_star <- get_V_star(W, y)
            #Last paragraph of section 3.2. To solve (6) in practice...
            A <- t(solve(L, t(V_star))) # A = V* t(inv(L))
            E <- solve(L, A) # E = inv(L) A = inv(L) V* t(inv(L))
            eigen_result <- eigen(E)
            cat("eigenvalues: ", eigen_result$values)
            e <- eigen_result$vectors[, 1]
            #e <- e / Matrix::norm(e, "2") #normalize; omit this by request of Dr. Jang

            xi_star <- solve(t(L), e) # t(L) xi* = e
            xi_hat <- hybrid_from_coef(format = W, xi_star = xi_star)


            return(
              list(
                xi = xi_hat,
                E = E,
                V_star = V_star,
                eigen_val = eigen_result$values
              )
            )

          }

          #######################################
)




library(Matrix)

setClass(
  "hybrid_pls_kidney",
  representation(
    eta = "predictor_hybrid",
    xi = "list",
    nu = "list",
    rho = "list",
    delta = "list",
    E = "list",
    V_star = "list",
    eigen_val = "list",
    first_eigen_val = "vector",
    J_Lambda_Jpp = "dgCMatrix",
    L_mat = "dtCMatrix",
    resid_y = "list",
    mse_W = 'vector',
    mse_y = 'vector',
    fitted_value_W = "list",
    fitted_value_y = "list",
    W_now = "list"
  ))

setMethod("predict_test", "hybrid_pls_kidney",
          function(pls_object, W_test) {
            #centering
            return(
              y_pred = hybrid_inner_prod(W_test, pls_object@eta)
            )
          })
