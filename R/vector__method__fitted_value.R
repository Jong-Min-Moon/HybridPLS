setMethod("fitted_value",
          signature(
            "vector","numeric","vector"
          ),
#############################################
          function(data, coeff, pls_score) {
            #  n * 1    n* 1      scalar
            return(pls_score * coeff)
          }

#############################################
)
