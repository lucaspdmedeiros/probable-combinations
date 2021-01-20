# Function that scales an interaction matrix using different methods.
# mat: a matrix
# remove_zero: logical; whether to remove rows and columns filled with zeroes
# scale: how to scale the off-diagonal elements; "none" does not scale elements,
# "col_sum" normalizes columns to sum one and "scaling_factor" scales off-diagonal 
# elements using a scaling factor
# p: the proportion of the scaling factor to use; the scaling factor is computed
# as the maximum the matrix tolarates before one of its eigenvalues turns zero
# tol: the tolerance to use to define a zero eigenvalue
# diag: how to set diagonal elements after scaling; "preserve" preserves
# the original values, "one" sets values to one, and "scale" keeps the
# values scaled

scale_matrix <- function(mat, remove_zero = TRUE, scale = "col_sum", 
                         p = 1, tol = 0.00000001, diag = "one") {
  # remove rows and columns filled with zeroes
  if (remove_zero == TRUE) {
    mat <- remove_nodes_no_links(mat)
  }
  # obtain diagonal
  diag_mat <- diag(mat)
  # scale matrix
  if (scale == "none") {
    scaled_mat <- mat
  }
  if (scale == "col_sum") {
    scaled_mat <- scale(mat, center = FALSE, 
                        scale = colSums(mat))
    attr(scaled_mat, "scaled:center") <- NULL 
    attr(scaled_mat, "scaled:scale") <- NULL 
  }
  if (scale == "scaling_factor") {
    scal_fac <- scaling_factor(mat = mat, tol = tol)
    if (is.na(scal_fac))
      return(NA)
    scaled_mat <- mat * scal_fac * p
  }
  # change diagonal
  if (diag == "preserve") {
    diag(scaled_mat) <- diag_mat
    return(scaled_mat)
  }
  if (diag == "one") {
    diag(scaled_mat) <- 1
    return(scaled_mat)
  }
  if (diag == "scale") {
    return(scaled_mat)
  }
}