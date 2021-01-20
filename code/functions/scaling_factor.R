# Function that computes the largest scaling factor that multiplies the 
# off-diagonal elements of a matrix and the matrix plus its transpose has 
# all eigenvalues positive (Saavedra et al 2017 Jour Anim Ecol)
# mat: a matrix
# tol: an eigenvalue below this value is considered zero

scaling_factor <- function(mat, tol = 0.00000001) {
  mu <- 0
  for (i in 1:10000) {
    mu <- mu + 0.001
    curr_mat <- mat * mu
    diag(curr_mat) <- diag(mat)
    min_eigen <- tryCatch(min(Re(eigen(curr_mat + t(curr_mat))$values)), error = function(e) NA)
    if (is.na(min_eigen))
      return(NA)
    if (min_eigen < tol)
      return(mu - 0.001)
  }
  return(mu)
}
