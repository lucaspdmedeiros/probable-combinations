# Function that computes the size of the feasibility domain of an interaction 
# matrix and outputs a list with the unscaled and scaled (by number of species) 
# size of the feasibility domain

Omega <- function(alpha) {
  # number of species
  S <- nrow(alpha)
  # implements quasi Monte Carlo integration of the multivariate normal distribution
  omega <- function(S, Sigma) {
    out <- pmvnorm(lower = rep(0, S), upper = rep(Inf, S), mean = rep(0, S), sigma = Sigma)[1]
    return(list(out, out^(1 / S)))
  }
  # compute matrix to invert
  B <- t(alpha) %*% alpha
  # if symmetric and positive definite, invert using Choleski decomposition
  if (isSymmetric(B)) {
    if (is.positive.definite(B)) {
      chol_dec <- chol(B)
      Sigma <- chol2inv(chol_dec)
      return(omega(S, Sigma))
    }
  }
  # attempts to invert B directly and returns 0 if fails
  if (!(class(try(solve(B), silent = TRUE)) == "matrix"))
    return(0)
  # compute inverse directly
  else {
    Sigma <- solve(B)
    return(omega(S, Sigma))
  }
}
