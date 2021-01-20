# Function that samples a niche overlap competition matrix for a given average niche overlap value
# S: number of species 
# rho: average niche overlap
# same_width: whether to use the same niche width (assumed to be 1) for all species

sample_niche_overlap_matrix <- function(S, rho, same_width) {
  # niche centers
  mu <- runif(n = S, min = 0, max = 1)
  # squared difference between niche centers
  #A <- outer(mu, mu, FUN = "-")^2
  A <- as.matrix(dist(mu))
  if (same_width) {
    # function that computes the squared difference between mean of
    # off-diagonal coefficients and desired average niche overlap
    f <- function(sigma2) {
      # interaction matrix 
      alpha <- exp(- A / (4 * sigma2))
      # squared difference
      out <- ((sum(alpha) - S) / ((S - 1) * S) - rho)^2
      return(out)
    }
    # find sigma squared that minimizes function above
    sigma2_opt <- optimize(f, interval = c(0, 100), maximum = FALSE, tol = 10^-14)$minimum
    # final interaction matrix
    alpha <- exp(- A / (4 * sigma2_opt))
  } else {
    # niche widths (for rho approximately 0.05)
    sigma2 <- runif(n = S, min = 0.005, max = 0.009)
    # final interaction matrix
    B <- outer(sigma2, sigma2, FUN = "+")
    alpha <- exp(- A / (2 * B))
  }
  return(alpha)
}
