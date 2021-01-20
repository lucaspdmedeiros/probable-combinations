# Function that builds a bipartite erdos-renyi matrix 
# with a given size and connectance

bip_erdos_renyi <- function(rows, columns, connect) {
  mat <- matrix(data = 1, nrow = rows, ncol = columns)
  prob_mat <- matrix(data = runif(prod(dim(mat)), 0, 1), 
                     nrow = rows, ncol = columns)
  mat[prob_mat > connect] <- 0
  while(any(apply(mat, 1, sum) == 0) | any(apply(mat, 2, sum) == 0)) {
    mat <- matrix(data = 1, nrow = rows, ncol = columns)
    prob_mat <- matrix(data = runif(prod(dim(mat)), 0, 1), 
                       nrow = rows, ncol = columns)
    mat[prob_mat > connect] <- 0
  }
  return(mat)
}