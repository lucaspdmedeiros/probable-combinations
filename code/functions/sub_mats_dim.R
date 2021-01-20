# Obtains submatrices of a given dimension from a matrix and returns
# the submatrices in a list.
# mat: a bipartite or monopartite matrix
# dim: the dimension (number of rows and/or columns) of the sampled matrices
# sampling_dim: where to sample nodes from; "row" means that only 
# rows will be sampled, "col" means that only columns will be sampled,
# and "both" means that both rows and columns will be sampled 
# (i.e., a principal submatrix is sampled)
# all: logical; should all possible matrices be returned or just a sample?
# n_sample: if all == FALSE, the number of samples that should be taken

sub_mats_dim <- function(mat, dim, sampling_dim = "both", 
                         all = FALSE, n_sample = 100) {
  n_row <- nrow(mat)
  n_col <- ncol(mat)
  if (all == TRUE) {
    if (sampling_dim == "row") {
      ind_list_row <- combn(seq(1, n_row, by = 1), dim, simplify = FALSE)
      sub_mat_list <- lapply(ind_list_row, function (ind, A) A[ind, ], A = mat)
    }
    if (sampling_dim == "col") {
      ind_list_col <- combn(seq(1, n_col, by = 1), dim, simplify = FALSE)
      sub_mat_list <- lapply(ind_list_col, function (ind, A) A[ , ind], A = mat)
    }
    if (sampling_dim == "both") {
      n <- min(c(n_row, n_col))
      ind_list <- combn(seq(1, n, by = 1), dim, simplify = FALSE)
      sub_mat_list <- lapply(ind_list, function (ind, A) A[ind, ind], A = mat)
    }
  } else {
    if (sampling_dim == "row") {
      ind_list_row <- replicate(n_sample, 
                                sample(x = 1:n_row, size = dim, replace = FALSE), 
                                simplify = FALSE)
      sub_mat_list <- lapply(ind_list_row, function (ind, A) A[ind, ], A = mat)
    }
    if (sampling_dim == "col") {
      ind_list_col <- replicate(n_sample, 
                                sample(x = 1:n_col, size = dim, replace = FALSE), 
                                simplify = FALSE)
      sub_mat_list <- lapply(ind_list_col, function (ind, A) A[ , ind], A = mat)
    }
    if (sampling_dim == "both") {
      n <- min(c(n_row, n_col))
      ind_list <- replicate(n_sample, 
                            sample(x = 1:n, size = dim, replace = FALSE), 
                            simplify = FALSE)
      sub_mat_list <- lapply(ind_list, function (ind, A) A[ind, ind], A = mat)
    }
  }
  return(sub_mat_list)
}