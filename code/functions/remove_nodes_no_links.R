# Removes nodes with no links (degree = 0) from a network.
#
# Args:
#   mat: bipartite adjacency matrix
#   warn: logical; whether a warning should be printed out with the rows/columns
#         that were removed
#           
# Returns:
#   A new matrix without the row(s) and/or column(s) corresponding to nodes with no links.
#   A warning is printed indicating which nodes have been removed.

remove_nodes_no_links = function(mat, warn = FALSE) {
  mat <- as.matrix(mat)
  # degree of row nodes
  k_row = apply(mat, 1, sum)
  # degree of columns nodes
  k_col = apply(mat, 2, sum)
  if (any(k_row == 0)) { 
    # node positions
    k_row_zeros = which(k_row == 0)
    # remove nodes
    mat = mat[-k_row_zeros, ] 
    if (warn == TRUE) {
      warning(paste("Row(s)", k_row_zeros, "have been removed from the matrix. "))
    }
  }
  if (any(k_col == 0)) { 
    # node positions
    k_col_zeros = which(k_col == 0)
    # remove nodes
    mat = mat[ ,-k_col_zeros]
    if (warn == TRUE) {
      warning(paste("Column(s)", k_col_zeros, "have been removed from the matrix. "))
    }
  }
  return(mat)
}