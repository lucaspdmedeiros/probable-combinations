# Function that creates a metanetwork using a list containing 
# different networks. The matrices in the list have to be 
# labeled so that the species are correctly pooled together.

build_metanet <- function(mats) {
  # obtain row labels for all matrices
  mats_row_names <- unlist(lapply(mats, rownames))
  unique_mats_row_names <- unique(mats_row_names)
  unique_mats_row_names <- sort(unique_mats_row_names)
  # obtain column labels for all matrices
  mats_col_names <- unlist(lapply(mats, colnames))
  unique_mats_col_names <- unique(mats_col_names)
  unique_mats_col_names <- sort(unique_mats_col_names)
  # create matrix with zeros
  mat_all <- matrix(0, nrow = length(unique(mats_row_names)),
                    ncol = length(unique(mats_col_names)))
  # label the matrix using all the obtained labels
  rownames(mat_all) <- unique_mats_row_names
  colnames(mat_all) <- unique_mats_col_names
  # fill up the metanetwork using the individual networks
  for (i in 1:nrow(mat_all)) {
    # extract species name
    sp_name <- rownames(mat_all)[i]
    for (j in 1:length(mats)) {
      # if this species is present in matrix j
      if (any(sp_name == rownames(mats[[j]]))) {
        # match the names of column species
        match_names <- match(names(mats[[j]][sp_name, ]), names(mat_all[sp_name, ]))
        # add the interactions to the metanetwork
        mat_all[sp_name, match_names] <- as.numeric(mat_all[sp_name, match_names] + mats[[j]][sp_name, ])
      }
    }
  }
  # making network binary
  mat_all[mat_all > 0] <- 1
  return(mat_all)
}
