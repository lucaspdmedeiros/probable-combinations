# Reconstruct species pool from observed local communities and computes the 
# probability of persistence of the observed communities and of potential 
# communities sampled from the species pool

# cleaning wd, loading functions and packages ------------------------------
rm(list = ls(all = TRUE))
if(!require(plyr)) {install.packages("plyr"); library(plyr)}
if(!require(mvtnorm)) {install.packages("mvtnorm"); library(mvtnorm)}
if(!require(matrixcalc)) {install.packages("matrixcalc"); library(matrixcalc)}
source("code/functions/build_metanet.R")
source("code/functions/remove_nodes_no_links.R")
source("code/functions/omega.R")
source("code/functions/sub_mats_dim.R")
source("code/functions/scale_matrix.R")
source("code/functions/scaling_factor.R")

# settings ------------------------------
# which dataset to use
dataset <- "year_stage"
# interaction inferrence method
inferrence <- "scaling_factor"

# read matrix files ------------------------------
# file paths
path <- "data/matrices"
all_files <- dir(path)
all_files_sub <- all_files[grep(".txt", all_files)]
file_paths <- paste(path, all_files_sub, sep = "/")
# read files
mats <- lapply(file_paths, read.table, header = TRUE, row.names = 1)
mats <- lapply(mats, as.matrix)
# change matrices to binary
mats <- lapply(mats, function(mat) {mat[mat > 0] <- 1; return(mat)})
# remove species with no interactions
mats <- lapply(mats, remove_nodes_no_links)

# build data frame to store results ------------------------------
net_names <- gsub("\\..*", "", all_files_sub)
names_list <- strsplit(net_names, "_")
plot <- ldply(names_list, function(x) x[1])$V1
new_plot <- plot
new_plot[plot == "1"] <- 1
new_plot[plot == "8"] <- 2
new_plot[plot == "10"] <- 3
new_plot[plot == "3"] <- 4
new_plot[plot == "7"] <- 5
new_plot[plot == "9"] <- 6
new_plot[plot == "4"] <- 7
new_plot[plot == "5"] <- 8
new_plot[plot == "12"] <- 9
plot <- as.factor(new_plot)
stage <- factor(ldply(names_list, function(x) x[2])$V1, 
                levels = c("initial", "middle", "late"))
year <- factor(ldply(names_list, function(x) x[3])$V1,
               levels = c("2007", "2008", "2009", "2010", 
                          "2011", "2012", "2013", "2014", "2016", "2017"))
results_df <- data.frame(stage = stage,
                         plot = plot,
                         year = year,
                         n_sp_real = rep(NA, length(all_files_sub)),
                         n_sp_pot = rep(NA, length(all_files_sub)),
                         Omega_pool = rep(NA, length(all_files_sub)),
                         omega_pool = rep(NA, length(all_files_sub)),
                         Omega_real = rep(NA, length(all_files_sub)),
                         omega_real = rep(NA, length(all_files_sub)),
                         Omega_percentile = rep(NA, length(all_files_sub)),
                         omega_percentile = rep(NA, length(all_files_sub)),
                         Omega_mean = rep(NA, length(all_files_sub)),
                         omega_mean = rep(NA, length(all_files_sub)),
                         comp_strength_real = rep(NA, length(all_files_sub)),
                         comp_strength_mean = rep(NA, length(all_files_sub)))
# define years and stages to use
years <- c("2007", "2008", "2009", "2010", "2011",
           "2012", "2013", "2014", "2016", "2017")
stages <- c("initial", "middle", "late")
# define groups of matrices to build species pool
if (dataset == "year_stage") {
  groups <- apply(expand.grid(stages, years), 1, paste, collapse = "_")
}
if (dataset == "year") {
  groups <- years
}
if (dataset == "stage") {
  groups <- stages
}

# loop over all species pools ------------------------------
for (i in 1:length(groups)) {
  # build species pool
  group <- groups[i]
  files_id <- grep(group, all_files_sub)
  metanet <- build_metanet(mats[files_id])
  
  # loop over all matrices within the species pool ------------------------------
  for (j in 1:length(files_id)) {
    # build matrix
    mat_id <- files_id[j]
    # print name of current matrix
    print(gsub("\\..*", "", all_files_sub[mat_id]))
    # extract plant labels from matrix
    plant_labels <- rownames(mats[[mat_id]])
    # create potential matrix
    potential_net <- metanet[plant_labels, ]
    # remove species with no interactions
    potential_net <- remove_nodes_no_links(potential_net)
    # create monopartite projection (herbivore competition matrix)
    proj_pot <- t(potential_net) %*% potential_net
    # compute number of species of potential matrix
    n_sp_pot <- nrow(proj_pot)
    # rescale matrix
    if (inferrence == "col_sum") {
      comp_net_pot <- scale_matrix(proj_pot, remove_zero = TRUE, 
                                   scale = "col_sum", p = 1,
                                   tol = 10^-8, diag = "one")
    }
    if (inferrence == "scaling_factor") {
      comp_net_pot <- scale_matrix(proj_pot, remove_zero = TRUE, 
                                   scale = "scaling_factor", p = 0.5,
                                   tol = 10^-8, diag = "preserve")
    }
    # test if matrix is positive definite
    if(!is.positive.definite(comp_net_pot + t(comp_net_pot)))
      stop("Species pool matrix is not positive definite")
    # create realized matrix
    real_net <- as.matrix(mats[[mat_id]])
    # create monopartite projection
    proj_real <- t(real_net) %*% real_net
    # rescale matrix
    if (inferrence == "col_sum") {
      comp_net_real <- scale_matrix(proj_real, remove_zero = TRUE, 
                                    scale = "col_sum", p = 1,
                                    tol = 10^-8, diag = "one")
    }
    if (inferrence == "scaling_factor") {
      comp_net_real <- scale_matrix(proj_real, remove_zero = TRUE, 
                                    scale = "scaling_factor", p = 0.5,
                                    tol = 10^-8, diag = "preserve")
    }
    # test if matrix is positive definite
    if(!is.positive.definite(comp_net_real + t(comp_net_real)))
      stop("Observed matrix is not positive definite")
    # compute number of species of realized matrix
    n_sp_real <- nrow(comp_net_real)
    # sample sub-matrices from potential matrix
    n_samples <- 1000
    proj_same_dim <- sub_mats_dim(mat = proj_pot, dim = n_sp_real, 
                                  sampling_dim = "both", all = FALSE, n_sample = n_samples)
    # rescale matrix
    if (inferrence == "col_sum") {
      comp_net_same_dim <- lapply(proj_same_dim, scale_matrix, remove_zero = TRUE,
                                  scale = "col_sum", p = 1, tol = 10^-8, diag = "one")
    }
    if (inferrence == "scaling_factor") {
      comp_net_same_dim <- lapply(proj_same_dim, scale_matrix, remove_zero = TRUE,
                                  scale = "scaling_factor", p = 0.5, tol = 10^-8, diag = "preserve")
    }
    # remove NA values
    comp_net_same_dim <- comp_net_same_dim[!is.na(comp_net_same_dim)]
    # test if matrices are positive definite
    if(!all(sapply(lapply(comp_net_same_dim, function(mat) mat + t(mat)),
                   is.positive.definite)))
      stop("At least one potential matrix is not positive definite")
    # transform values to negative
    comp_net_pot_neg <- - comp_net_pot
    comp_net_real_neg <- - comp_net_real
    comp_net_same_dim <- lapply(comp_net_same_dim, function(mat) - mat)
    # compute size of feasibility domains
    omega_pot_list <- Omega(comp_net_pot_neg)
    Omega_pot <- (2^n_sp_pot) * omega_pot_list[[1]]
    omega_pot <- 2 * omega_pot_list[[2]]
    omega_real_list <- Omega(comp_net_real_neg)
    Omega_real <- (2^n_sp_real) * omega_real_list[[1]]
    omega_real <- 2 * omega_real_list[[2]]
    omega_same_dim_list <- sapply(comp_net_same_dim, Omega)
    Omega_same_dim <- (2^n_sp_real) * unlist(omega_same_dim_list[1, ])
    omega_same_dim <- 2 * unlist(omega_same_dim_list[2, ])
    # compute omega percentile
    Omega_ecdf <- ecdf(Omega_same_dim)
    Omega_percentile <- Omega_ecdf(Omega_real)
    omega_ecdf <- ecdf(omega_same_dim)
    omega_percentile <- omega_ecdf(omega_real)
    # compute mean omega
    Omega_mean <- mean(Omega_same_dim, na.rm = TRUE)
    omega_mean <- mean(omega_same_dim, na.rm = TRUE)
    # compute mean competition strength
    comp_strength_real <- mean(comp_net_real[row(comp_net_real) != col(comp_net_real)])
    comp_strength_same_dim <- - sapply(comp_net_same_dim, function(mat) mean(mat[row(mat) != col(mat)]))
    comp_strength_same_dim_mean <- mean(comp_strength_same_dim, na.rm = TRUE)
    
    # build results data frame ------------------------------
    # add results to data frame
    results_df$n_sp_real[mat_id] <- n_sp_real
    results_df$n_sp_pot[mat_id] <- n_sp_pot
    results_df$Omega_pool[mat_id] <- Omega_pot
    results_df$omega_pool[mat_id] <- omega_pot
    results_df$Omega_percentile[mat_id] <- Omega_percentile
    results_df$omega_percentile[mat_id] <- omega_percentile
    results_df$Omega_real[mat_id] <- Omega_real
    results_df$omega_real[mat_id] <- omega_real
    results_df$Omega_mean[mat_id] <- Omega_mean
    results_df$omega_mean[mat_id] <- omega_mean
    results_df$comp_strength_real[mat_id] <- comp_strength_real
    results_df$comp_strength_mean[mat_id] <- comp_strength_same_dim_mean
  }
}

# save result file ------------------------------
write.csv(results_df, 
          paste("results/fig3/results_2007-2017_", dataset, "_metanetwork_", inferrence, ".csv", sep = ""), 
          row.names = FALSE)
