# Generates multiple species pools according to a given framework 
# (e.g., niche framework) and computes the probability of persistence 
# of each species pool

# cleaning wd, loading functions and packages ------------------------------
rm(list = ls(all = TRUE))
set.seed(42)
if(!require(plyr)) {install.packages("plyr"); library(plyr)}
if(!require(deSolve)) {install.packages("deSolve"); library(deSolve)}
if(!require(mvtnorm)) {install.packages("mvtnorm"); library(mvtnorm)}
if(!require(matrixcalc)) {install.packages("matrixcalc"); library(matrixcalc)}
source("code/functions/lotka_volterra.R")
source("code/functions/lv_pruning.R")
source("code/functions/omega.R")
source("code/functions/sub_mats_dim.R")
source("code/functions/scale_matrix.R")
source("code/functions/scaling_factor.R")
source("code/functions/sample_niche_overlap_matrix.R")
source("code/functions/hypersphere_sampling.R")
source("code/functions/bip_erdos_renyi.R")
source("code/functions/remove_nodes_no_links.R")

# settings ------------------------------
# species pool size
n_sp_pool <- 20
# number of species pools
n_pool <- 500
# how to build species pool matrices
build_pool <- "niche"
# number of rows (bipartite framework)
rows <- 20
# number of columns (bipartite framework)
columns <- 20
# connectance (bipartite framework)
connect <- c(0.1, 0.2, 0.3)
# average niche overlap (niche framework)
rho <- c(0.025, 0.05, 0.1)
# whether to normalize matrices
normalize <- FALSE
# whether to use same niche width for all species
same_width <- TRUE
# competition strength
u_mean <- 1
# mean of coefficients
coef_mean <- 0
# sd of coefficients
coef_sd <- 1
# defining parameter to loop over
if (build_pool == "niche")
  param <- rho
if (build_pool == "bipartite")
  param <- connect
if (build_pool == "equal")
  param <- u_mean
if (build_pool == "random")
  param <- coef_sd

# sample metaneworks ------------------------------
# build data frame to store results
results_df <- data.frame(comm = rep(1:n_pool, length(param)),
                         param = as.factor(rep(param, each = n_pool)),
                         Omega_pool = rep(NA, length(param) * n_pool),
                         omega_pool = rep(NA, length(param) * n_pool),
                         comp_strength_pool = rep(NA, length(param) * n_pool))
# loop over parameters and species pools
for (i in 1:length(param)) {
  print(paste("param", param[i]))
  for (j in 1:n_pool) {
    print(paste("community", j))
    # build matrix of species pool
    if (build_pool == "niche")
      comp_net_pool <- sample_niche_overlap_matrix(n_sp_pool, param[i], same_width)
    if (build_pool == "bipartite") {
      comp_net_pool <- bip_erdos_renyi(rows, columns, param[i])
      comp_net_pool <- t(comp_net_pool) %*% comp_net_pool
      while (any(apply(comp_net_pool, 1, sum) == 0) | any(apply(comp_net_pool, 2, sum) == 0)) {
        comp_net_pool <- bip_erdos_renyi(rows, columns, param[i])
        comp_net_pool <- t(comp_net_pool) %*% comp_net_pool
      }
    }
    if (build_pool == "equal") {
      u_mean <- param[i]
      u <- u_mean / n_sp_pool
      alpha <- u_mean
      comp_net_pool <- matrix(u, nrow = n_sp_pool, ncol = n_sp_pool)
      diag(comp_net_pool) <- alpha
    }
    if (build_pool == "random") {
      comp_net_pool <- matrix(rnorm(n_sp_pool * n_sp_pool, coef_mean, param[i]), 
                              nrow = n_sp_pool, ncol = n_sp_pool)
      leading_eigen <- max(Re(eigen(comp_net_pool + t(comp_net_pool))$values))
      d <- - leading_eigen - 0.0001
      comp_net_pool <- comp_net_pool - diag(d, nrow = nrow(comp_net_pool))
    }
    # add species labels
    rownames(comp_net_pool) <- paste("sp", 1:nrow(comp_net_pool), sep = "")
    colnames(comp_net_pool) <- paste("sp", 1:nrow(comp_net_pool), sep = "")
    # normalize matrix
    if (normalize)
      comp_net_pool <- scale_matrix(comp_net_pool, remove_zero = TRUE, 
                                    scale = "col_sum", p = 1,
                                    tol = 10^-8, diag = "one")
    # test if matrix is positive definite
    if(!is.positive.definite(comp_net_pool + t(comp_net_pool)))
      stop("Initial matrix is not positive definite")
    # setting coefficients to negative
    comp_net_pool_neg <- - comp_net_pool
    # compute omega of species pool matrix
    omega_pool_list <- Omega(comp_net_pool_neg)
    Omega_pool <- (2^n_sp_pool) * omega_pool_list[[1]]
    omega_pool <- 2 * omega_pool_list[[2]]
    if (build_pool == "random") {
      Omega_pool <- omega_pool_list[[1]]
      omega_pool <- omega_pool_list[[2]]
    }
    # compute competition strength of species pool matrix
    comp_strength_pool <- mean(comp_net_pool[row(comp_net_pool) != col(comp_net_pool)])
    # add results to data frame
    results_df$Omega_pool[results_df$param == param[i] & results_df$comm == j] <- Omega_pool
    results_df$omega_pool[results_df$param == param[i] & results_df$comm == j] <- omega_pool
    results_df$comp_strength_pool[results_df$param == param[i] & results_df$comm == j] <- comp_strength_pool
  }
}
# saving results
write.csv(results_df, paste("results/fig2/", build_pool, "/results_species_pool_theoretical_matrices_n", 
                            n_sp_pool, "_", build_pool, "pool", "_ncomm", n_pool,
                            ".csv", sep = ""), row.names = FALSE)
