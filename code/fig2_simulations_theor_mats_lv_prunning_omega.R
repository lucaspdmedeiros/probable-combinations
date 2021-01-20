# Run Lotka-Volterra simulations for species pools built according to different 
# theoretical frameworks and computes the size of the feasibility domain of 
# realized combinations of species and of potential combinations of species

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

# simulation settings ------------------------------
# species pool size
n_sp_pool <- 20
# number of species pools
n_pool <- 10
# how to build species pool matrices
build_pool <- "niche"
# which LV model to use
formalism <- "K"
# number of rows (bipartite framework)
rows <- 20
# number of columns (bipartite framework)
columns <- 20
# connectance (bipartite framework)
connect <- c(0.1, 0.2, 0.3)
# average niche overlap (niche framework)
rho <- c(0.025, 0.05, 0.1)
# whether to use same niche width for all species
same_width <- TRUE
# competition strength
u_mean <- 1
# mean of coefficients
coef_mean <- 0
# sd of coefficients
coef_sd <- 1
# whether to normalize matrices
normalize <- FALSE
# number of LV simulations
n_sim <- 100
# number of potential pruned communities
n_samples <- 500
# defining parameter to loop over
if (build_pool == "niche")
  param <- rho
if (build_pool == "bipartite")
  param <- connect
if (build_pool == "equal")
  param <- u_mean
if (build_pool == "random")
  param <- coef_sd

# loop over parameters and species pools
for (i in 1:length(param)) {
  # sample metaneworks ------------------------------
  print(paste("parameter", param[i]))
  # data frame to store results
  all_results_df <- data.frame()
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
    mean_comp_strength_pool <- mean(comp_net_pool[row(comp_net_pool) != col(comp_net_pool)])
    
    # run simulations of LV dynamics ------------------------------
    # LV simulations settings
    times <- seq(0, 100, 0.01)
    N <- rep(1, n_sp_pool)
    K <- replicate(n_sim, hypersphere_sampling(n_sp_pool, positive = TRUE, within = FALSE), simplify = FALSE)
    if (build_pool == "random") {
      K <- replicate(n_sim, hypersphere_sampling(n_sp_pool, positive = FALSE, within = FALSE), simplify = FALSE)
    }
    r <- K
    # simulate and obtain indexes of surviving species
    which_surv_sp <- mapply(lv_pruning, r = r, K = K,
                            MoreArgs = list(N = N, A = comp_net_pool_neg, times = times, 
                                            extinct_tol = 10^-2, formalism = formalism),
                            SIMPLIFY = FALSE)
    # obtain pruned competition matrices
    comp_net_pruned <- lapply(which_surv_sp, function(mat, sp) mat[sp, sp],
                              mat = comp_net_pool)
    # compute size of pruned matrices
    n_sp_pruned <- sapply(comp_net_pruned, nrow)
    # remove communities with zero or one species
    n_sp_pruned[sapply(n_sp_pruned, is.null)] <- NA
    comp_net_pruned <- comp_net_pruned[!is.na(n_sp_pruned)]
    n_sp_pruned <- unlist(n_sp_pruned[!is.na(n_sp_pruned)])
    # normalize matrices
    if (normalize)
      comp_net_pruned <- lapply(comp_net_pruned, scale_matrix, remove_zero = TRUE,
                                scale = "col_sum", p = 1, tol = 10^-8, diag = "one")
    # test if matrices are positive definite
    if(!all(sapply(lapply(comp_net_pruned, function(mat) mat + t(mat)),
                   is.positive.definite)))
      stop("At least one realized matrix is not positive definite")
    # setting coefficients to negative
    comp_net_pruned <- lapply(comp_net_pruned, function(mat) - mat)
    # compute size of feasibility domains of pruned matrices
    omega_pruned_list <- sapply(comp_net_pruned, Omega)
    Omega_pruned <- (2^n_sp_pruned) * unlist(omega_pruned_list[1, ])
    omega_pruned <- 2 * unlist(omega_pruned_list[2, ])
    if (build_pool == "random") {
      Omega_pruned <- unlist(omega_pruned_list[1, ])
      omega_pruned <- unlist(omega_pruned_list[2, ])
    }
    # compute mean competition strength
    mean_comp_strength <- - sapply(comp_net_pruned, function(mat) mean(mat[row(mat) != col(mat)]))
    
    # build data frame to store results
    results_df <- data.frame(comm = rep(j, length(n_sp_pruned)),
                             param = rep(param[i], length(n_sp_pruned)),
                             n_sp_pruned = n_sp_pruned, 
                             n_sp_pool = n_sp_pool,
                             mean_comp_strength_pool = mean_comp_strength_pool,
                             mean_comp_strength = mean_comp_strength,
                             Omega_pool = rep(Omega_pool, length(n_sp_pruned)),
                             omega_pool = rep(omega_pool, length(n_sp_pruned)),
                             Omega_obs = Omega_pruned, 
                             omega_obs = omega_pruned, 
                             Omega_percentile = rep(NA, length(n_sp_pruned)),
                             omega_percentile = rep(NA, length(n_sp_pruned)),
                             Omega_expect_comm = rep(NA, length(n_sp_pruned)),
                             omega_expect_comm = rep(NA, length(n_sp_pruned)),
                             comp_strength_expect = rep(NA, length(n_sp_pruned)))
    
    # compute size of feasibility domain of sub-matrices ------------------------------
    # loop over sizes
    for (k in 1:length(n_sp_pruned)) {
      print(k)
      # compute number of species of pruned network
      n_sp_real <- n_sp_pruned[k]
      # sample sub-communities from potential network
      comp_net_same_dim <- sub_mats_dim(mat = comp_net_pool, dim = n_sp_real, 
                                        sampling_dim = "both", all = FALSE, n_sample = n_samples)
      # normalize matrices
      if (normalize)
        comp_net_same_dim <- lapply(comp_net_same_dim, scale_matrix, remove_zero = TRUE,
                                    scale = "col_sum", p = 1, tol = 10^-8, diag = "one")
      # test if matrices are positive definite
      if(!all(sapply(lapply(comp_net_same_dim, function(mat) mat + t(mat)), 
                     is.positive.definite)))
        stop("At least one potential matrix is not positive definite")
      # setting coefficients to negative
      comp_net_same_dim <- lapply(comp_net_same_dim, function(mat) - mat)
      # compute size of feasibility domain of sampled networks
      omega_same_dim_list <- sapply(comp_net_same_dim, Omega)
      Omega_same_dim <- (2^n_sp_real) * unlist(omega_same_dim_list[1, ])
      omega_same_dim <- 2 * unlist(omega_same_dim_list[2, ])
      if (build_pool == "random") {
        Omega_same_dim <- unlist(omega_same_dim_list[1, ])
        omega_same_dim <- unlist(omega_same_dim_list[2, ])
      }
      # compute omega percentile of pruned matrix
      Omega_ecdf <- ecdf(Omega_same_dim)
      Omega_percentile <- Omega_ecdf(Omega_pruned[k])
      omega_ecdf <- ecdf(omega_same_dim)
      omega_percentile <- omega_ecdf(omega_pruned[k])
      # compute mean omega
      Omega_expect_comm <- mean(Omega_same_dim, na.rm = TRUE)
      omega_expect_comm <- mean(omega_same_dim, na.rm = TRUE)
      # compute mean competition strength
      mean_comp_strength_same_dim <- -sapply(comp_net_same_dim, function(mat) mean(mat[row(mat) != col(mat)]))
      comp_strength_expect <- mean(mean_comp_strength_same_dim)
      
      # build results data frame ------------------------------
      # add results of percentiles
      results_df$Omega_percentile[k] <- Omega_percentile
      results_df$omega_percentile[k] <- omega_percentile
      # add results of omega
      results_df$Omega_expect_comm[k] <- Omega_expect_comm
      results_df$omega_expect_comm[k] <- omega_expect_comm
      # add results of competition strength
      results_df$comp_strength_expect[k] <- comp_strength_expect
    }
    all_results_df <- rbind(all_results_df, results_df)
  }
  # saving all results
  write.csv(all_results_df, paste("results/fig2/", build_pool, "/results_", formalism, "_lv_pruning_theoretical_matrices_n", 
                                  n_sp_pool, "_", build_pool, "pool", "_ncomm", n_pool, 
                                  "_nsamples", n_samples, "_param", param[i], 
                                  ".csv", sep = ""), row.names = FALSE)
}
