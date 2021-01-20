# Simulates Lotka-Volterra dynamics for multiple environments (K vectors) 
# and a single initial condition and records the identities of surviving species

# cleaning wd, loading functions and packages ------------------------------
rm(list = ls(all = TRUE))
if(!require(deSolve)) {install.packages("deSolve"); library(deSolve)}
source("code/functions/lotka_volterra.R")
source("code/functions/lv_pruning.R")
source("code/functions/hypersphere_sampling.R")
set.seed(42)

# reading interaction matrix file ------------------------------
mat_name <- "saavedra_etal_2017_7B.txt"
A <- as.matrix(read.table(paste("data/fig1/", mat_name, sep = "")))
# extracting number of species
n_sp <- nrow(A)
rownames(A) <- paste("sp", 1:n_sp, sep = "")
colnames(A) <- paste("sp", 1:n_sp, sep = "")

# setting simulation parameters ------------------------------
# initial condition
N <- rep(1, n_sp)
# define r vectors
r <- rep(0.05, n_sp)
# number of simulations
n_sim <- 2000
# time steps
times <- seq(0, 200, 0.01)
# r or k formalism
formalism <- "K"

# performing simulations ------------------------------
df <- data.frame()
# sample K vectors on the hypersphere
K <- replicate(n_sim, hypersphere_sampling(n_sp, positive = TRUE, within = FALSE), simplify = FALSE)
# run simulations and extract surviving species
surv_sp <- mapply(lv_pruning, K = K,
                  MoreArgs = list(A = A, N = N, times = times, r = r,
                                  formalism = formalism, extinct_tol = 0.01), 
                  SIMPLIFY = FALSE) 
# collapse names
surv_sp <- lapply(surv_sp, function(x) paste(x, collapse = "_"))
# build data frame with results
k_df <- data.frame(matrix(unlist(K), nrow = length(K), byrow = TRUE))
names(k_df) <- paste("K", 1:nrow(A), sep = "")
N_df <- as.data.frame(matrix(rep(N, each = nrow(k_df)), ncol = nrow(A)))
names(N_df) <- colnames(A)
df <- cbind(k_df, N_df)
df$surv_sp <- unlist(surv_sp)

# save results ------------------------------
write.csv(df, paste("results/fig1/lv_surv_sp_", substr(mat_name, 1, nchar(mat_name) - 4), ".csv", sep = ""), 
          row.names = FALSE)
