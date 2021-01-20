# Function that solves the Lotka-Volterra dynamics and returns the 
# indexes of surviving species (i.e., abundance larger than threshold)
# N: vector of initial population sizes
# r: vector of intrinsic growth rates
# K: vector of carrying capacities
# A: square interaction matrix
# times: sequence of times points to numerically integrate ODE
# extinct_tol: species with a final population size smaller than
# this value are considered extinct
# formalism: whether to use r or K Lotka-Volterra formalism

lv_pruning <- function(N, r, K, A, times, extinct_tol = 0.00000001, formalism) {
  # solve Lotka-Volterra dynamics
  out <- lotka_volterra(N, r, K, A, times, formalism)
  eq <- out[nrow(out), -1]
  # return pruned system
  which_surv_sp <- as.numeric(which(eq > extinct_tol))
  return(which_surv_sp)
}