# Solves the system of ordinary differential equations
# given by the generalized Lotka-Volterra dynamics and
# returns the state variables over time

lotka_volterra <- function(N, r, K, A, times, formalism) {
  if (formalism == "r") {
    # list of parameters
    pars <- list(r = r, A = A)
    # function that returns rate of change
    model <- function(t, N, pars) {
      dN_dt <- N * (pars$r + c(pars$A %*% N))
      return(list(dN_dt))
    }
    # numerical integration
    out <- ode(y = N, times = times, func = model, parms = pars)
  }
  if (formalism == "K") {
    pars <- list(r = r, K = K, A = A)
    model <- function(t, N, pars) {
      dN_dt <- N * (pars$r / pars$K) * (pars$K + c(pars$A %*% N))
      return(list(dN_dt))
    }
    out <- ode(y = N, times = times, func = model, parms = pars)
  }
  if (formalism == "K_typeII") {
    pars <- list(r = r, K = K, A = A)
    model <- function(t, N, pars) {
      dN_dt <- pars$r * N * (1 + (1 / pars$K) * c(pars$A %*% diag(1 / (1 + N)) %*% N))
      return(list(dN_dt))
    }
    out <- ode(y = N, times = times, func = model, parms = pars)
  }
  if (formalism == "K_stochastic") {
    pars <- list(r = r, K = K, A = A)
    # defining deterministic part
    f <- function(u, p, t) {
      deterministic <- u * (p$r / p$K) * (p$K + c(p$A %*% u))
      return(deterministic)
    }
    # defining stochastic part
    g <- function(u, p, t) {
      s <- rep(1 / sqrt(length(p$K)), length(p$K))
      stochastic <- s * u * (p$r / p$K) * (p$K - c(p$A %*% u))
      return(stochastic)
    }
    # integration time steps
    time_step <- times[2] - times[1]
    # numerical integration
    sol <- sde.solve(f = f, g = g, u0 = N, tspan = range(times), p = pars, saveat = time_step)
    out <- as.data.frame(cbind(sol$t, sol$u))
    names(out) <- paste("time", 1:length(pars$K))
  }
  return(out)
}
