# Read LV simulation results and make plots for the SI version of
# Figure 2 (alternative species pools and LV models)

# clean wd, load packages and functions ------------------------------
rm(list = ls(all = TRUE))
if(!require(ggplot2)) {install.packages("ggplot2"); library(ggplot2)}
if(!require(latex2exp)) {install.packages("latex2exp"); library(latex2exp)}
if(!require(reshape2)) {install.packages("reshape2"); library(reshape2)}
if(!require(dplyr)) {install.packages("dplyr"); library(dplyr)}
if(!require(viridis)) {install.packages("viridis"); library(viridis)}
if(!require(ggpubr)) {install.packages("ggpubr"); library(ggpubr)}
if(!require(deSolve)) {install.packages("deSolve"); library(deSolve)}
source("code/functions/sample_niche_overlap_matrix.R")
source("code/functions/scale_matrix.R")
source("code/functions/scaling_factor.R")
source("code/functions/lv_pruning.R")
source("code/functions/lotka_volterra.R")
source("code/functions/hypersphere_sampling.R")

# define plot scenario ------------------------------
pool <- "niche_stochastic"

# random species pool ------------------------------
if (pool == "random") {
  # read results file
  results_df <- read.csv("results/fig2/random/results_K_lv_pruning_theoretical_matrices_n20_randompool_ncomm10_nsamples500_param1.csv")
  # plotting
  fig <- ggplot(data = results_df) +
    geom_point(aes(x = n_sp_pruned/n_sp_pool, y = omega_percentile), 
               size = 2.2, alpha = 0.7) +
    geom_hline(yintercept = 0.5, size = 1.2, linetype = "dashed") +
    scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
    xlab(TeX("Fraction of surviving species from the pool $\\left( \\frac{s}{S} \\right)$")) +
    ylab(expression(atop("Percentile rank of expected probability",
                         "of species persistence" ~ (omega(alpha))))) +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA, size = 1),
          axis.text.y = element_text(size = 14),
          axis.title = element_text(size = 18),
          axis.text.x = element_text(size = 14),
          strip.text.x = element_text(size = 14),
          strip.text.y = element_text(size = 14),
          legend.position = "none",
          plot.margin = unit(c(30, 30, 30, 30), "pt"))
  # save complete figure
  ggsave("figs/fig2/fig2_random.pdf", fig, 
         width = 16, height = 13, units = "cm")
}

# niche with different widths species pool ------------------------------
if (pool == "niche_different_width") {
  # read results file
  results_df <- read.csv("results/fig2/niche_different_width/results_K_lv_pruning_theoretical_matrices_n20_nichepool_ncomm10_nsamples500_param0.05.csv")
  # plotting
  fig <- ggplot(data = results_df) +
    geom_point(aes(x = n_sp_pruned/n_sp_pool, y = omega_percentile), 
               size = 2.2, alpha = 0.7) +
    geom_hline(yintercept = 0.5, size = 1.2, linetype = "dashed") +
    scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
    xlab(TeX("Fraction of surviving species from the pool $\\left( \\frac{s}{S} \\right)$")) +
    ylab(expression(atop("Percentile rank of expected probability",
                         "of species persistence" ~ (omega(alpha))))) +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA, size = 1),
          axis.text.y = element_text(size = 14),
          axis.title = element_text(size = 18),
          axis.text.x = element_text(size = 14),
          strip.text.x = element_text(size = 14),
          strip.text.y = element_text(size = 14),
          legend.position = "none",
          plot.margin = unit(c(30, 30, 30, 30), "pt"))
  # save complete figure
  ggsave("figs/fig2/fig2_niche_different_width.pdf", fig, 
         width = 16, height = 13, units = "cm")
}

# niche species pool ------------------------------
if (pool == "niche_normalized" | pool == "niche_typeII" | pool == "niche_stochastic") {
  # read species pool file
  if (pool == "niche_normalized") {
    results_pool_df <- read.csv("results/fig2/niche_normalized/results_species_pool_theoretical_matrices_n20_nichepool_ncomm500.csv")
  } else {
    results_pool_df <- read.csv("results/fig2/niche/results_species_pool_theoretical_matrices_n20_nichepool_ncomm500.csv")
  }
  results_pool_df$param[results_pool_df$param == 0.025] <- "Low"
  results_pool_df$param[results_pool_df$param == 0.05] <- "Intermediate"
  results_pool_df$param[results_pool_df$param == 0.1] <- "High"
  results_pool_df$param <- factor(results_pool_df$param, 
                                  levels = c("Low", "Intermediate", "High"))
  # read results files
  path <- paste("results/fig2/", pool, sep = "")
  all_files <- dir(path)
  all_files_sub <- all_files[grep("lv_pruning", all_files)]
  file_paths <- paste(path, all_files_sub, sep = "/")
  results_list <- lapply(file_paths, read.csv, header = TRUE)
  names(results_list) <- gsub('.{4}$', '', sub("_.*", "", sub(".*param", "", all_files_sub)))
  results_df <- bind_rows(results_list)
  results_df$param <- as.numeric(results_df$param)
  results_df$param[results_df$param == 0.025] <- "Low"
  results_df$param[results_df$param == 0.05] <- "Intermediate"
  results_df$param[results_df$param == 0.1] <- "High"
  results_df$param <- factor(results_df$param, 
                             levels = c("Low", "Intermediate", "High"))
  # species pool plot
  fig_A <- ggplot() +
    geom_boxplot(data = results_pool_df, aes(x = param, y = omega_pool, group = param, color = param), 
                 size = 1, outlier.size = 2, outlier.alpha = 0.7) +
    scale_color_viridis(discrete = TRUE) +
    xlab(TeX("Mean competition strength ($\\rho(A)$)")) +
    ylab(expression(atop("Expected probability of",
                         "species persistence" ~ (omega(A))))) +
    ggtitle("Species pool (A)") + 
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA, size = 1),
          axis.text.y = element_text(size = 14),
          axis.title = element_text(size = 18),
          axis.text.x = element_text(size = 14),
          plot.title = element_text(size = 20, hjust = 0.5),
          legend.position = "none",
          plot.margin = unit(c(200, 50, 200, 50), "pt"))
  # observed vs potential omega plot
  fig_B <- ggplot(data = results_df) +
    geom_point(aes(x = n_sp_pruned/n_sp_pool, y = omega_percentile, color = param), 
               size = 2.5, alpha = 0.7) +
    geom_hline(yintercept = 0.5, size = 1, linetype = "dashed") +
    facet_wrap(~ param, nrow = 3) +
    scale_color_viridis(discrete = TRUE) +
    scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
    ggtitle(TeX("Realized combination ($\\alpha$)")) +
    xlab(TeX("Fraction of surviving species from the pool $\\left( \\frac{s}{S} \\right)$")) +
    ylab(expression(atop("Percentile rank of expected probability",
                         "of species persistence" ~ (omega(alpha))))) +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA, size = 1),
          axis.text.y = element_text(size = 14),
          axis.title = element_text(size = 18),
          axis.text.x = element_text(size = 14),
          strip.text.x = element_text(size = 14),
          strip.text.y = element_text(size = 14),
          plot.title = element_text(size = 20, hjust = 0.5),
          legend.position = "none",
          plot.margin = unit(c(10, 20, 10, 20), "pt"))
  # save complete figure
  fig <- ggarrange(ggarrange(fig_A, fig_B, ncol = 2, labels = c("A", "B"), 
                             font.label = list(size = 28, face = "bold")),
                   heights = c(1, 1.8))
  ggsave(paste("figs/fig2/fig2_", pool, ".pdf", sep = ""), 
         fig, width = 32, height = 26, units = "cm")
}

# bipartite species pool ------------------------------
if (pool == "bipartite") {
  # read species pool file
  results_pool_df <- read.csv(paste("results/fig2/", pool, 
                                    "/results_species_pool_theoretical_matrices_n20_bipartitepool_ncomm500.csv",
                                    sep = ""))
  results_pool_df$param[results_pool_df$param == 0.1] <- "Low"
  results_pool_df$param[results_pool_df$param == 0.2] <- "Intermediate"
  results_pool_df$param[results_pool_df$param == 0.3] <- "High"
  results_pool_df$param <- factor(results_pool_df$param, 
                                  levels = c("Low", "Intermediate", "High"))
  # read results files
  path <- paste("results/fig2/", pool, sep = "")
  all_files <- dir(path)
  all_files_sub <- all_files[grep("lv_pruning", all_files)]
  file_paths <- paste(path, all_files_sub, sep = "/")
  results_list <- lapply(file_paths, read.csv, header = TRUE)
  names(results_list) <- gsub('.{4}$', '', sub("_.*", "", sub(".*param", "", all_files_sub)))
  results_df <- bind_rows(results_list)
  results_df$param <- as.numeric(results_df$param)
  results_df$param[results_df$param == 0.1] <- "Low"
  results_df$param[results_df$param == 0.2] <- "Intermediate"
  results_df$param[results_df$param == 0.3] <- "High"
  results_df$param <- factor(results_df$param,
                             levels = c("Low", "Intermediate", "High"))
  # species pool plot
  fig_A <- ggplot() +
    geom_boxplot(data = results_pool_df, aes(x = param, y = omega_pool, group = param, color = param), 
                 size = 1, outlier.size = 2, outlier.alpha = 0.7) +
    scale_color_viridis(discrete = TRUE) +
    xlab(TeX("Mean connectance ($p_c$)")) +
    ylab(expression(atop("Expected probability of",
                         "species persistence" ~ (omega(A))))) +
    ggtitle("Species pool (A)") + 
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA, size = 1),
          axis.text.y = element_text(size = 14),
          axis.title = element_text(size = 18),
          axis.text.x = element_text(size = 14),
          plot.title = element_text(size = 20, hjust = 0.5),
          legend.position = "none",
          plot.margin = unit(c(200, 50, 200, 50), "pt"))
  # observed vs potential omega plot
  fig_B <- ggplot(data = results_df) +
    geom_point(aes(x = n_sp_pruned/n_sp_pool, y = omega_percentile, color = param), 
               size = 2.5, alpha = 0.7) +
    geom_hline(yintercept = 0.5, size = 1, linetype = "dashed") +
    facet_wrap(~ param, nrow = 3) +
    scale_color_viridis(discrete = TRUE) +
    scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
    ggtitle(TeX("Realized combination ($\\alpha$)")) +
    xlab(TeX("Fraction of surviving species from the pool $\\left( \\frac{s}{S} \\right)$")) +
    ylab(expression(atop("Percentile rank of expected probability",
                         "of species persistence" ~ (omega(alpha))))) +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA, size = 1),
          axis.text.y = element_text(size = 14),
          axis.title = element_text(size = 18),
          axis.text.x = element_text(size = 14),
          strip.text.x = element_text(size = 14),
          strip.text.y = element_text(size = 14),
          plot.title = element_text(size = 20, hjust = 0.5),
          legend.position = "none",
          plot.margin = unit(c(10, 20, 10, 20), "pt"))
  # save complete figure
  fig <- ggarrange(ggarrange(fig_A, fig_B, ncol = 2, labels = c("A", "B"), 
                             font.label = list(size = 28, face = "bold")),
                   heights = c(1, 1.8))
  ggsave("figs/fig2/fig2_bipartite.pdf", fig, width = 32, height = 26, units = "cm")
}
