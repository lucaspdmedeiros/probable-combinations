# Read LV simulation results and make plots for Figure 2

# clean wd and set seed ------------------------------
rm(list = ls(all = TRUE))
# set random seed
set.seed(8)

# load packages and functions ------------------------------
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

# heatmap of species pool interaction matrix ------------------------------
# define matrix size
n_sp <- 20
# species pool type
pool <- "niche"
# whether to normalize matrices
normalize <- FALSE
# define average niche overlap
rho <- 0.1
# sample matrices
mat <- sample_niche_overlap_matrix(n_sp, rho, same_width = TRUE)
# add species labels
rownames(mat) <- paste("sp", 1:nrow(mat), sep = "")
colnames(mat) <- paste("sp", 1:nrow(mat), sep = "")
# normalize matrix
if (normalize) {
  mat <- scale_matrix(mat, remove_zero = TRUE, 
                      scale = "col_sum", p = 1,
                      tol = 10^-8, diag = "one")
}
# building data frame by melting matrix
molten_mat <- melt(mat)
names(molten_mat) <- c("row", "column", "interaction")
# plotting
fig_A <- ggplot(data = subset(molten_mat, interaction < 1), 
                aes(x = column, y = row, fill = interaction)) + 
  geom_tile() +
  coord_equal() +
  scale_fill_viridis(option = "plasma", name = "Competition\nstrength",
                     limits = c(min(mat[row(mat)!=col(mat)]), 
                                max(mat[row(mat)!=col(mat)]))) +
  ggtitle("Species pool (A)") + 
  ylim(rev(levels(molten_mat$row))) +
  scale_x_discrete(position = "top") +
  ylab("") +
  xlab("") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1),
        plot.title = element_text(size = 18, hjust = 0.5),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.key.size = unit(0.5, "cm"),
        axis.text.x.top = element_text(size = 9, angle = 90, vjust = 0.5),
        axis.text.y = element_text(size = 9),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        plot.margin = unit(c(10, 10, 10, 10), "pt"))

# heatmap of pruned interaction matrix ------------------------------
# create simulation parameters
N <- rep(1, n_sp)
r <- rep(1, n_sp)
K <- hypersphere_sampling(n_sp, positive = TRUE, within = FALSE)
times <- seq(0, 100, 0.01)
# setting coefficients negative
mat <- - mat
# simulate and obtain indexes of surviving species
which_surv_sp <- lv_pruning(N, r, K, mat, times, extinct_tol = 10^-2, formalism = "K")
# obtain pruned matrix
mat_pruned <- - mat[which_surv_sp, which_surv_sp]
# building data frame by melting matrix
molten_mat <- melt(mat_pruned)
names(molten_mat) <- c("row", "column", "interaction")
# plotting
fig_B <- ggplot(data = subset(molten_mat, interaction < 1), 
                aes(x = column, y = row, fill = interaction)) + 
  geom_tile() +
  coord_equal() +
  scale_fill_viridis(option = "plasma", name = "Competition\nstrength",
                     limits = c(min(mat_pruned[row(mat_pruned)!=col(mat_pruned)]), 
                                max(mat_pruned[row(mat_pruned)!=col(mat_pruned)]))) +
  ggtitle(TeX("Realized combination ($\\alpha$)")) +
  ylim(rev(levels(molten_mat$column))) +
  scale_x_discrete(position = "top") +
  ylab("") +
  xlab("") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1),
        plot.title = element_text(size = 18, hjust = 0.5),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.key.size = unit(0.5, "cm"),
        axis.text.x.top = element_text(size = 10, angle = 90, vjust = 0.5),
        axis.text.y = element_text(size = 10),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        plot.margin = unit(c(30, 30, 30, 30), "pt"))

# probability of persistence for different competition levels ------------------------------
# read results file
results_pool_df <- read.csv(paste("results/fig2/", pool, 
                                  "/results_species_pool_theoretical_matrices_n20_nichepool_ncomm500.csv",
                                  sep = ""))
results_pool_df$param[results_pool_df$param == 0.025] <- "Low"
results_pool_df$param[results_pool_df$param == 0.05] <- "Intermediate"
results_pool_df$param[results_pool_df$param == 0.1] <- "High"
results_pool_df$param <- factor(results_pool_df$param, 
                                levels = c("Low", "Intermediate", "High"))
# plotting
fig_C <- ggplot() +
  geom_boxplot(data = results_pool_df, aes(x = param, y = omega_pool, group = param, color = param), 
               size = 1, outlier.size = 2, outlier.alpha = 0.7) +
  scale_color_viridis(discrete = TRUE) +
  xlab(TeX("Mean competition strength ($\\rho(A)$)")) +
  ylab(expression(atop("Expected probability of",
                       "species persistence" ~ (omega(A))))) +
  ggtitle("Species pool (A)") + 
  scale_y_continuous(limits = c(0, 0.75), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1),
        axis.text.y = element_text(size = 14),
        axis.title = element_text(size = 18),
        axis.text.x = element_text(size = 14),
        plot.title = element_text(size = 20, hjust = 0.5),
        legend.position = "none",
        plot.margin = unit(c(80, 50, 80, 50), "pt"))

# probability of persistence of pruned matrices for different competition levels ------------------------------
# read results files
path <- paste("results/fig2/", pool, sep = "")
all_files <- dir(path)
all_files_sub <- all_files[grep("lv_pruning", all_files)]
file_paths <- paste(path, all_files_sub, sep = "/")
results_list <- lapply(file_paths, read.csv, header = TRUE)
names(results_list) <- gsub('.{4}$', '', sub("_.*", "", sub(".*param", "", all_files_sub)))
results_df <- bind_rows(results_list)
results_df$param <- as.numeric(results_df$param)

# observed vs potential omega ------------------------------
# computing relative difference from mean
results_df$omega_comm_std <- (results_df$omega_obs - results_df$omega_expect_comm) / 
  results_df$omega_expect_comm
results_df$param[results_df$param == 0.025] <- "Low"
results_df$param[results_df$param == 0.05] <- "Intermediate"
results_df$param[results_df$param == 0.1] <- "High"
results_df$param <- factor(results_df$param,
                           levels = c("Low", "Intermediate", "High"))
# plotting
fig_D <- ggplot(data = results_df) +
  geom_point(aes(x = n_sp_pruned/n_sp_pool, y = omega_percentile, color = param), 
             size = 2.5, alpha = 0.7) +
  geom_hline(yintercept = 0.5, size = 1.5, linetype = "dashed") +
  facet_wrap(~ param, nrow = 3) +
  scale_color_viridis(discrete = TRUE) +
  ggtitle(TeX("Realized combination ($\\alpha$)")) +
  scale_x_continuous(limits = c(0.25, 1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
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
        plot.title = element_text(size = 20, hjust = 0.5),
        legend.position = "none",
        plot.margin = unit(c(10, 20, 10, 20), "pt"))

# save complete figure ------------------------------
fig <- ggarrange(ggarrange(fig_A, fig_B, ncol = 2, labels = c("A", ""), 
                           font.label = list(size = 28, face = "bold")),
                 ggarrange(fig_C, fig_D, ncol = 2, labels = c("B", "C"),
                           font.label = list(size = 28, face = "bold")), 
                 heights = c(1, 1.8), nrow = 2)
ggsave("figs/fig2/fig2.pdf", fig, width = 32, height = 26, units = "cm")
