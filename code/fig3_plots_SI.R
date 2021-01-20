# Read empirical data results and make plots for the SI version of
# Figure 3 (alternative methods of reconstructing species pools and
# inferring interactions)

# clean wd, load packages and functions ------------------------------
rm(list = ls(all = TRUE))
if(!require(ggplot2)) {install.packages("ggplot2"); library(ggplot2)}
if(!require(latex2exp)) {install.packages("latex2exp"); library(latex2exp)}
if(!require(reshape2)) {install.packages("reshape2"); library(reshape2)}
if(!require(dplyr)) {install.packages("dplyr"); library(dplyr)}
if(!require(viridis)) {install.packages("viridis"); library(viridis)}
if(!require(ggpubr)) {install.packages("ggpubr"); library(ggpubr)}
if(!require(deSolve)) {install.packages("deSolve"); library(deSolve)}
source("code/functions/scale_matrix.R")
source("code/functions/scaling_factor.R")
source("code/functions/build_metanet.R")
source("code/functions/remove_nodes_no_links.R")

# read matrix files ------------------------------
# all files
all_files <- dir("data/matrices/")
# define whether to use year, stage, or year and stage
dataset <- "year_stage"
# which interaction inferrence method to use
inferrence <- "scaling_factor"
# files for stage and year
if (dataset == "year_stage") {
  nets <- all_files[grep("middle_2011", dir("data/matrices/"))]
}
# files for year
if (dataset == "year") {
  nets <- all_files[grep("2011", dir("data/matrices/"))]
}
# files for stage
if (dataset == "stage") {
  nets <- all_files[grep("middle", dir("data/matrices/"))]
}
# complete file paths
net_files <- paste("data/matrices/", nets, sep = "")
# read files
mats <- lapply(net_files, read.table, header = TRUE, row.names = 1)
mats <- lapply(mats, as.matrix)
# change matrices to binary
mats <- lapply(mats, function(mat) {mat[mat > 0] <- 1; return(mat)})
# remove species with no interactions
mats <- lapply(mats, remove_nodes_no_links)

# plot species pool matrix ------------------------------
# building species pool
metanet <- build_metanet(mats)
# building observed matrix for stage and year
mat <- mats[[which(nets == "7_middle_2011.txt")]]
# extract plant labels from matrix
plant_labels <- rownames(mat)
# create potential matrix
potential_net <- metanet[plant_labels, ]
# remove species with no interactions
potential_net <- remove_nodes_no_links(potential_net)
# create monopartite projection (herbivore competition matrix)
proj_mat <- t(potential_net) %*% potential_net
# plot for column sum inferrence
if (inferrence == "col_sum") {
  # rescale projection
  comp_mat <- scale_matrix(proj_mat, remove_zero = TRUE, 
                           scale = "col_sum", p = 1,
                           tol = 10^-8, diag = "one")
  # building data frame by melting matrix
  molten_mat <- melt(comp_mat)
  names(molten_mat) <- c("row", "column", "interaction")
  # plotting
  fig_A <- ggplot(data = subset(molten_mat, interaction < 1), 
                  aes(x = column, y = row, fill = interaction)) + 
    geom_tile() +
    coord_equal() +
    scale_fill_viridis(option = "plasma", name = "Competition\nstrength",
                       limits = c(min(comp_mat[row(comp_mat)!=col(comp_mat)]), 
                                  max(comp_mat[row(comp_mat)!=col(comp_mat)]))) +
    ylim(rev(levels(molten_mat$row))) +
    scale_x_discrete(position = "top") +
    ggtitle("Species pool (A) for plot 5") + 
    ylab("Herbivores") +
    xlab("Herbivores") +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA, size = 1),
          legend.title = element_text(size = 16),
          legend.text = element_text(size = 14),
          legend.key.size = unit(0.5, "cm"),
          plot.title = element_text(size = 20, hjust = 0.5),
          axis.title = element_text(size = 18),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.x.top = element_blank(),
          axis.text.y = element_blank(),
          plot.margin = unit(c(10, 10, 10, 10), "pt"))
}
# plot for scaling factor inferrence
if (inferrence == "scaling_factor") {
  # rescale projection
  comp_mat <- scale_matrix(proj_mat, remove_zero = TRUE, 
                           scale = "scaling_factor", p = 0.5,
                           tol = 10^-8, diag = "preserve")
  # building data frame by melting matrix
  molten_mat <- melt(comp_mat)
  names(molten_mat) <- c("row", "column", "interaction")
  # plotting
  fig_A <- ggplot(data = molten_mat, 
                  aes(x = column, y = row, fill = interaction)) + 
    geom_tile() +
    coord_equal() +
    scale_fill_viridis(option = "plasma", name = "Competition\nstrength") +
    ylim(rev(levels(molten_mat$row))) +
    scale_x_discrete(position = "top") +
    ggtitle("Species pool (A) for plot 5") + 
    ylab("Herbivores") +
    xlab("Herbivores") +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA, size = 1),
          legend.title = element_text(size = 16),
          legend.text = element_text(size = 14),
          legend.key.size = unit(0.5, "cm"),
          plot.title = element_text(size = 20, hjust = 0.5),
          axis.title = element_text(size = 18),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.x.top = element_blank(),
          axis.text.y = element_blank(),
          plot.margin = unit(c(10, 10, 10, 10), "pt"))
}

# plot observed competition matrix ------------------------------
# create monopartite projection (herbivore competition matrix)
proj_mat <- t(mat) %*% mat
# plot for column sum inferrence
if (inferrence == "col_sum") {
  # rescale projection
  comp_mat <- scale_matrix(proj_mat, remove_zero = TRUE, 
                           scale = "col_sum", p = 1,
                           tol = 10^-8, diag = "one")
  # building data frame by melting matrix
  molten_mat <- melt(comp_mat)
  names(molten_mat) <- c("row", "column", "interaction")
  # plotting
  fig_B <- ggplot(data = subset(molten_mat, interaction < 1), 
                  aes(x = column, y = row, fill = interaction)) + 
    geom_tile() +
    coord_equal() +
    scale_fill_viridis(option = "plasma", name = "Competition\nstrength",
                       limits = c(min(comp_mat[row(comp_mat)!=col(comp_mat)]), 
                                  max(comp_mat[row(comp_mat)!=col(comp_mat)]))) +
    ylim(rev(levels(molten_mat$row))) +
    scale_x_discrete(position = "top") +
    ggtitle(expression(atop("Observed combination" ~ (alpha),
                            "for plot 5"))) +
    ylab("Herbivores") +
    xlab("Herbivores") +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA, size = 1),
          legend.title = element_text(size = 16),
          legend.text = element_text(size = 14),
          legend.key.size = unit(0.5, "cm"),
          plot.title = element_text(size = 20, hjust = 0.5),
          axis.title = element_text(size = 18),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.x.top = element_blank(),
          axis.text.y = element_blank(),
          plot.margin = unit(c(30, 30, 30, 30), "pt"))
}
# plot for scaling_factor inferrence
if (inferrence == "scaling_factor") {
  # rescale projection
  comp_mat <- scale_matrix(proj_mat, remove_zero = TRUE, 
                           scale = "scaling_factor", p = 0.5,
                           tol = 10^-8, diag = "preserve")
  # building data frame by melting matrix
  molten_mat <- melt(comp_mat)
  names(molten_mat) <- c("row", "column", "interaction")
  # plotting
  fig_B <- ggplot(data = molten_mat, 
                  aes(x = column, y = row, fill = interaction)) + 
    geom_tile() +
    coord_equal() +
    scale_fill_viridis(option = "plasma", name = "Competition\nstrength") +
    ylim(rev(levels(molten_mat$row))) +
    scale_x_discrete(position = "top") +
    ggtitle(expression(atop("Observed combination" ~ (alpha),
                            "for plot 5"))) +
    ylab("Herbivores") +
    xlab("Herbivores") +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA, size = 1),
          legend.title = element_text(size = 16),
          legend.text = element_text(size = 14),
          legend.key.size = unit(0.5, "cm"),
          plot.title = element_text(size = 20, hjust = 0.5),
          axis.title = element_text(size = 18),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.x.top = element_blank(),
          axis.text.y = element_blank(),
          plot.margin = unit(c(30, 30, 30, 30), "pt"))
}

# read results file ------------------------------
# file for stage and year
if (dataset == "year_stage") {
  # file for column sum inferrence
  if (inferrence == "col_sum") {
    results_df <- read.csv(file = "results/fig3/results_2007-2017_year_stage_metanetwork_col_sum.csv", 
                           as.is = TRUE)
  }
  # file for scaling factor inferrence
  if (inferrence == "scaling_factor") {
    results_df <- read.csv(file = "results/fig3/results_2007-2017_year_stage_metanetwork_scaling_factor.csv", 
                           as.is = TRUE)
  }
}
# file for year
if (dataset == "year") {
  results_df <- read.csv(file = "results/fig3/results_2007-2017_year_metanetwork_col_sum.csv",
                         as.is = TRUE)
}
# file for stage
if (dataset == "stage") {
  results_df <- read.csv(file = "results/fig3/results_2007-2017_stage_metanetwork_col_sum.csv",
                         as.is = TRUE)
}
# changing factors
results_df$stage[results_df$stage == "initial"] <- "Initial"
results_df$stage[results_df$stage == "middle"] <- "Middle"
results_df$stage[results_df$stage == "late"] <- "Late"
results_df$stage <- factor(results_df$stage, 
                           levels = c("Initial", "Middle", "Late"))
results_df$year <- factor(results_df$year,
                          levels = c("2007", "2008", "2009", "2010", 
                                     "2011", "2012", "2013", "2014", "2016", "2017"))
results_df$plot <- factor(results_df$plot,
                          levels = c("1", "2", "3", "4", "5", "6",
                                     "7", "8", "9"))
# computing relative omega and competition strength
results_df$omega_std <- (results_df$omega_real - results_df$omega_mean) / 
  results_df$omega_mean
results_df$comp_strength_std <- (results_df$comp_strength_real - results_df$comp_strength_mean) / 
  results_df$comp_strength_mean

# probability of persistence of species pool ------------------------------
# plotting
fig_C <- ggplot() +
  geom_boxplot(data = results_df, aes(x = stage, y = omega_pool, group = stage, color = stage), 
               size = 1, outlier.size = 2.5) +
  scale_color_viridis(discrete = TRUE) +
  xlab("Successional stage") +
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
        plot.margin = unit(c(80, 50, 80, 50), "pt"))

# probability of persistence of observed combinations ------------------------------
# plotting
fig_D <- ggplot(data = results_df) +
  geom_point(aes(x = n_sp_real/n_sp_pot, y = omega_percentile, color = stage), 
             size = 3.5) +
  geom_hline(yintercept = 0.5, size = 1, linetype = "dashed") +
  facet_wrap(~ stage, nrow = 3) +
  scale_color_viridis(discrete = TRUE) +
  ggtitle(TeX("Observed combination ($\\alpha$)")) +
  xlab(TeX("Fraction of observed species from the pool $\\left( \\frac{s}{S} \\right)$")) +
  ylab(expression(atop("Percentile rank of expected probability",
                       "of species persistence" ~ (omega(alpha))))) +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1),
        legend.position = "none",
        strip.text.x = element_text(size = 14),
        strip.text.y = element_text(size = 14),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title = element_text(size = 18),
        plot.title = element_text(size = 20, hjust = 0.5),
        plot.margin = unit(c(10, 20, 10, 20), "pt"))

# save complete figure ------------------------------
fig <- ggarrange(ggarrange(fig_A, fig_B, ncol = 2, labels = c("A", ""), 
                           font.label = list(size = 28, face = "bold")),
                 ggarrange(fig_C, fig_D, ncol = 2, labels = c("B", "C"),
                           font.label = list(size = 28, face = "bold")), 
                 heights = c(1, 1.8), nrow = 2)
ggsave(paste("figs/fig3/fig3_", dataset, "_", inferrence, ".pdf", sep = ""), 
       fig, width = 32, height = 26, units = "cm")
