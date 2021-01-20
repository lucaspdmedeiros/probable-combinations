# Computes the fraction of LV simulations in which each individual 
# species and each combination of species survived and plots the results 

# cleaning wd, loading functions and packages ------------------------------
rm(list = ls(all = TRUE))
if(!require(plotly)) {install.packages("plotly"); library(plotly)}
if(!require(tidyverse)) {install.packages("tidyverse"); library(tidyverse)}
if(!require(tidyr)) {install.packages("tidyr"); library(tidyr)}
if(!require(plyr)) {install.packages("plyr"); library(plyr)}
if(!require(viridis)) {install.packages("viridis"); library(viridis)}
if(!require(ggpubr)) {install.packages("ggpubr"); library(ggpubr)}
if(!require(mvtnorm)) {install.packages("mvtnorm"); library(mvtnorm)}
if(!require(matrixcalc)) {install.packages("matrixcalc"); library(matrixcalc)}
if(!require(processx)) {install.packages("processx"); library(processx)}
source("code/functions/omega.R")
# whether to save plots
save_plots <- TRUE

# reading results file ------------------------------
mat_name <- "saavedra_etal_2017_7B"
df <- read.csv(paste("results/fig1/lv_surv_sp_", mat_name, ".csv", sep = ""))
sp_names <- paste("sp", 1:3, sep = "")
# change combination names
df$surv_sp <- as.character(df$surv_sp)
df$surv_sp[df$surv_sp == "1_2"] <- "1, 2"
df$surv_sp[df$surv_sp == "1_3"] <- "1, 3"
df$surv_sp[df$surv_sp == "2_3"] <- "2, 3"
df$surv_sp[df$surv_sp == "1_2_3"] <- "1, 2, 3"
df$surv_sp <- factor(df$surv_sp, levels = c("1", "2", "3", "1, 2", "1, 3", "2, 3", "1, 2, 3"))

# Figure 1B ------------------------------
# creating color palette
pal_set1 <- c("#E41A1C", "#377EB8", "#FFFF33", "#984EA3", "#FF7F00", "#4DAF4A", "#A65628")
# plot distribution of simulation outcomes on sphere
fig_B <- plot_ly(df, x = ~K2, y = ~K1, z = ~K3, color = ~surv_sp, colors = pal_set1,
                 showlegend = FALSE)
fig_B <- fig_B %>% add_markers(marker = list(size = 6, 
                                             line = list(color = "black", width = 1)))
fig_B <- fig_B %>% layout(scene = list(xaxis = list(title = "K<sub>2</sub>",
                                                    titlefont = list(size = 32), 
                                                    tickfont = list(size = 18)),
                                       yaxis = list(title = "K<sub>1</sub>",
                                                    titlefont = list(size = 32), 
                                                    tickfont = list(size = 18)),
                                       zaxis = list(title = "K<sub>3</sub>",
                                                    titlefont = list(size = 32), 
                                                    tickfont = list(size = 18))))
# save plot
if (save_plots) {
  orca(fig_B, "figs/fig1/fig1_B.pdf", format = "pdf",
       width = 800, height = 800)
}

# Figure 1C ------------------------------
# building data frame for plotting with sp combinations
sp_counts <- table(df$surv_sp)
comb_df <- data.frame(matrix(unlist(sp_counts), nrow = length(sp_counts), byrow = TRUE))
colnames(comb_df) <- "number_simulations"
comb_df$fraction_simulations <- comb_df$number_simulations / sum(comb_df$number_simulations)
comb_df$surviving_sp <- names(sp_counts)
# changing factor levels
comb_df$surviving_sp <- factor(comb_df$surviving_sp, 
                               levels = c("1", "2", "3", "1, 2", "1, 3", "2, 3", "1, 2, 3"))
comb_df$n_sp <- factor(c(1, 1, 1, 2, 2, 2, 3))
# plot probability of number of surviving species
fig_C <- ggplot(comb_df, aes(x = n_sp, y = fraction_simulations, fill = surviving_sp)) +
  geom_bar(stat = "identity", width = 0.5, size = 0.7, color = "black") +
  xlab("Number of surviving species") +
  ylab("Fraction of simulations") +
  scale_fill_manual(values = pal_set1) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1),
        axis.text.y = element_text(size = 20),
        plot.title = element_text(size = 26, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 24),
        axis.text.x = element_text(size = 20),
        legend.position = "none")
# save plot
if (save_plots) {
  ggsave("figs/fig1/fig1_C.pdf", fig_C, width = 13, height = 13, units = "cm")
}

# Figure 1D ------------------------------
# data frame for subsets with 2 species
comb_df_2_sp <- subset(comb_df, surviving_sp == "1, 2" | surviving_sp ==  "1, 3" | 
                         surviving_sp == "2, 3")
# computing omegas (size of feasibility domain)
mat <- as.matrix(read.table("data/fig1/saavedra_etal_2017_7B.txt"))
mat_1_2 <- mat[c(1, 2), c(1, 2)]
omega_1_2 <- 2 * Omega(mat_1_2)[[2]]
mat_1_3 <- mat[c(1, 3), c(1, 3)]
omega_1_3 <- 2 * Omega(mat_1_3)[[2]]
mat_2_3 <- mat[c(2, 3), c(2, 3)]
omega_2_3 <- 2 * Omega(mat_2_3)[[2]]
comb_df_2_sp$omega <- c(omega_1_2, omega_1_3, omega_2_3)
# plot for simulation probabilities
pal_2_sp <- c("#984EA3", "#FF7F00", "#4DAF4A")
fig <- ggplot(comb_df_2_sp, aes(x = surviving_sp, y = fraction_simulations, color = surviving_sp)) +
  geom_point(size = 5) +
  geom_hline(yintercept = mean(comb_df_2_sp$fraction_simulations), linetype = "dashed", 
             size = 1.5) +
  xlab("2-species combination") +
  ylab("Probability of occurrence") +
  scale_color_manual(values = pal_2_sp) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1),
        axis.text.y = element_text(size = 20),
        plot.title = element_text(size = 26, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 24),
        axis.text.x = element_text(size = 20),
        legend.position = "none")
# plot for analytical probabilities (omega)
fig_D <- ggplot(comb_df_2_sp, aes(x = surviving_sp, y = omega, color = surviving_sp)) +
  geom_hline(yintercept = median(comb_df_2_sp$omega), linetype = "dashed", 
             size = 1.5) +
  geom_point(size = 5) +
  xlab("2-species combination") +
  ylab("Expected probability of\nspecies persistence") +
  scale_color_manual(values = pal_2_sp) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1),
        axis.text.y = element_text(size = 20),
        plot.title = element_text(size = 26, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 24),
        axis.text.x = element_text(size = 20),
        legend.position = "none")
if (save_plots) {
  ggsave("figs/fig1/fig1_D.pdf", fig_D, width = 13, height = 13, units = "cm")
}
