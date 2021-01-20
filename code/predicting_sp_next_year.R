# Simulates Lotka-Volterra dynamics parameterized with the empirical data from
# a given year to predict the species composition of the next year and returns
# a plot of the true positive rate vs false positive rate

# cleaning wd, loading functions and packages ------------------------------
rm(list = ls(all = TRUE))
if(!require(plyr)) {install.packages("plyr"); library(plyr)}
if(!require(mvtnorm)) {install.packages("mvtnorm"); library(mvtnorm)}
if(!require(matrixcalc)) {install.packages("matrixcalc"); library(matrixcalc)}
if(!require(deSolve)) {install.packages("deSolve"); library(deSolve)}
if(!require(ggplot2)) {install.packages("ggplot2"); library(ggplot2)}
if(!require(ggpubr)) {install.packages("ggpubr"); library(ggpubr)}
source("code/functions/lotka_volterra.R")
source("code/functions/lv_pruning.R")
source("code/functions/omega.R")
source("code/functions/scale_matrix.R")
source("code/functions/scaling_factor.R")
source("code/functions/remove_nodes_no_links.R")

# read network files ------------------------------
# file paths
path <- "data/matrices"
all_files <- dir(path)
all_files_sub <- all_files[grep(".txt", all_files)]
file_paths <- paste(path, all_files_sub, sep = "/")
# read files
mats <- lapply(file_paths, read.table, header = TRUE, row.names = 1)
mats <- lapply(mats, as.matrix)
# change networks to binary
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
# vectors to store results
tpr <- c()
fpr <- c()
mat_name_curr <- c()
n_sp_focus <- c()
mat_name_next <- c()
n_sp_pred <- c()

# loop over all plots ------------------------------
for (i in 1:length(unique(plot))) {
  files_id <- which(plot == i)

  # loop over all matrices from the plot ------------------------------
  for (j in 1:(length(files_id)-1)) {
    # build matrices
    mat_id <- files_id[j]
    mat_curr <- as.matrix(mats[[mat_id]])
    mat_next <- as.matrix(mats[[mat_id + 1]])
    # matrix to simulate LV
    print(gsub("\\..*", "", all_files_sub[mat_id]))
    # matrix to predict species persistence
    print(gsub("\\..*", "", all_files_sub[mat_id + 1]))
    cat("\n")
    # create monopartite projection (herbivore competition matrix)
    proj_mat_curr <- t(mat_curr) %*% mat_curr
    proj_mat_next <- t(mat_next) %*% mat_next
    # compute number of species of both matrices
    n_sp_curr <- nrow(proj_mat_curr)
    n_sp_focus <- c(n_sp_focus, n_sp_curr)
    n_sp_next <- nrow(proj_mat_next)
    n_sp_pred <- c(n_sp_pred, n_sp_next)
    
    # setting simulation parameters ------------------------------
    # initial conditions (do not affect results)
    N <- rep(1, n_sp_curr)
    # defining growth rates
    r <- rep(1, n_sp_curr)
    # time steps
    times <- seq(0, 100, 0.01)
    # r or k formalism
    formalism <- "K"
    # saving names
    mat_name_curr <- c(mat_name_curr, gsub("\\..*", "", all_files_sub[mat_id]))
    mat_name_next <- c(mat_name_next, gsub("\\..*", "", all_files_sub[mat_id + 1]))
    # rescale matrix
    comp_mat_curr <- scale_matrix(proj_mat_curr, remove_zero = TRUE, 
                                  scale = "col_sum", p = 1,
                                  tol = 10^-8, diag = "one")
    # check if matrix is positive definite
    if(!is.positive.definite(comp_mat_curr + t(comp_mat_curr)))
      stop("Matrix is not positive definite")
    
    # performing simulations ------------------------------
    # species IDs
    sp_curr <- rownames(comp_mat_curr)
    df <- data.frame(sp = sp_curr, prediction = rep(NA, n_sp_curr))
    # carrying capacity is proportional to number of resources consumed by herbivore
    degree <- apply(mat_curr, 2, sum)
    K <- degree / sum(degree)
    # run simulations and extract surviving species
    surv_sp <- lv_pruning(N = N, r = r, K = K, A = -comp_mat_curr, times = times, 
                          extinct_tol = 0.02, formalism = formalism)

    # making predictions ------------------------------
    # predictions
    df$prediction <- "extinct"
    df$prediction[surv_sp] <- "persist"
    sp_persist <- as.character(df$sp[df$prediction == "persist"])
    sp_extinct <- as.character(df$sp[df$prediction == "extinct"])
    # observations
    sp_next <- rownames(proj_mat_next)
    # number of true positives
    tp <- sum(!is.na(match(sp_persist, sp_next)))
    # number of false positives
    fp <- sum(is.na(match(sp_persist, sp_next)))
    # number of true negatives
    tn <- sum(is.na(match(sp_extinct, sp_next)))
    # number of false negatives
    fn <- sum(!is.na(match(sp_extinct, sp_next)))
    # true positive rate
    tpr <- c(tpr, tp / (tp + fn))
    fpr <- c(fpr, fp / (fp + tn))
  }
}

# result data frame ------------------------------
results_df <- data.frame(focus_year = mat_name_curr, n_sp_focus_year = n_sp_focus,
                         prediction_year = mat_name_next, n_sp_pred_year = n_sp_pred,
                         tpr = tpr, fpr = fpr)
results_df$predictor <- rep("good", nrow(results_df))
results_df$predictor[results_df$tpr < results_df$fpr] <- "bad"
results_df$predictor[results_df$tpr == results_df$fpr] <- "neutral"
sum(results_df$predictor == "good") / nrow(results_df)

# create and save plot ------------------------------
col_pal <- c("#E41A1C", "#377EB8", "#000000")
fig <- ggplot(data = results_df, mapping = aes(x = fpr, y = tpr, color = predictor)) +
  geom_point(size = 3.5, alpha = 0.8) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", size = 1.2) +
  scale_color_manual(values = col_pal) +
  xlab("False positive rate") +
  ylab("True positive rate") +
  scale_x_continuous(limits = c(0, 1)) + 
  scale_y_continuous(limits = c(0, 1)) + 
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1),
        axis.text.y = element_text(size = 16),
        axis.title = element_text(size = 20),
        axis.text.x = element_text(size = 16),
        legend.position = "none")
ggsave("figs/fig_SI_predicting_sp_next_year.pdf", fig, width = 14, height = 13, units = "cm")
