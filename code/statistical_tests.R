# Computes statistics and performs the binomial tests using the theoretical 
# and empirical results

# cleaning wd, loading functions and packages ------------------------------
rm(list = ls(all = TRUE))
if(!require(dplyr)) {install.packages("dplyr"); library(dplyr)}

# theoretical data ------------------------------
# define pool
pool <- "niche"
# read file for regional pool
results_pool_df <- read.csv(paste("results/fig2/", pool, 
                                  "/results_species_pool_theoretical_matrices_n20_", pool, "pool_ncomm500.csv",
                                  sep = ""))
# computing omega for species pool
tapply(results_pool_df$omega_pool, results_pool_df$param, mean)
tapply(results_pool_df$omega_pool, results_pool_df$param, sd)
# read files for realized combinations
path <- paste("results/fig2/", pool, sep = "")
all_files <- dir(path)
all_files_sub <- all_files[grep("lv_pruning", all_files)]
file_paths <- paste(path, all_files_sub, sep = "/")
results_list <- lapply(file_paths, read.csv, header = TRUE)
names(results_list) <- gsub('.{4}$', '', sub("_.*", "", sub(".*param", "", all_files_sub)))
results_df <- bind_rows(results_list)
results_df$param <- as.numeric(results_df$param)
# computing fraction of surviving species
tapply(results_df$n_sp_pruned / results_df$n_sp_pool, results_df$param, mean)
tapply(results_df$n_sp_pruned / results_df$n_sp_pool, results_df$param, sd)
# computing percentile ranks above 0.5
uniq_param <- unique(results_df$param)
mean(results_df$omega_percentile > 0.5)
tapply(results_df$omega_percentile > 0.5, results_df$param, mean)
# performing binomial tests
binom.test(x = sum(results_df$omega_percentile[results_df$param == 0.025] > 0.5),
           n = length(results_df$omega_percentile[results_df$param == 0.025]),
           alternative = "greater")
binom.test(x = sum(results_df$omega_percentile[results_df$param == 0.05] > 0.5),
           n = length(results_df$omega_percentile[results_df$param == 0.05]),
           alternative = "greater")
binom.test(x = sum(results_df$omega_percentile[results_df$param == 0.1] > 0.5),
           n = length(results_df$omega_percentile[results_df$param == 0.1]),
           alternative = "greater")

# empirical data ------------------------------
# define pool
emp_pool <- "year_stage"
# read files 
results_df <- read.csv(file = paste("results/fig3/results_2007-2017_", emp_pool, "_metanetwork_col_sum.csv", 
                                    sep = ""))
results_df$stage <- factor(results_df$stage, 
                           levels = c("initial", "middle", "late"))
results_df$year <- factor(results_df$year,
                          levels = c("2007", "2008", "2009", "2010", 
                                     "2011", "2012", "2013", "2014", "2016", "2017"))
results_df$plot <- factor(results_df$plot,
                          levels = c("1", "2", "3", "4", "5", "6",
                                     "7", "8", "9"))
# computing omega for species pool
tapply(results_df$omega_pool, results_df$stage, mean)
tapply(results_df$omega_pool, results_df$stage, sd)
# compute the analysis of variance
summary(aov(results_df$omega_pool ~ results_df$stage))
# computing fraction of observed species
tapply(results_df$n_sp_real / results_df$n_sp_pot, results_df$stage, mean)
tapply(results_df$n_sp_real / results_df$n_sp_pot, results_df$stage, sd)
# compute the analysis of variance
summary(aov(results_df$n_sp_real / results_df$n_sp_pot ~ results_df$stage))
# computing percentile ranks above 0.5
uniq_param <- unique(results_df$stage)
mean(results_df$omega_percentile > 0.5)
tapply(results_df$omega_percentile > 0.5, results_df$stage, mean)
# performing binomial tests
binom.test(x = sum(results_df$omega_percentile[results_df$stage == "initial"] > 0.5),
           n = length(results_df$omega_percentile[results_df$stage == "initial"]),
           alternative = "greater")
binom.test(x = sum(results_df$omega_percentile[results_df$stage == "middle"] > 0.5),
           n = length(results_df$omega_percentile[results_df$stage == "middle"]),
           alternative = "greater")
binom.test(x = sum(results_df$omega_percentile[results_df$stage == "late"] > 0.5),
           n = length(results_df$omega_percentile[results_df$stage == "late"]),
           alternative = "greater")
