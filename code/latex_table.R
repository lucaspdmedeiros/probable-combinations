# Create and save a latex table from the data frame with empirical analyzes results

# cleaning wd, loading functions and packages ------------------------------
rm(list = ls(all = TRUE))
if(!require(xtable)) {install.packages("xtable"); library(xtable)}

# create latex table ------------------------------
results_df <- read.csv(file = "results/fig3/results_2007-2017_year_stage_metanetwork_col_sum.csv",
                       as.is = TRUE)
results_df$stage[results_df$stage == "initial"] <- "Initial"
results_df$stage[results_df$stage == "middle"] <- "Middle"
results_df$stage[results_df$stage == "late"] <- "Late"
results_df$stage <- factor(results_df$stage, 
                           levels = c("Initial", "Middle", "Late"))
results_df$year <- factor(results_df$year,
                          levels = c("2007", "2008", "2009", "2010", 
                                     "2011", "2012", "2013", "2014", "2016", "2017"))
results_df$plot <- factor(results_df$plot,
                          levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9"))
# reorder data frame
results_df <- results_df[order(results_df$stage, results_df$plot, results_df$year), ]
row.names(results_df) <- 1:88
results_df <- results_df[ , c("stage", "plot", "year", "n_sp_real", "n_sp_pot", "omega_pool", 
                              "omega_real", "omega_percentile")]
# creating latex table
tex_table = xtable(results_df, digits = c(0, 0, 0, 0, 0, 0, 3, 3, 3))
# saving file
print.xtable(tex_table, type = "latex", file = "results/Table_S1.tex")
