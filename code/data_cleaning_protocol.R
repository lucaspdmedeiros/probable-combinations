# Cleans original field data before performing analyzes;
# check data/data_cleaning_protocol.txt for details

# clean wd
rm(list = ls(all = TRUE))
# define year
year <- "2007"

# read data frame
data <- read.csv(paste("data/", year, "/", year, "_cleaned_full.csv", sep = ""), as.is = TRUE)
data[1:5, 1:15]

# remove columns filled with NAs
cols_na <- apply(data, 2, function(x) all(is.na(x)))
data <- data[ , !cols_na]
#data <- data[ , 1:(ncol(data) - 1)]
# remove rows that contain NAs
data <- na.omit(data)
# check if data frame has any NA left
sum(is.na(data))

# removing dots from plant names
data$PLANT <- gsub("..", " ", data$PLANT, fixed = TRUE)
data$PLANT <- gsub(".", " ", data$PLANT, fixed = TRUE)
# changing double spaces to single spaces
data$PLANT <- gsub("  ", " ", data$PLANT, fixed = TRUE)
# removing trailing and leading spaces
data$PLANT <- trimws(data$PLANT)

# building single networks by merging monthly data for each plot/successional stage
data$PLOT <- as.factor(data$PLOT)
data$SUCCESSION <- as.factor(data$SUCCESSION)
mat_list <- split(data, list(data$PLOT, data$SUCCESSION),
                  drop = TRUE)

# changing species names and merging duplicates
for (i in 1:length(mat_list)) {
  
  # extracting matrix i
  mat <- mat_list[[i]]
  
  # removing species with no name
  if (length(which(mat$PLANT == "")) > 0)
    mat <- mat[- which(mat$PLANT == ""), ]
  if (length(which(mat$PLANT == "(en blanco)")) > 0)
    mat <- mat[- which(mat$PLANT == "(en blanco)"), ]
  if (length(which(mat$PLANT == "NO SAMPLING")) > 0)
    mat <- mat[- which(mat$PLANT == "NO SAMPLING"), ]
  
  # changing wrong plant species names
  mat$PLANT[mat$PLANT == "Acacia macrocantha"] <- "Acacia macracantha"
  mat$PLANT[mat$PLANT == "Buchonsia palmeri"] <- "Bunchosia palmeri"
  mat$PLANT[mat$PLANT == "Caesapinia sp"] <- "Caesalpinia sp"
  mat$PLANT[mat$PLANT == "Coccoloba liebmanni"] <- "Coccoloba liebmannii"
  mat$PLANT[mat$PLANT == "Cordia sp ?"] <- "Cordia sp"
  mat$PLANT[mat$PLANT == "Diospyros aequivis"] <- "Diospyros aequoris"
  mat$PLANT[mat$PLANT == "Lonchocarpus sp ?"] <- "Lonchocarpus sp"
  mat$PLANT[mat$PLANT == "Achatocarpus gracilis H Walt"] <- "Achatocarpus gracilis"
  mat$PLANT[mat$PLANT == "Lonchocarpus (4)"] <- "Lonchocarpus sp 4"
  mat$PLANT[mat$PLANT == "Hpp Hemiangium excelsum"] <- "Hemiangium excelsum"
  mat$PLANT[mat$PLANT == "Hpp Hemiangium excelsium"] <- "Hemiangium excelsum"
  mat$PLANT[mat$PLANT == "Serjania brachycarpa A Gray"] <- "Serjania brachycarpa"
  mat$PLANT[mat$PLANT == "aff Calliandra emarginata"] <- "Calliandra emarginata"
  mat$PLANT[mat$PLANT == "aff Lonchocarpus sp L"] <- "Lonchocarpus sp L"
  mat$PLANT[mat$PLANT == "Celtis iguanaea (Jacq ) Sarg"] <- "Celtis iguanaea"
  mat$PLANT[mat$PLANT == "Randia thurberi S Wats"] <- "Randia thurberi"
  mat$PLANT[mat$PLANT == "Cordia aff gerascanthus"] <- "Cordia gerascanthus"
  mat$PLANT[mat$PLANT == "Achatocarpusgracilis"] <- "Achatocarpus gracilis"
  mat$PLANT[mat$PLANT == "Bourreria purpussii"] <- "Bourreria purpusii"
  mat$PLANT[mat$PLANT == "Cordia aff Gerascanthus"] <- "Cordia gerascanthus"
  mat$PLANT[mat$PLANT == "Cordia eleagnoides"] <- "Cordia elaeagnoides"
  mat$PLANT[mat$PLANT == "Cordia sp1"] <- "Cordia sp 1"
  mat$PLANT[mat$PLANT == "Dyphisa occidentalis"] <- "Diphysa occidentalis"
  mat$PLANT[mat$PLANT == "Forchammeria pallida"] <- "Forchhammeria pallida"
  mat$PLANT[mat$PLANT == "Forchammeria sessiliflora"] <- "Forchhammeria sessiliflora"
  mat$PLANT[mat$PLANT == "Forchhammeria sesiiflora"] <- "Forchhammeria sessiliflora"
  mat$PLANT[mat$PLANT == "Gettarda elliptica"] <- "Guettarda elliptica"
  mat$PLANT[mat$PLANT == "Hemiangium excelsium"] <- "Hemiangium excelsum"
  mat$PLANT[mat$PLANT == "Ipomea wolcottiana"] <- "Ipomoea wolcottiana"
  mat$PLANT[mat$PLANT == "Malphigia emilae"] <- "Malpighia emiliae" 
  mat$PLANT[mat$PLANT == "Malphigia emiliae"] <- "Malpighia emiliae"
  mat$PLANT[mat$PLANT == "Myrospermun frutensces"] <- "Myrospermum frutescens"
  mat$PLANT[mat$PLANT == "Myrospermun frutescens"] <- "Myrospermum frutescens"
  mat$PLANT[mat$PLANT == "Phyllanthus mociniatus"] <- "Phyllanthus mocinianus"
  mat$PLANT[mat$PLANT == "Phyllanthus mociniaus"] <- "Phyllanthus mocinianus"
  mat$PLANT[mat$PLANT == "Stemmadenia donell-smithii"] <- "Stemmadenia donnell-smithii"
  mat$PLANT[mat$PLANT == "Thounia paucidentata"] <- "Thouinia paucidentata"
  mat$PLANT[mat$PLANT == "Trichilia triflora"] <- "Trichilia trifolia"
  mat$PLANT[mat$PLANT == "Zanthoxyllum sp 2"] <- "Zanthoxylum sp 2"
  mat$PLANT[mat$PLANT == "Zanthoxyllum sp 3"] <- "Zanthoxylum sp 3"
  mat$PLANT[mat$PLANT == "Zanthoxyllum caribaeum"] <- "Zanthoxylum caribaeum"
  
  # changing wrong herbivore species names
  colnames(mat)[colnames(mat) == "O114....ANTES.O460."] <- "O114"
  colnames(mat)[colnames(mat) == "O460"] <- "O114"
  colnames(mat)[colnames(mat) == "O16....ANTES.O478."] <- "O16"
  colnames(mat)[colnames(mat) == "O478"] <- "O16"
  colnames(mat)[colnames(mat) == "O20....ANTES.O189."] <- "O20"
  colnames(mat)[colnames(mat) == "O189"] <- "O20"
  colnames(mat)[colnames(mat) == "O216....ANTES.O58."] <- "O216"
  colnames(mat)[colnames(mat) == "O58"] <- "O216"
  colnames(mat)[colnames(mat) == "O23....ANTES.O161."] <- "O23"
  colnames(mat)[colnames(mat) == "O161"] <- "O23"
  colnames(mat)[colnames(mat) == "O26....ANTES.O608."] <- "O26"
  colnames(mat)[colnames(mat) == "O608"] <- "O26"
  colnames(mat)[colnames(mat) == "O29....ANTES.O6."] <- "O29"
  colnames(mat)[colnames(mat) == "O6"] <- "O29"
  colnames(mat)[colnames(mat) == "O32....ANTES.O498."] <- "O32"
  colnames(mat)[colnames(mat) == "O498"] <- "O32"
  colnames(mat)[colnames(mat) == "O326....ANTES.O284."] <- "O326"
  colnames(mat)[colnames(mat) == "O284"] <- "O326"
  colnames(mat)[colnames(mat) == "O440....ANTES.O565.Y.O569."] <- "O440"
  colnames(mat)[colnames(mat) == "O565.Y.O569"] <- "O440"
  colnames(mat)[colnames(mat) == "O46....ANTES.O481."] <- "O46" 
  colnames(mat)[colnames(mat) == "O481"] <- "O46"
  colnames(mat)[colnames(mat) == "O75....ANTES.O131."] <- "O75"
  colnames(mat)[colnames(mat) == "O131"] <- "O75"
  colnames(mat)[colnames(mat) == "O79....ANTES.O356."] <- "O79"
  colnames(mat)[colnames(mat) == "O356"] <- "O79"
  colnames(mat)[colnames(mat) == "O9....ANTES.O2."] <- "O9"
  colnames(mat)[colnames(mat) == "O2"] <- "O9"
  colnames(mat)[colnames(mat) == "O99....ANTES.O363."] <- "O99"
  colnames(mat)[colnames(mat) == "O363"] <- "O99"
  colnames(mat)[colnames(mat) == "O525.1"] <- "O525"
  colnames(mat)[colnames(mat) == "O.216"] <- "O216"
  colnames(mat)[colnames(mat) == "O330."] <- "O330"
  colnames(mat)[colnames(mat) == "Hylesia"] <- "O299"
  colnames(mat)[colnames(mat) == "O20....O256."] <- "O20"
  colnames(mat)[colnames(mat) == "O256"] <- "O20"
  colnames(mat)[colnames(mat) == "O66.O337"] <- "O66"
  
  # merging duplicates (rows containing data from the same plant species)
  row_names <- as.character(unique(mat$PLANT))
  col_names <- as.character(colnames(mat)[5:ncol(mat)])
  row_mat <- matrix(NA, nrow = length(row_names), ncol = length(col_names))
  rownames(row_mat) <- row_names
  colnames(row_mat) <- col_names
  for (j in 1:length(row_names)) {
    row_vec <- subset(mat, PLANT == row_names[j])
    if (length(dim(row_vec)) > 1) {
      row_vec <- row_vec[ , 5:length(row_vec)]
      row_vec_sum <- apply(row_vec, 2, sum)
      row_vec_sum[row_vec_sum > 0] <- 1
    } else {
      row_vec <- row_vec[5:length(row_vec)]
      row_vec_sum <- row_vec
    }
    row_mat[row_names[j], ] <- row_vec_sum
  }
  # merging duplicates (columns containing data from the same herbivore species)
  uniq_col_names <- unique(col_names)
  col_mat <- matrix(NA, nrow = length(row_names), ncol = length(uniq_col_names))
  rownames(col_mat) <- row_names
  colnames(col_mat) <- uniq_col_names
  for (j in 1:length(uniq_col_names)) {
    col_vec <- row_mat[ , which(colnames(row_mat) == uniq_col_names[j])]
    if (length(dim(col_vec)) > 1) {
      col_vec_sum <- apply(col_vec, 1, sum)
      col_vec_sum[col_vec_sum > 0] <- 1
    } else {
      col_vec_sum <- col_vec
    }
    col_mat[ , uniq_col_names[j]] <- col_vec_sum
  }
  mat <- col_mat
  
  # ordering species names
  mat <- mat[order(rownames(mat)), order(colnames(mat))]
  # adding results to matrix list
  mat_list[[i]] <- mat
}

# saving interaction matrices 
mat_names <- gsub(".", "_", names(mat_list), fixed = TRUE)
mat_names <- gsub("EARLY", "initial", mat_names, fixed = TRUE)
mat_names <- gsub("LATE", "middle", mat_names, fixed = TRUE)
mat_names <- gsub("FOREST", "late", mat_names, fixed = TRUE)
mat_id <- grep("INITIAL", mat_names, invert = TRUE)

for (i in 1:length(mat_id)) {
  mat_name <- mat_names[mat_id[i]]
  write.table(mat_list[[mat_id[i]]], 
              file = paste("data/", year, "/", mat_name, "_", year, ".txt", sep = ""),
              sep = "\t", col.names = NA)
}
