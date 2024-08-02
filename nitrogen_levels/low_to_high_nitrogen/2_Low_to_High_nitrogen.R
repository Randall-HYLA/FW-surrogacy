# This script creates the target curves of ABF and CAZ, that would be use for the analysis in the script nitrogen.R

#The final files of this script that will be used in the other script are available in the folder data.

library(raster);library(dplyr);library(zonator);library(vegan)

setwd("") #Set proper working directory

#Load nitrogen matrix
nitro_levels <- readRDS("nitro_levels_matrix.rds")

#Load matrix of FW species
threatened_fw_pa_matrix50 <- readRDS("data/fw_species_t_matrix_zeros.rds")

#Generate surrogacy curve
fw_nitro_matrix <- left_join(threatened_fw_pa_matrix50, nitro_levels, by = c("x","y"))
head(fw_nitro_matrix)

fw_nitro_matrix <- fw_nitro_matrix %>%
  dplyr::select(r_nitrogen_moll_50, everything())

fw_nitro_matrix <- fw_nitro_matrix[complete.cases(fw_nitro_matrix$r_nitrogen_moll_50), ]

# Find column indices where the column sums are not equal to 0
non_zero_col_indices <- which(colSums(fw_nitro_matrix != 0) > 0)
# Subset the dataframe to include only non-zero sum columns
fw_nitro_matrix <- fw_nitro_matrix[, non_zero_col_indices]

####ORGANIZE FW_NITRO_MATRIX FROM LOW TO HIGH VALUES

fw_nitro_matrix <- fw_nitro_matrix %>% 
  arrange((r_nitrogen_moll_50)) #organize from low to high values

fw_nitro_matrix$r_nitrogen_moll_50 <- as.factor(fw_nitro_matrix$r_nitrogen_moll_50 )

num_repeats <- 5

#Create five df randomized 
set.seed(123) # For reproducibility

shuffled_list <- purrr::map(1:num_repeats, function(i) {
  fw_nitro_matrix %>%
    group_by(r_nitrogen_moll_50) %>%
    arrange(r_nitrogen_moll_50, -sample(n())) %>%
    ungroup() %>%
    dplyr::select(-r_nitrogen_moll_50, -x, -y) 
})

View(shuffled_list[[1]])
View(shuffled_list[[2]])

# Generate accumulation curves of FW_species based on water stress (low to high) which will be used as surrogate curve

num_samples <- length(shuffled_list)
accu_surrogacy_list <- list()
n=0

for (i in 1:num_samples) {
  acc_curve_temp <- specaccum(shuffled_list[[i]], "collector") #Method "collector" adds sites in the order they happen to be in the data
  accu_surrogacy_list[[i]] <- acc_curve_temp
  n=n+1
  print(n)
}


