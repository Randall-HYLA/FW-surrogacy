library(vegan)

setwd("")

nitro_levels <- readRDS("../nitro_levels_matrix.rds")

#Load matrix of FW species
threatened_fw_pa_matrix50 <- readRDS("data/fw_species_t_matrix_zeros.rds")

#CREATE SURROGATE CURVE ####
fw_nitro_matrix <- left_join(threatened_fw_pa_matrix50, nitro_levels, by = c("x","y"))
head(fw_nitro_matrix)

fw_nitro_matrix <- fw_nitro_matrix %>%
  dplyr::select(r_nitrogen_moll_50, everything())

fw_nitro_matrix <- fw_nitro_matrix[complete.cases(fw_nitro_matrix$r_nitrogen_moll_50), ]

# Find column indices where the column sums are not equal to 0
non_zero_col_indices <- which(colSums(fw_nitro_matrix != 0) > 0)
# Subset the dataframe to include only non-zero sum columns
fw_nitro_matrix <- fw_nitro_matrix[, non_zero_col_indices]

####ORGANIZE FW_NITRO_MATRIX FROM HIGH TO LOW VALUES

fw_nitro_matrix <- fw_nitro_matrix %>% 
  arrange(desc((r_nitrogen_moll_50))) #organize from high to low values

fw_nitro_matrix$r_nitrogen_moll_50 <- as.factor(fw_nitro_matrix$r_nitrogen_moll_50 )

View(fw_nitro_matrix)

num_repeats <- 5

#Create five df randomized 

set.seed(123) # For reproducibility
shuffled_list <- purrr::map(1:num_repeats, function(i) {
  fw_nitro_matrix %>%
    dplyr::group_by(r_nitrogen_moll_50) %>%
    dplyr::arrange(r_nitrogen_moll_50, -sample(n())) %>%
    ungroup() %>%
    dplyr::select(-r_nitrogen_moll_50, -x, -y) 
})

View(shuffled_list[[1]])
View(shuffled_list[[2]])

num_samples <- length(shuffled_list)
accu_surrogacy_list <- list()
n=0

for (i in 1:num_samples) {
  acc_curve_temp <- specaccum(shuffled_list[[i]], "collector") #Method "collector" adds sites in the order they happen to be in the data
  accu_surrogacy_list[[i]] <- acc_curve_temp
  n=n+1
  print(n)
}