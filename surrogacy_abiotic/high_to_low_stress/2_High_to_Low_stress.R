# Surrogacy analysis of FW_species as targets and water stress (high to low stress) as surrogates

library(sf); library(dplyr); library(foreach); library(forcats); 
library(raster); library(data.table); library(ggplot2); library(vegan);
library(pracma); library(Rmisc)
  
setwd("") # Set proper working directory

# Load water stress matrix
water_stress <- readRDS("data/water_stress_matrix.rds")

# Load matrix of FW species
threatened_fw_pa_matrix50 <- readRDS("data/fw_species_t_matrix_zeros.rds")

# CREATE SURROGATE CURVE 
fw_stress_matrix <- left_join(threatened_fw_pa_matrix50, water_stress, by = c("x","y"))

water_stress$water_stress_raster_res50

fw_stress_matrix <- fw_stress_matrix %>%
  dplyr::select(water_stress_raster_res50, everything())

fw_stress_matrix <- fw_stress_matrix[complete.cases(fw_stress_matrix$water_stress_raster_res50), ]

# Find column indices where the column sums are not equal to 0
non_zero_col_indices <- which(colSums(fw_stress_matrix != 0) > 0)
# Subset the dataframe to include only non-zero sum columns
fw_stress_matrix <- fw_stress_matrix[, non_zero_col_indices]

# Calculate Row Sums
row_sums <- rowSums(fw_stress_matrix[, 4:3607])

# Identify Rows that sum up to 0
rows_with_0 <- which(row_sums == 0)

# Rows with no Species
matrix_with_zero_rows <- fw_stress_matrix[rows_with_0, ]

# Calculate the row sums
row_sums <- rowSums(fw_stress_matrix[, 4:3607])
new_matrix <- fw_stress_matrix[, 1:3]

# Add a new column with the row sums to the DataFrame
new_matrix$RowSums <- row_sums
summary(new_matrix)

new_matrix <- new_matrix %>%
  mutate(Category = ifelse(RowSums > 0, ">=1", "0"))

(p1 <- ggplot(new_matrix, aes(x = x, y = y, fill = Category)) +
  geom_tile() +
  coord_equal() +
  scale_fill_manual(values = c(">=1" = "blue", "0" = "gray")) +  # Customize fill colors
  labs(x = "X Axis", y = "Y Axis", fill = "Value") +  # Add axis labels and legend title
  theme_minimal())  # Customize the plot's theme

# Organize water stress dataset from low to high values
fw_stress_matrix <- fw_stress_matrix %>% 
  arrange(desc(water_stress_raster_res50)) #Arrange from high to low

fw_stress_matrix$water_stress_raster_res50 <- as.factor(fw_stress_matrix$water_stress_raster_res50)

# Create five randomized datasets keeping the low to high values order 
num_repeats <- 5

# Create five df randomized 
set.seed(123) # For reproducibility

shuffled_list <- purrr::map(1:num_repeats, function(i) {
  fw_stress_matrix %>%
    group_by(water_stress_raster_res50) %>%
    arrange(water_stress_raster_res50, -sample(n())) %>%
    ungroup() %>%
    dplyr::select(-water_stress_raster_res50, -x, -y) 
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

# Have a look at curves:
plot(accu_surrogacy_list[[5]], ci.type="poly", col="blue", lwd=2, ci.lty=0, ci.col="lightblue")

# Extract values from accumulation curves
name_column <- accu_surrogacy_list[[1]]$richness
accu_surrogacy_df <- data.frame(curve1 = name_column)

# Iterate through each object in the list
for (i in 2:5) {
  extracted_column <- accu_surrogacy_list[[i]]$richness
  accu_surrogacy_df <- cbind(accu_surrogacy_df, extracted_column)
  colnames(accu_surrogacy_df)[i] <- paste("curve", i, sep = "")
}  

# Calculate mean per row
accu_surrogacy_df$mean_surrogacy <- rowMeans(accu_surrogacy_df)

accu_surrogacy_df  <- accu_surrogacy_df  %>%
  tibble::rownames_to_column(var = "sites")

accu_surrogacy_df$sites <- as.numeric(accu_surrogacy_df$sites)

ggplot(data= accu_surrogacy_df, aes(x=sites,y=mean_surrogacy)) +
  geom_point()

accu_surrogacy_df <- accu_surrogacy_df %>%
  dplyr::select(-sites, -mean_surrogacy)

rm(accu_surrogacy_list)

# Generate random curves

# Set the species matrix as beginning
fw_stress_matrix <- left_join(threatened_fw_pa_matrix50, water_stress, by = c("x","y"))
water_stress$water_stress_raster_res50

fw_stress_matrix <- fw_stress_matrix %>%
  dplyr::select(water_stress_raster_res50, everything())

fw_stress_matrix <- fw_stress_matrix[complete.cases(fw_stress_matrix$water_stress_raster_res50), ]

# Find column indices where the column sums are not equal to 0
non_zero_col_indices <- which(colSums(fw_stress_matrix != 0) > 0)

# Subset the dataframe to include only non-zero sum columns
fw_stress_matrix <- fw_stress_matrix[, non_zero_col_indices]

# Filter the water stress var and the coordinates
fw_stress_matrix_f <- fw_stress_matrix %>%
  select(-water_stress_raster_res50, -x, -y)

# Generate 10 random datasets to create random curves
num_samples <- 10
sample_size <- 45694

sampled_data <- list()

for (i in 1:num_samples) {
  sampled_indices <- sample(nrow(fw_stress_matrix_f), sample_size)
  sampled_data[[i]] <- fw_stress_matrix_f[sampled_indices, ]
}

View(sampled_data[[3]])

# Create accumulation curves for random curve
num_samples <- length(sampled_data)

# This iteration will take 20 min aprox
accu_random_list <- list()
n=0

for (i in 1:num_samples) {
  acc_curve_temp <- specaccum(sampled_data[[i]], "collector")
  accu_random_list[[i]] <- acc_curve_temp
  n=n+1
  print(n)
}

# Extract values from accumulation curves
name_column <- accu_random_list[[1]]$richness
accu_random_df <- data.frame(curve1 = name_column)

# Iterate through each object in the list
for (i in 2:10) {
  extracted_column <- accu_random_list[[i]]$richness
  accu_random_df <- cbind(accu_random_df, extracted_column)
  colnames(accu_random_df)[i] <- paste("curve", i, sep = "")
 }  

# Calculate the number of rows in the data
num_rows <- nrow(accu_random_df)
# Change row names to a sequence from 1 to num_rows
rownames(accu_random_df) <- 1:num_rows

# Calculate mean per row
accu_random_df$mean_random <- rowMeans(accu_random_df)

accu_random_df  <- accu_random_df  %>%
  tibble::rownames_to_column(var = "sites")

accu_random_df$sites <- as.numeric(accu_random_df$sites)
str(accu_random_df)

ggplot(data= accu_random_df, aes(x=sites,y=mean_random)) +
  geom_point()

accu_random_df <- accu_random_df %>%
  dplyr::select(-sites, -mean_random)


# At this point we can start calculating the surrogacy, we will used the target curves generate in 1_zonation_target_curves.R

# Load curves
curve_CAZ_fw_species_optimal <- readRDS("data/curve_CAZ_fw_species_optimal.rds")
curve_ABF_fw_species_optimal <- readRDS("data/curve_ABF_fw_species_optimal.rds")

# Build sai 
source("functions_modified.R")

# Generate summary outputs

# CAZ
accu_surrogacy_df <- accu_surrogacy_df[1:45682, ]
curve_CAZ_fw_species_optimal <- curve_CAZ_fw_species_optimal[1:45682, ]
curve_CAZ_fw_species_random <- curve_CAZ_fw_species_random[1:45682, ]

# Generate a list to create the surrogate
CAZ_fw_species_stress_sai_curves <- combine_sai_curves(
  target = curve_CAZ_fw_species_optimal,
  surrogate = accu_surrogacy_df,
  random = accu_random_df)

CAZ_fw_comb_threat_sai <- calculate_sai(CAZ_fw_species_stress_sai_curves)[[1]]
CAZ_fw_comb_threat_sai_plot <- plot_sai_curves(CAZ_fw_species_stress_sai_curves)

ggsave(CAZ_fw_comb_threat_sai_plot,
       filename = "CAZ_fw_comb_threat_sai_plot.tiff",
       device = "tiff", path = "./",
       width = 95, height = 95, units="mm",dpi = 300)  

# ABF
accu_surrogacy_df <- accu_surrogacy_df[1:45682, ]
curve_ABF_fw_species_optimal <- curve_ABF_fw_species_optimal[1:45682, ]
curve_ABF_fw_species_random <- curve_ABF_fw_species_random[1:45682, ]

# Generate a list to create the surrogate
ABF_fw_species_stress_sai_curves <- combine_sai_curves(
  target = curve_ABF_fw_species_optimal,
  surrogate = accu_surrogacy_df,
  random = accu_random_df) 

ABF_fw_comb_threat_sai <- calculate_sai(ABF_fw_species_stress_sai_curves)[[1]]
ABF_fw_comb_threat_sai_plot <- plot_sai_curves(ABF_fw_species_stress_sai_curves)

ggsave(ABF_fw_comb_threat_sai_plot,
       filename = "ABF_fw_comb_threat_sai_plot.tiff",
       device = "tiff", path = "./",
       width = 95, height = 95, units="mm",dpi = 300)  
