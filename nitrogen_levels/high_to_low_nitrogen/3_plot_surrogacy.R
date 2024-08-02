library(raster);library(dplyr);library(zonator);library(pracma);library(Rmisc);library(ggplot2)

#LOAD SURROGACY ACCU CURVE LIST AND MAKE DF
plot(accu_surrogacy_list[[5]], ci.type="poly", col="blue", lwd=2, ci.lty=0, ci.col="lightblue")

#Extract values from accumulation curves
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
str(accu_surrogacy_df)

ggplot(data= accu_surrogacy_df, aes(x=sites,y=mean_surrogacy)) +
  geom_point()

accu_surrogacy_df <- accu_surrogacy_df %>%
  dplyr::select(-sites, -mean_surrogacy)

rm(accu_surrogacy_list)
head(accu_surrogacy_df)

# Build sai 
source("../functions_modified.R")

# Generate summary outputs #####

#CAZ
accu_surrogacy_df <- accu_surrogacy_df[1:52309, ]
curve_CAZ_fw_species_optimal <- curve_CAZ_fw_species_optimal[1:52309, ]
curve_CAZ_fw_species_random <- curve_CAZ_fw_species_random[1:52309, ]

## Generate a list to create the surrogate
CAZ_fw_species_stress_sai_curves <- combine_sai_curves(
  target = curve_CAZ_fw_species_optimal,
  surrogate = accu_surrogacy_df,
  random = curve_CAZ_fw_species_random) #accu_random_df

CAZ_fw_comb_threat_sai <- calculate_sai(CAZ_fw_species_stress_sai_curves)[[1]]
CAZ_fw_comb_threat_sai_plot <- plot_sai_curves(CAZ_fw_species_stress_sai_curves)

#ABF
accu_surrogacy_df <- accu_surrogacy_df[1:52309, ]
curve_ABF_fw_species_optimal <- curve_ABF_fw_species_optimal[1:52309, ]
curve_ABF_fw_species_random <- curve_ABF_fw_species_random[1:52309, ]

## Generate a list to create the surrogate
ABF_fw_species_stress_sai_curves <- combine_sai_curves(
  target = curve_ABF_fw_species_optimal,
  surrogate = accu_surrogacy_df,
  random = curve_ABF_fw_species_random) #accu_random_df

ABF_fw_comb_threat_sai <- calculate_sai(ABF_fw_species_stress_sai_curves)[[1]]
ABF_fw_comb_threat_sai_plot <- plot_sai_curves(ABF_fw_species_stress_sai_curves)
