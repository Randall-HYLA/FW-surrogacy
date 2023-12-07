# This script creates the target curves of ABF and CAZ, that would be use for the analysis in script Low_to_High_stress.R

# The final files of this script that will be used in the other script are available in the folder data.

library(raster);library(dplyr);library(zonator)

setwd("") # Set proper working directory

# Load water stress matrix
water_stress <- readRDS("data/water_stress_matrix.rds")

# Load matrix of FW species
threatened_fw_pa_matrix50 <- readRDS("data/fw_species_t_matrix_zeros.rds")

fw_stress_matrix <- left_join(threatened_fw_pa_matrix50, 
                              water_stress, by = c("x","y"))

fw_stress_matrix <- na.omit(fw_stress_matrix)

fw_stress_matrix <- fw_stress_matrix %>% 
  dplyr::arrange((fw_stress_matrix)) # Organize from low to high values

fw_stress_matrix <- fw_stress_matrix[, -ncol(fw_stress_matrix)]

# Find column indices where the column sums are not equal to 0
non_zero_col_indices <- which(colSums(fw_stress_matrix != 0) > 0)

# Subset the dataframe to include only non-zero sum columns
fw_stress_matrix <- fw_stress_matrix[, non_zero_col_indices] #3606 spps

# Generate rasters from cell-by-species matrices

# FW SPECIES (TARGET)  
head(fw_stress_matrix)
names(fw_stress_matrix) <- gsub(x = names(fw_stress_matrix), 
                                         pattern = "//.", replacement = "_")

purrr::map(names(fw_stress_matrix)[-c(1:2)], function(sp){
  r <- rasterFromXYZ(fw_stress_matrix %>% dplyr::select(x, y, sp))
  proj4string(r) <- "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs"
  writeRaster(r, filename = paste0("data/input_rasters/50/fw_species/", sp, ".tif"), format = "GTiff")
})

# Run Zonation 

# Install R zonator package, if not already installed
if (!("zonator" %in% installed.packages()[,"Package"])) install.packages("zonator")
# Create .dat template files for running core area Zonation (CAZ) or additive benefit function (ABF)
# Create template files by copying "template.dat" file included with zonator installation

zonator_path <- find.package("zonator")
setwd(zonator_path)
file.copy("extdata/template.dat", "extdata/template_CAZ.dat")
file.copy("extdata/template.dat", "extdata/template_ABF.dat")
# Update removal rule (rule 2 is the Additive Benefit Function) within "template_ABF.dat" file
template_ABF <- readLines("extdata/template_ABF.dat")
template_ABF[2] <- "removal rule = 2"
writeLines(template_ABF, "extdata/template_ABF.dat")

# Create directory to store outputs
library(zonator)
getwd()
setwd("C:/FW_analysis/abiotic")
dir.create("output")
dir.create("output/zonation_runs")
dir.create("output/zonation_runs/50")

source("functions.R")

# Run zonation

run_zonation(id = "fw_species", runs = 5, algorithm = "CAZ", out_dir = "output/zonation_runs/50", input_rasters_dir = "data/input_rasters/50")

run_zonation(id = "fw_species", runs = 5, algorithm = "ABF", out_dir = "output/zonation_runs/50", input_rasters_dir = "data/input_rasters/50")

# LOAD RANKS

# CAZ algorithm
ranks_CAZ_fw_species <- extract_zonation_ranks(taxon = "fw_species", algorithm = "CAZ", out_dir = "output/zonation_runs/50")

# ABF algorithm
ranks_ABF_fw_species <- extract_zonation_ranks(taxon = "fw_species", algorithm = "ABF", out_dir = "output/zonation_runs/50")

# Generate Surrogacy Curves
source("functions_modified.R")

# CAZ

# FW species as target

# Optimal curve
curve_CAZ_fw_species_optimal <- get_sai_curves_target(
  target_matrix = fw_stress_matrix, 
  target_ranks = ranks_CAZ_fw_species, 
  algorithm = "CAZ", 
  runs = 5)

saveRDS(curve_CAZ_fw_species_optimal, "output/zonation_runs/50/curve_CAZ_fw_species_optimal.rds")

curve_CAZ_fw_species_optimal <- readRDS("output/zonation_runs/50/curve_CAZ_fw_species_optimal.rds")

# ABF

# FW species as target

# Optimal curve
 
curve_ABF_fw_species <- get_sai_curves_target(
  target_matrix = fw_stress_matrix, 
  target_ranks = ranks_ABF_fw_species, 
  algorithm = "ABF", 
  runs = 5)

saveRDS(curve_ABF_fw_species, "output/zonation_runs/50/curve_ABF_fw_species_optimal.rds")

curve_ABF_fw_species_optimal <- readRDS("output/zonation_runs/50/curve_ABF_fw_species_optimal.rds")

# Files have been created that will be used in the R script 2_Low_to_High_stress.R to calculate the surrogacy.