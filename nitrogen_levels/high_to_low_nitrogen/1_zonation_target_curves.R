## Set things up
library(raster);library(dplyr);library(zonator)
setwd("C:/surrogacy_world_bank/high_to_low")

#Create new input raster file directories
dir.create("data")
dir.create("data/input_rasters")
#' #### Create directories for 50km resolution rasters
dir.create("data/input_rasters/50")
dir.create("data/input_rasters/50/fw_species")

#Load matrix of FW species
nitro_levels <- readRDS("../nitro_levels_matrix.rds")
threatened_fw_pa_matrix50 <- readRDS("G:/surrogacy_fw_species/Surrogacy_FW_PART2/data/fw_species_t_matrix_zeros.rds")

#CREATE SURROGATE CURVE ####
fw_nitro_matrix <- left_join(threatened_fw_pa_matrix50, nitro_levels, by = c("x","y"))
fw_nitro_matrix <- na.omit(fw_nitro_matrix)
fw_nitro_matrix <- fw_nitro_matrix %>% arrange(desc(r_nitrogen_moll_50)) #organize from high TO LOW values
View(fw_nitro_matrix)
fw_nitro_matrix$r_nitrogen_moll_50 <- NULL

# Find column indices where the column sums are not equal to 0
non_zero_col_indices <- which(colSums(fw_nitro_matrix != 0) > 0)

# Subset the dataframe to include only non-zero sum columns
fw_nitro_matrix <- fw_nitro_matrix[, non_zero_col_indices] #3606 spps

## Generate rasters from cell-by-species matrices

### FW SPECIES (TARGET)  ####
head(fw_nitro_matrix)
names(fw_nitro_matrix) <- gsub(x = names(fw_nitro_matrix), 
                                         pattern = "\\.", replacement = "_")

purrr::map(names(fw_nitro_matrix)[-c(1:2)], function(sp){
  r <- rasterFromXYZ(fw_nitro_matrix %>% dplyr::select(x, y, sp))
  proj4string(r) <- "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs"
  writeRaster(r, filename = paste0("data/input_rasters/50/fw_species/", sp, ".tif"), format = "GTiff")
})

## Run Zonation ####

##### Install R zonator package, if not already installed
if (!("zonator" %in% installed.packages()[,"Package"])) install.packages("zonator")
### Create .dat template files for running core area Zonation (CAZ) or additive benefit function (ABF)
## Create template files by copying "template.dat" file included with zonator installation

zonator_path <- find.package("zonator")
setwd(zonator_path)
file.copy("extdata/template.dat", "extdata/template_CAZ.dat")
file.copy("extdata/template.dat", "extdata/template_ABF.dat")
## Update removal rule (rule 2 is the Additive Benefit Function) within "template_ABF.dat" file
template_ABF <- readLines("extdata/template_ABF.dat")
template_ABF[2] <- "removal rule = 2"
writeLines(template_ABF, "extdata/template_ABF.dat")

#' ### Create directory to store outputs
library(zonator)
getwd()
dir.create("output")
dir.create("output/zonation_runs")
dir.create("output/zonation_runs/50")

source("../functions.R")

#Run zonation
run_zonation(id = "fw_species", runs = 5, algorithm = "CAZ", out_dir = "output/zonation_runs/50", input_rasters_dir = "data/input_rasters/50")

run_zonation(id = "fw_species", runs = 5, algorithm = "ABF", out_dir = "output/zonation_runs/50", input_rasters_dir = "data/input_rasters/50")

#LOAD RANKS

# CAZ algorithm
ranks_CAZ_fw_species <- extract_zonation_ranks(taxon = "fw_species", algorithm = "CAZ", out_dir = "output/zonation_runs/50")

# ABF algorithm
ranks_ABF_fw_species <- extract_zonation_ranks(taxon = "fw_species", algorithm = "ABF", out_dir = "output/zonation_runs/50")

# Generate Surrogacy Curves ####
source("../functions_modified.R")

# CAZ

# FW species as target

# Optimal curve
curve_CAZ_fw_species_optimal <- get_sai_curves_target(
  target_matrix = fw_nitro_matrix, 
  target_ranks = ranks_CAZ_fw_species, 
  algorithm = "CAZ", 
  runs = 5)

# ABF
# FW species as target
# Optimal curve
 
curve_ABF_fw_species <- get_sai_curves_target(
  target_matrix = fw_nitro_matrix, 
  target_ranks = ranks_ABF_fw_species, 
  algorithm = "ABF", 
  runs = 5)

#RANDOM
getwd()
#CAZ
ranks_CAZ_random_fw_species <- random_ranks(ranks_CAZ_fw_species, runs = 100)
#ABF
ranks_ABF_random_fw_species <- random_ranks(ranks_ABF_fw_species, runs = 100)

vert_classes <- c("fw_species")

map_extent <- purrr::map(vert_classes, 
                         extract_zonation_maps, 
                         algorithm = "CAZ", 
                         out_dir = "output/zonation_runs", 
                         suffix = "062023", 
                         resolution = "50") %>% 
  purrr::set_names(vert_classes)

#terrestrial_extent <- map_extent[["combined"]] > 0 | map_extent[["fw_species"]] > 0 
terrestrial_extent <- map_extent[["fw_species"]]

terrestrial_extent <- as.data.frame(terrestrial_extent, xy = TRUE) %>% 
  dplyr::mutate(cellID = 1:length(terrestrial_extent[])) %>% 
  dplyr::filter(complete.cases(.))

ranks_random_terrestrial <- random_ranks(terrestrial_extent$cellID, runs = 100)

## Generate random expectations ####
#'
#' ### FW species as target

#' #### CAZ
#' ##### Terrestrial extent
curve_CAZ_fw_species_random <- 
  get_sai_curves_random(target_matrix = fw_nitro_matrix, 
                        random_ranks = ranks_random_terrestrial, 
                        reference_coordinates = ranks_CAZ_fw_species$reference_coordinates, 
                        algorithm = "CAZ", 
                        runs = 10)

#' #### ABF
#' ##### Terrestrial extent
curve_ABF_fw_species_random <- 
  get_sai_curves_random(target_matrix = fw_nitro_matrix, 
                        random_ranks = ranks_random_terrestrial, 
                        reference_coordinates = ranks_ABF_fw_species$reference_coordinates, 
                        algorithm = "ABF", 
                        runs = 10)

