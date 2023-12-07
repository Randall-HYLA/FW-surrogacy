# ---
# Title: FW Species Assessment Surrogacy in Conservation For Threatened Species - 50km resolution analyses
# ---

# Set things up
library(raster);library(dplyr);library(zonator)

setwd("") #Set proper working directory

# Load custom functions
source("functions_code/functions.R")

# Create new input raster file directories
dir.create("data/input_rasters")

# Create directories for 50km resolution rasters
dir.create("data/input_rasters/50")
dir.create("data/input_rasters/50/fw_species")
dir.create("data/input_rasters/50/reptiles")
dir.create("data/input_rasters/50/birds")
dir.create("data/input_rasters/50/mammals")
dir.create("data/input_rasters/50/amphibians")
dir.create("data/input_rasters/50/combined") # Reptiles, Birds, Mammals and Amphibians

dir.create("data/input_rasters/50/decapods")
dir.create("data/input_rasters/50/fishes")
dir.create("data/input_rasters/50/odonata")

# Load cell-by-species matrices

threatened_reptiles_pa_matrix50 <- readRDS("data/reptile_t_matrix_zeros.rds")
threatened_birds_pa_matrix50 <- readRDS("data/birds_t_matrix_zeros.rds")
threatened_mammals_pa_matrix50 <- readRDS("data/mammals_t_matrix_zeros.rds")
threatened_amphibians_pa_matrix50 <- readRDS("data/amphibians_t_matrix_zeros.rds")

threatened_odonata_pa_matrix50 <- readRDS("data/odonata_t_matrix_zeros.rds")
threatened_decapods_pa_matrix50 <- readRDS("data/decapods_t_matrix_zeros.rds")
threatened_fishes_pa_matrix50 <- readRDS("data/fishes_t_matrix_zeros.rds")

threatened_fw_pa_matrix50 <- readRDS("data/fw_species_t_matrix_zeros.rds") #Combined FW species (FW fishes, odonata, decapods and decapods)

# 1. Generate rasters from cell-by-species matrices ####

# Targets

# FW SPECIES COMBINED (TARGET) (ODONATA + FW FISHES + DECAPODS)
names(threatened_fw_pa_matrix50) <- gsub(x = names(threatened_fw_pa_matrix50), 
                                         pattern = "\\.", replacement = "_")

purrr::map(names(threatened_fw_pa_matrix50)[-c(1:2)], function(sp){
  r <- rasterFromXYZ(threatened_fw_pa_matrix50 %>% dplyr::select(x, y, sp))
  proj4string(r) <- "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs"
  writeRaster(r, filename = paste0("data/input_rasters/50/fw_species/", sp, ".tif"), format = "GTiff")
})

# Check projection
raster("data/input_rasters/50/fw_species/Aborichthys_garoensis.tif")

# ODONATA (TARGET)  
names(threatened_odonata_pa_matrix50) <- gsub(x = names(threatened_odonata_pa_matrix50), 
                                              pattern = "\\.", replacement = "_")

purrr::map(names(threatened_odonata_pa_matrix50)[-c(1:2)], function(sp){
  r <- rasterFromXYZ(threatened_odonata_pa_matrix50 %>% dplyr::select(x, y, sp))
  proj4string(r) <- "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs"
  writeRaster(r, filename = paste0("data/input_rasters/50/odonata/", sp, ".tif"), format = "GTiff")
})

# Check projection
raster("data/input_rasters/50/odonata/Crenigomphus_abyssinicus.tif")

# DECAPODS (TARGET)

names(threatened_decapods_pa_matrix50) <- gsub(x = names(threatened_decapods_pa_matrix50), 
                                               pattern = "\\.", replacement = "_")

purrr::map(names(threatened_decapods_pa_matrix50)[-c(1:2)], function(sp){
  r <- rasterFromXYZ(threatened_decapods_pa_matrix50 %>% dplyr::select(x, y, sp))
  proj4string(r) <- "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs"
  writeRaster(r, filename = paste0("data/input_rasters/50/decapods/", sp, ".tif"), format = "GTiff")
})

# We will also used combined amphibians, birds reptiles and mammals and we will named as combined

# Copy out amphibians, birds reptiles and mammals into a "combined" folder
file.copy(from = c(list.files("data/input_rasters/50/birds", full.names = TRUE),
                   list.files("data/input_rasters/50/mammals", full.names = TRUE),
                   list.files("data/input_rasters/50/reptiles", full.names = TRUE),
                   list.files("data/input_rasters/50/amphibians", full.names = TRUE)), 
          to = "data/input_rasters/50/combined" 
)

# FW FISHES (TARGET)

names(threatened_fishes_pa_matrix50) <- gsub(x = names(threatened_fishes_pa_matrix50), 
                                             pattern = "\\.", replacement = "_")

purrr::map(names(threatened_fishes_pa_matrix50)[-c(1:2)], function(sp){
  r <- rasterFromXYZ(threatened_fishes_pa_matrix50 %>% dplyr::select(x, y, sp))
  proj4string(r) <- "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs"
  writeRaster(r, filename = paste0("data/input_rasters/50/fishes/", sp, ".tif"), format = "GTiff")
})

# Check projection
raster("data/input_rasters/50/fw_species/Acipenser_baerii.tif")

# Surrogate

# REPTILES (SURROGATE)
 
names(threatened_reptiles_pa_matrix50) <- gsub(x = names(threatened_reptiles_pa_matrix50), 
                                         pattern = "\\.", replacement = "_")

purrr::map(names(threatened_reptiles_pa_matrix50)[-c(1:2)], function(sp){
  r <- rasterFromXYZ(threatened_reptiles_pa_matrix50 %>% dplyr::select(x, y, sp))
  proj4string(r) <- "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs"
  writeRaster(r, filename = paste0("data/input_rasters/50/reptiles/", sp, ".tif"), format = "GTiff")
})

# Check projection
raster("data/input_rasters/50/reptiles/Abronia_anzuetoi.tif")

# BIRDS (SURROGATE) 
names(threatened_birds_pa_matrix50) <- gsub(x = names(threatened_birds_pa_matrix50), 
pattern = "\\.", replacement = "_")

purrr::map(names(threatened_birds_pa_matrix50)[-c(1:2)], function(sp){
  r <- rasterFromXYZ(threatened_birds_pa_matrix50 %>% dplyr::select(x, y, sp))
  proj4string(r) <- "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs"
  writeRaster(r, filename = paste0("data/input_rasters/50/birds/", sp, ".tif"), format = "GTiff")
})

#Check projection
raster("data/input_rasters/50/birds/Acanthiza_katherina.tif")

# MAMMALS (SURROGATE)
names(threatened_mammals_pa_matrix50) <- gsub(x = names(threatened_mammals_pa_matrix50), 
                                            pattern = "\\.", replacement = "_")

purrr::map(names(threatened_mammals_pa_matrix50)[-c(1:2)], function(sp){
  r <- rasterFromXYZ(threatened_mammals_pa_matrix50 %>% dplyr::select(x, y, sp))
  proj4string(r) <- "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs"
  writeRaster(r, filename = paste0("data/input_rasters/50/mammals/", sp, ".tif"), format = "GTiff")
})

#Check projection
raster("data/input_rasters/50/mammals/Acinonyx_jubatus.tif")

# AMPHIBIANS (SURROGATE)
names(threatened_amphibians_pa_matrix50) <- gsub(x = names(threatened_amphibians_pa_matrix50), 
                                              pattern = "\\.", replacement = "_")

purrr::map(names(threatened_amphibians_pa_matrix50)[-c(1:2)], function(sp){
  r <- rasterFromXYZ(threatened_amphibians_pa_matrix50 %>% dplyr::select(x, y, sp))
  proj4string(r) <- "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs"
  writeRaster(r, filename = paste0("data/input_rasters/50/amphibians/", sp, ".tif"), format = "GTiff")
})

#Check projection
raster("data/input_rasters/50/amphibians/Hyperolius_nimbae.tif")

# 2. Run Zonation ####

# Install R zonator package, if not already installed
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

# Create directory to store outputs
setwd("F:/surrogacy_fw_species/Surrogacy_FW_PART2")
getwd()
dir.create("output")
dir.create("output/zonation_runs")
dir.create("output/zonation_runs/50")

# ZONATION ABF AND CAZ METHOD 
# (ABF stands for Additive Benefit Function algorithm and CAZ stands for Core Area Zonation algorithm)

#Surrogates:

# REPTILES
run_zonation(id = "reptiles", runs = 5, algorithm = "CAZ", out_dir = "output/zonation_runs/50", input_rasters_dir = "data/input_rasters/50")
run_zonation(id = "reptiles", runs = 5, algorithm = "ABF", out_dir = "output/zonation_runs/50", input_rasters_dir = "data/input_rasters/50")

# BIRDS
run_zonation(id = "birds", runs = 5, algorithm = "CAZ", out_dir = "output/zonation_runs/50", input_rasters_dir = "data/input_rasters/50")
run_zonation(id = "birds", runs = 5, algorithm = "ABF", out_dir = "output/zonation_runs/50", input_rasters_dir = "data/input_rasters/50")

# AMPHIBIANS
run_zonation(id = "amphibians", runs = 5, algorithm = "CAZ", out_dir = "output/zonation_runs/50", input_rasters_dir = "data/input_rasters/50")

run_zonation(id = "amphibians", runs = 5, algorithm = "ABF", out_dir = "output/zonation_runs/50", input_rasters_dir = "data/input_rasters/50")

# MAMMALS
run_zonation(id = "mammals", runs = 5, algorithm = "CAZ", out_dir = "output/zonation_runs/50", input_rasters_dir = "data/input_rasters/50")
run_zonation(id = "mammals", runs = 5, algorithm = "ABF", out_dir = "output/zonation_runs/50", input_rasters_dir = "data/input_rasters/50")

# COMBINED (REPTILES, AMPHIBIANS, BIRDS AND MAMMALS)

run_zonation(id = "combined", runs = 5, algorithm = "CAZ", out_dir = "output/zonation_runs/50", input_rasters_dir = "data/input_rasters/50")

run_zonation(id = "combined", runs = 5, algorithm = "ABF", out_dir = "output/zonation_runs/50", input_rasters_dir = "data/input_rasters/50")

#Targets:

# FW SPECIES 
run_zonation(id = "fw_species", runs = 5, algorithm = "CAZ", out_dir = "output/zonation_runs/50", input_rasters_dir = "data/input_rasters/50")

run_zonation(id = "fw_species", runs = 5, algorithm = "ABF", out_dir = "output/zonation_runs/50", input_rasters_dir = "data/input_rasters/50")

# ODONATA

run_zonation(id = "odonata", runs = 5, algorithm = "CAZ", out_dir = "output/zonation_runs/50", input_rasters_dir = "data/input_rasters/50")
run_zonation(id = "odonata", runs = 5, algorithm = "ABF", out_dir = "output/zonation_runs/50", input_rasters_dir = "data/input_rasters/50")

# DECAPODS

run_zonation(id = "decapods", runs = 5, algorithm = "CAZ", out_dir = "output/zonation_runs/50", input_rasters_dir = "data/input_rasters/50")
run_zonation(id = "decapods", runs = 5, algorithm = "ABF", out_dir = "output/zonation_runs/50", input_rasters_dir = "data/input_rasters/50")

#FW FISHES

run_zonation(id = "fishes", runs = 5, algorithm = "CAZ", out_dir = "output/zonation_runs/50", input_rasters_dir = "data/input_rasters/50")
run_zonation(id = "fishes", runs = 5, algorithm = "ABF", out_dir = "output/zonation_runs/50", input_rasters_dir = "data/input_rasters/50")

# At this point move to other codes if you want to proceed with analysis of individual groups of decapods, FW fishes and odonates as targets, the following is for FW_SPECIES Combined (decapods + FW fishes and odonates) as target

#3. Extract Zonation ranks ####

# CAZ algorithm
ranks_CAZ_fw_species <- extract_zonation_ranks(taxon = "fw_species", algorithm = "CAZ", out_dir = "output/zonation_runs/50")

ranks_CAZ_reptiles <- extract_zonation_ranks(taxon = "reptiles", algorithm = "CAZ", out_dir = "output/zonation_runs/50")

ranks_CAZ_birds <- extract_zonation_ranks(taxon = "birds", algorithm = "CAZ", out_dir = "output/zonation_runs/50")

ranks_CAZ_mammals <- extract_zonation_ranks(taxon = "mammals", algorithm = "CAZ", out_dir = "output/zonation_runs/50")

ranks_CAZ_amphibians <- extract_zonation_ranks(taxon = "amphibians", algorithm = "CAZ", out_dir = "output/zonation_runs/50")

ranks_CAZ_combined <- extract_zonation_ranks(taxon = "combined", algorithm = "CAZ", out_dir = "output/zonation_runs/50")

# ABF algorithm
ranks_ABF_fw_species <- extract_zonation_ranks(taxon = "fw_species", algorithm = "ABF", out_dir = "output/zonation_runs/50")

ranks_ABF_reptiles <- extract_zonation_ranks(taxon = "reptiles", algorithm = "ABF", out_dir = "output/zonation_runs/50")

ranks_ABF_birds <- extract_zonation_ranks(taxon = "birds", algorithm = "ABF", out_dir = "output/zonation_runs/50")

ranks_ABF_mammals <- extract_zonation_ranks(taxon = "mammals", algorithm = "ABF", out_dir = "output/zonation_runs/50")

ranks_ABF_amphibians <- extract_zonation_ranks(taxon = "amphibians", algorithm = "ABF", out_dir = "output/zonation_runs/50")

ranks_ABF_combined <- extract_zonation_ranks(taxon = "combined", algorithm = "ABF", out_dir = "output/zonation_runs/50")

#4. Generate ranks based on randomized cell sequences ####

ranks_random_fw_species <- random_ranks(ranks_CAZ_fw_species, runs = 100)

ranks_random_comb <- random_ranks(ranks_CAZ_combined, runs = 100)
ranks_random_rept <- random_ranks(ranks_CAZ_reptiles, runs = 100)
ranks_random_bird <- random_ranks(ranks_CAZ_birds, runs = 100)
ranks_random_mamm <- random_ranks(ranks_CAZ_mammals, runs = 100)
ranks_random_amph <- random_ranks(ranks_CAZ_amphibians, runs = 100)

#5. Calculate terrestrial extents over which to assess surrogacy

# Create vertebrate classes object
#vert_classes <- c("fw_species", "combined")
vert_classes <- c("birds", "reptiles", "mammals", "amphibians", "fw_species", "combined")

# Extract Zonation rank maps

# Create a new folder with the following name xxxx_projected_062023

map_extent <- purrr::map(vert_classes, extract_zonation_maps, algorithm = "CAZ", out_dir = "output/zonation_runs", suffix = "062023", resolution = "50") %>% 
  purrr::set_names(vert_classes)

map_extent$amphibians

# Extract combined terrestrial extent over which are any threatened species occur, regardless of class

terrestrial_extent <- map_extent[["birds"]] > 0 | map_extent[["fw_species"]] > 0 | map_extent[["mammals"]] > 0 | map_extent[["amphibians"]] > 0 | map_extent[["reptiles"]] > 0

terrestrial_extent <- as.data.frame(terrestrial_extent, xy = TRUE) %>% 
  dplyr::mutate(cellID = 1:length(terrestrial_extent[])) %>% 
  dplyr::filter(complete.cases(.))

ranks_random_terrestrial <- random_ranks(terrestrial_extent$cellID, runs = 100)

# FW species + Combined extent
ranks_random_fw_species_comb <- random_ranks(as.data.frame(map_extent[["combined"]] > 0 | map_extent[["fw_species"]] > 0, xy = TRUE) %>% dplyr::mutate(cellID = 1:length(x[])) %>% dplyr::filter(complete.cases(.)) %>% pull(cellID), runs = 100)

# FW species + Birds extent
ranks_random_fw_species_bird <- random_ranks(as.data.frame(map_extent[["birds"]] > 0 | map_extent[["fw_species"]] > 0, xy = TRUE) %>% 
                                         dplyr::mutate(cellID = 1:length(x[])) %>% 
                                         dplyr::filter(complete.cases(.)) %>% pull(cellID), runs = 100)

# FW species + Mammals extent
ranks_random_fw_species_mammals <- random_ranks(as.data.frame(map_extent[["mammals"]] > 0 | map_extent[["fw_species"]] > 0, xy = TRUE) %>% 
                                               dplyr::mutate(cellID = 1:length(x[])) %>% 
                                               dplyr::filter(complete.cases(.)) %>% pull(cellID), runs = 100)

# FW species + Amphibians extent
ranks_random_fw_species_amphibians <- random_ranks(as.data.frame(map_extent[["amphibians"]] > 0 | map_extent[["fw_species"]] > 0, xy = TRUE) %>% 
                                                  dplyr::mutate(cellID = 1:length(x[])) %>% 
                                                  dplyr::filter(complete.cases(.)) %>% pull(cellID), runs = 100)

# FW species + Reptiles extent
ranks_random_fw_species_reptiles <- random_ranks(as.data.frame(map_extent[["reptiles"]] > 0 | map_extent[["fw_species"]] > 0, xy = TRUE) %>% 
                                                     dplyr::mutate(cellID = 1:length(x[])) %>% 
                                                     dplyr::filter(complete.cases(.)) %>% pull(cellID), runs = 100)

# 6. Generate Surrogacy Curves ####

# CAZ

# FW species as Target

# Optimal curve (FW species)
curve_CAZ_fw_species_threat_fw_species_RAW <- get_sai_curves_target(
  target_matrix = threatened_fw_pa_matrix50, 
  target_ranks = ranks_CAZ_fw_species, 
  algorithm = "CAZ", 
  runs = 5)

saveRDS(curve_CAZ_fw_species_threat_fw_species, "output/zonation_runs/50/curve_CAZ_fw_species_threat_fw_species.rds")
curve_CAZ_fw_species_threat_fw_species <- readRDS("output/zonation_runs/50/curve_CAZ_fw_species_threat_fw_species.rds")

# Combine (amphibians + birds + mammals + reptiles) as surrogate
curve_CAZ_fw_species_threat_comb <- get_sai_curves_surrogate(target_matrix = threatened_fw_pa_matrix50, 
                                                             surrogate_ranks = ranks_CAZ_combined, 
                                                             algorithm = "CAZ",
                                                             runs = 5)

saveRDS(curve_CAZ_fw_species_threat_comb, "output/zonation_runs/50/curve_CAZ_fw_species_threat_comb.rds")
curve_CAZ_fw_species_threat_comb <- readRDS("output/zonation_runs/50/curve_CAZ_fw_species_threat_comb.rds")

# Birds as surrogate
curve_CAZ_fw_species_threat_bird <- get_sai_curves_surrogate(target_matrix = threatened_fw_pa_matrix50, 
                                                             surrogate_ranks = ranks_CAZ_birds, 
                                                             algorithm = "CAZ",
                                                             runs = 5)

saveRDS(curve_CAZ_fw_species_threat_bird, "output/zonation_runs/50/curve_CAZ_fw_species_threat_bird.rds")
curve_CAZ_fw_species_threat_bird <- readRDS("output/zonation_runs/50/curve_CAZ_fw_species_threat_bird.rds")

# Mammals as surrogate
curve_CAZ_fw_species_threat_mamm <- get_sai_curves_surrogate(target_matrix = threatened_fw_pa_matrix50, 
                                                             surrogate_ranks = ranks_CAZ_mammals, 
                                                             algorithm = "CAZ",
                                                             runs = 5)

saveRDS(curve_CAZ_fw_species_threat_mamm, "output/zonation_runs/50/curve_CAZ_fw_species_threat_mamm.rds")
curve_CAZ_fw_species_threat_mamm <- readRDS("output/zonation_runs/50/curve_CAZ_fw_species_threat_mamm.rds")

# Amphibians as surrogate
curve_CAZ_fw_species_threat_amph <- get_sai_curves_surrogate(target_matrix = threatened_fw_pa_matrix50, 
                                                             surrogate_ranks = ranks_CAZ_amphibians, 
                                                             algorithm = "CAZ",
                                                             runs = 5)

saveRDS(curve_CAZ_fw_species_threat_amph, "output/zonation_runs/50/curve_CAZ_fw_species_threat_amph.rds")
curve_CAZ_fw_species_threat_amph <- readRDS("output/zonation_runs/50/curve_CAZ_fw_species_threat_amph.rds")

# Reptiles as surrogate
curve_CAZ_fw_species_threat_rept <- get_sai_curves_surrogate(target_matrix = threatened_fw_pa_matrix50, 
                                                             surrogate_ranks = ranks_CAZ_reptiles, 
                                                             algorithm = "CAZ",
                                                             runs = 5)

saveRDS(curve_CAZ_fw_species_threat_rept, "output/zonation_runs/50/curve_CAZ_fw_species_threat_rept.rds")
curve_CAZ_fw_species_threat_rept <- readRDS("output/zonation_runs/50/curve_CAZ_fw_species_threat_rept.rds")

# ABF

# FW species as target

# Optimal curve (FW species)

curve_ABF_fw_species_threat_fw_species <- get_sai_curves_target(
  target_matrix = threatened_fw_pa_matrix50, 
  target_ranks = ranks_ABF_fw_species, 
  algorithm = "ABF", 
  runs = 5)

saveRDS(curve_ABF_fw_species_threat_fw_species, "output/zonation_runs/50/curve_ABF_fw_species_threat_fw_species.rds")
curve_ABF_fw_species_threat_fw_species <- readRDS("output/zonation_runs/50/curve_ABF_fw_species_threat_fw_species.rds")

# Combine (amphibians + birds + mammals + reptiles) as surrogate
curve_ABF_fw_species_threat_comb <- get_sai_curves_surrogate(target_matrix = threatened_fw_pa_matrix50, 
                                                             surrogate_ranks = ranks_ABF_combined, 
                                                             algorithm = "ABF", 
                                                             runs = 5)

saveRDS(curve_ABF_fw_species_threat_comb, 
        "output/zonation_runs/50/curve_ABF_fw_species_threat_comb.rds")
curve_ABF_fw_species_threat_comb <- readRDS("output/zonation_runs/50/curve_ABF_fw_species_threat_comb.rds")

# Birds as surrogate
curve_ABF_fw_species_threat_bird <- get_sai_curves_surrogate(target_matrix = threatened_fw_pa_matrix50, 
                                                             surrogate_ranks = ranks_ABF_birds, 
                                                             algorithm = "ABF",
                                                             runs = 5)

saveRDS(curve_ABF_fw_species_threat_bird, "output/zonation_runs/50/curve_ABF_fw_species_threat_bird.rds")
curve_ABF_fw_species_threat_bird <- readRDS("output/zonation_runs/50/curve_ABF_fw_species_threat_bird.rds")

# Mammals as surrogate
curve_ABF_fw_species_threat_mamm <- get_sai_curves_surrogate(target_matrix = threatened_fw_pa_matrix50, 
                                                             surrogate_ranks = ranks_ABF_mammals, 
                                                             algorithm = "ABF",
                                                             runs = 5)

saveRDS(curve_ABF_fw_species_threat_mamm, "output/zonation_runs/50/curve_ABF_fw_species_threat_mamm.rds")
curve_ABF_fw_species_threat_mamm <- readRDS("output/zonation_runs/50/curve_ABF_fw_species_threat_mamm.rds")

# Amphibians as surrogate
curve_ABF_fw_species_threat_amph <- get_sai_curves_surrogate(target_matrix = threatened_fw_pa_matrix50, 
                                                             surrogate_ranks = ranks_ABF_amphibians, 
                                                             algorithm = "ABF",
                                                             runs = 5)

saveRDS(curve_ABF_fw_species_threat_amph, "output/zonation_runs/50/curve_ABF_fw_species_threat_amph.rds")
curve_ABF_fw_species_threat_amph <- readRDS("output/zonation_runs/50/curve_ABF_fw_species_threat_amph.rds")

# Reptiles as surrogate
curve_ABF_fw_species_threat_rept <- get_sai_curves_surrogate(target_matrix = threatened_fw_pa_matrix50, 
                                                             surrogate_ranks = ranks_ABF_reptiles, 
                                                             algorithm = "ABF",
                                                             runs = 5)

saveRDS(curve_ABF_fw_species_threat_rept, "output/zonation_runs/50/curve_ABF_fw_species_threat_rept.rds")
curve_ABF_fw_species_threat_rept <- readRDS("output/zonation_runs/50/curve_ABF_fw_species_threat_rept.rds")

# 7. Generate random expectations ####

# FW species as target

# CAZ

# Terrestrial extent
curve_CAZ_fw_species_threat_random_terr_ext <- 
  get_sai_curves_random(target_matrix = threatened_fw_pa_matrix50, 
                        random_ranks = ranks_random_terrestrial, 
                        reference_coordinates = ranks_CAZ_fw_species$reference_coordinates, 
                        algorithm = "CAZ", 
                        runs = 10)

saveRDS(curve_CAZ_fw_species_threat_random_terr_ext, 
        "output/zonation_runs/50/curve_CAZ_fw_species_threat_random_terr_ext.rds")
curve_CAZ_fw_species_threat_random_terr_ext <- readRDS(
  "output/zonation_runs/50/curve_CAZ_fw_species_threat_random_terr_ext.rds")

# Birds as surrogate
curve_CAZ_fw_species_threat_random_bird <- get_sai_curves_random(target_matrix = threatened_fw_pa_matrix50, random_ranks = ranks_random_fw_species_bird, reference_coordinates = ranks_CAZ_fw_species$reference_coordinates, algorithm = "CAZ", runs = 10)

saveRDS(curve_CAZ_fw_species_threat_random_bird, "output/zonation_runs/50/curve_CAZ_fw_species_threat_random_bird.rds")
curve_CAZ_fw_species_threat_random_bird <- readRDS("output/zonation_runs/50/curve_CAZ_fw_species_threat_random_bird.rds")

# Mammals as surrogate
curve_CAZ_fw_species_threat_random_mamm <- get_sai_curves_random(target_matrix =  threatened_fw_pa_matrix50, random_ranks = ranks_random_fw_species_mammals, reference_coordinates = ranks_CAZ_fw_species$reference_coordinates, algorithm = "CAZ", runs = 10)

saveRDS(curve_CAZ_fw_species_threat_random_mamm, "output/zonation_runs/50/curve_CAZ_fw_species_threat_random_mamm.rds")
curve_CAZ_fw_species_threat_random_mamm <- readRDS("output/zonation_runs/50/curve_CAZ_fw_species_threat_random_mamm.rds")

# Amphibians as surrogate
curve_CAZ_fw_species_threat_random_amph <- get_sai_curves_random(target_matrix =threatened_fw_pa_matrix50, random_ranks = ranks_random_fw_species_amphibians, reference_coordinates = ranks_CAZ_fw_species$reference_coordinates, algorithm = "CAZ", runs = 10)

saveRDS(curve_CAZ_fw_species_threat_random_amph, "output/zonation_runs/50/curve_CAZ_fw_species_threat_random_amph.rds")
curve_CAZ_fw_species_threat_random_amph <- readRDS("output/zonation_runs/50/curve_CAZ_fw_species_threat_random_amph.rds")

# Reptiles as surrogate
curve_CAZ_fw_species_threat_random_rept <- get_sai_curves_random(target_matrix =threatened_fw_pa_matrix50, random_ranks = ranks_random_fw_species_reptiles, reference_coordinates = ranks_CAZ_fw_species$reference_coordinates, algorithm = "CAZ", runs = 10)

saveRDS(curve_CAZ_fw_species_threat_random_rept, "output/zonation_runs/50/curve_CAZ_fw_species_threat_random_rept.rds")
curve_CAZ_fw_species_threat_random_rept <- readRDS("output/zonation_runs/50/curve_CAZ_fw_species_threat_random_rept.rds")

# Combine (amphibians + birds + mammals + reptiles) as surrogate

curve_CAZ_fw_species_threat_random_comb <- get_sai_curves_random(target_matrix = threatened_fw_pa_matrix50,
                                                                 random_ranks = ranks_random_fw_species_comb, 
                                                                 reference_coordinates = ranks_CAZ_fw_species$reference_coordinates, 
                                                                 algorithm = "CAZ", 
                                                                 runs = 10)

saveRDS(curve_CAZ_fw_species_threat_random_comb, 
        "output/zonation_runs/50/curve_CAZ_fw_species_threat_random_comb.rds")
curve_CAZ_fw_species_threat_random_comb <- readRDS(
  "output/zonation_runs/50/curve_CAZ_fw_species_threat_random_comb.rds")

# ABF

# Terrestrial extent
curve_ABF_fw_species_threat_random_terr_ext <- get_sai_curves_random(target_matrix = threatened_fw_pa_matrix50, 
                                                                     random_ranks = ranks_random_terrestrial, 
                                                                     reference_coordinates = ranks_ABF_fw_species$reference_coordinates, 
                                                                     algorithm = "ABF", 
                                                                     runs = 10)

saveRDS(curve_ABF_fw_species_threat_random_terr_ext , 
        "output/zonation_runs/50/curve_ABF_fw_species_threat_random_terr_ext.rds")
curve_ABF_fw_species_threat_random_terr_ext  <- readRDS(
  "output/zonation_runs/50/curve_ABF_fw_species_threat_random_terr_ext.rds")

# Birds as surrogate
curve_ABF_fw_species_threat_random_bird <- get_sai_curves_random(target_matrix = threatened_fw_pa_matrix50, random_ranks = ranks_random_fw_species_bird, reference_coordinates = ranks_ABF_fw_species$reference_coordinates, algorithm = "ABF", runs = 10)

saveRDS(curve_ABF_fw_species_threat_random_bird, "output/zonation_runs/50/curve_ABF_fw_species_threat_random_bird.rds")
curve_ABF_fw_species_threat_random_bird <- readRDS("output/zonation_runs/50/curve_ABF_fw_species_threat_random_bird.rds")

# Mammals as surrogate
curve_ABF_fw_species_threat_random_mamm <- get_sai_curves_random(target_matrix =  threatened_fw_pa_matrix50, random_ranks = ranks_random_fw_species_mammals, reference_coordinates = ranks_ABF_fw_species$reference_coordinates, algorithm = "ABF", runs = 10)

saveRDS(curve_ABF_fw_species_threat_random_mamm, "output/zonation_runs/50/curve_ABF_fw_species_threat_random_mamm.rds")
curve_ABF_fw_species_threat_random_mamm <- readRDS("output/zonation_runs/50/curve_ABF_fw_species_threat_random_mamm.rds")

# Amphibians as surrogate
curve_ABF_fw_species_threat_random_amph <- get_sai_curves_random(target_matrix =threatened_fw_pa_matrix50, random_ranks = ranks_random_fw_species_amphibians, reference_coordinates = ranks_ABF_fw_species$reference_coordinates, algorithm = "ABF", runs = 10)

saveRDS(curve_ABF_fw_species_threat_random_amph, "output/zonation_runs/50/curve_ABF_fw_species_threat_random_amph.rds")
curve_ABF_fw_species_threat_random_amph <- readRDS("output/zonation_runs/50/curve_ABF_fw_species_threat_random_amph.rds")

# Reptiles as surrogate
curve_ABF_fw_species_threat_random_rept <- get_sai_curves_random(target_matrix =threatened_fw_pa_matrix50, random_ranks = ranks_random_fw_species_reptiles, reference_coordinates = ranks_ABF_fw_species$reference_coordinates, algorithm = "ABF", runs = 10)

saveRDS(curve_ABF_fw_species_threat_random_rept, "output/zonation_runs/50/curve_ABF_fw_species_threat_random_rept.rds")
curve_ABF_fw_species_threat_random_rept <- readRDS("output/zonation_runs/50/curve_ABF_fw_species_threat_random_rept.rds")

# Combine (amphibians + birds + mammals + reptiles) as surrogate

curve_ABF_fw_species_threat_random_comb <- get_sai_curves_random(target_matrix = threatened_fw_pa_matrix50, 
                                                                 random_ranks = ranks_random_fw_species_comb, 
                                                                 reference_coordinates = ranks_ABF_fw_species$reference_coordinates, 
                                                                 algorithm = "ABF", runs = 10)

saveRDS(curve_ABF_fw_species_threat_random_comb, "output/zonation_runs/50/curve_ABF_fw_species_threat_random_comb.rds")
curve_ABF_rept_threat_random_comb <- readRDS("output/zonation_runs/50/curve_ABF_fw_species_threat_random_comb.rds")

# 8. Generate summary outputs ####
library(pracma);library(Rmisc);library(ggplot2)

# CAZ

# FW SPECIES as target

# Combined as surrogate
CAZ_fw_species_comb_threat_sai_curves <- combine_sai_curves(
  target = curve_CAZ_fw_species_threat_fw_species[1:nrow(curve_CAZ_fw_species_threat_random_terr_ext), ],
  surrogate = curve_CAZ_fw_species_threat_comb[1:nrow(curve_CAZ_fw_species_threat_random_terr_ext), ],
  random = curve_CAZ_fw_species_threat_random_terr_ext)

CAZ_fw_comb_threat_sai <- calculate_sai(CAZ_fw_species_comb_threat_sai_curves)[[1]]
CAZ_fw_comb_threat_sai_plot <- plot_sai_curves(CAZ_fw_species_comb_threat_sai_curves)

ggsave(CAZ_fw_comb_threat_sai_plot,
       filename = "CAZ_COMBINED_PLOT.pdf",
       device = "pdf", path = "./",
       width = 95, height = 95, units="mm",dpi = 300)  

# Birds as surrogate
CAZ_fw_species_birds_threat_sai_curves <- combine_sai_curves(
  target = curve_CAZ_fw_species_threat_fw_species[1:nrow(curve_CAZ_fw_species_threat_random_terr_ext), ],
  surrogate = curve_CAZ_fw_species_threat_bird[1:nrow(curve_CAZ_fw_species_threat_random_terr_ext), ],
  random = curve_CAZ_fw_species_threat_random_terr_ext)

CAZ_fw_birds_threat_sai <- calculate_sai(CAZ_fw_species_birds_threat_sai_curves)[[1]]
CAZ_fw_birds_threat_sai_plot <- plot_sai_curves(CAZ_fw_species_birds_threat_sai_curves)

ggsave(CAZ_fw_birds_threat_sai_plot,
       filename = "CAZ_BIRDS_PLOT.pdf",
       device = "pdf", path = "./",
       width = 95, height = 95, units="mm",dpi = 300)  

# Amphibians as surrogate
CAZ_fw_species_amphibians_threat_sai_curves <- combine_sai_curves(
  target = curve_CAZ_fw_species_threat_fw_species[1:nrow(curve_CAZ_fw_species_threat_random_terr_ext), ],
  surrogate = curve_CAZ_fw_species_threat_amph[1:nrow(curve_CAZ_fw_species_threat_random_terr_ext), ],
  random = curve_CAZ_fw_species_threat_random_terr_ext)


CAZ_fw_amphibians_threat_sai <- calculate_sai(CAZ_fw_species_amphibians_threat_sai_curves)[[1]]
CAZ_fw_amphibians_threat_sai_plot <- plot_sai_curves(CAZ_fw_species_amphibians_threat_sai_curves)

ggsave(CAZ_fw_amphibians_threat_sai_plot,
       filename = "CAZ_AMPHIBIANS_PLOT.pdf",
       device = "pdf", path = "./",
       width = 95, height = 95, units="mm",dpi = 300)  

# Reptiles as surrogate
CAZ_fw_species_reptiles_threat_sai_curves <- combine_sai_curves(
  target = curve_CAZ_fw_species_threat_fw_species[1:nrow(curve_CAZ_fw_species_threat_random_terr_ext), ],
  surrogate = curve_CAZ_fw_species_threat_rept[1:nrow(curve_CAZ_fw_species_threat_random_terr_ext), ],
  random = curve_CAZ_fw_species_threat_random_terr_ext)

CAZ_fw_reptiles_threat_sai <- calculate_sai(CAZ_fw_species_reptiles_threat_sai_curves)[[1]]
CAZ_fw_reptiles_threat_sai_plot <- plot_sai_curves(CAZ_fw_species_reptiles_threat_sai_curves)

ggsave(CAZ_fw_reptiles_threat_sai_plot,
       filename = "CAZ_REPTILES_PLOT.pdf",
       device = "pdf", path = "./",
       width = 95, height = 95, units="mm",dpi = 300)  

# Mammals as surrogate
CAZ_fw_species_mammals_threat_sai_curves <- combine_sai_curves(
  target = curve_CAZ_fw_species_threat_fw_species[1:nrow(curve_CAZ_fw_species_threat_random_terr_ext), ],
  surrogate = curve_CAZ_fw_species_threat_mamm[1:nrow(curve_CAZ_fw_species_threat_random_terr_ext), ],
  random = curve_CAZ_fw_species_threat_random_terr_ext)

CAZ_fw_mammals_threat_sai <- calculate_sai(CAZ_fw_species_mammals_threat_sai_curves)[[1]]
CAZ_fw_mammals_threat_sai_plot <- plot_sai_curves(CAZ_fw_species_mammals_threat_sai_curves)

ggsave(CAZ_fw_mammals_threat_sai_plot,
       filename = "CAZ_MAMMALS_PLOT.pdf",
       device = "pdf", path = "./",
       width = 95, height = 95, units="mm",dpi = 300)  

# ABF

# FW SPECIES as target

# Combined as surrogate
ABF_fw_species_comb_threat_sai_curves <- combine_sai_curves(
  target = curve_ABF_fw_species_threat_fw_species[1:nrow(curve_ABF_fw_species_threat_random_terr_ext), ],
  surrogate = curve_ABF_fw_species_threat_comb[1:nrow(curve_ABF_fw_species_threat_random_terr_ext), ],
  random = curve_ABF_fw_species_threat_random_terr_ext)

ABF_fw_comb_threat_sai <- calculate_sai(ABF_fw_species_comb_threat_sai_curves)[[1]]
ABF_fw_comb_threat_sai_plot <- plot_sai_curves(ABF_fw_species_comb_threat_sai_curves)

ggsave(ABF_fw_comb_threat_sai_plot,
       filename = "ABF_COMBINED_PLOT.pdf",
       device = "pdf", path = "./",
       width = 95, height = 95, units="mm",dpi = 300)  

# Birds as surrogate
ABF_fw_species_birds_threat_sai_curves <- combine_sai_curves(
  target = curve_ABF_fw_species_threat_fw_species[1:nrow(curve_ABF_fw_species_threat_random_terr_ext), ],
  surrogate = curve_ABF_fw_species_threat_bird[1:nrow(curve_ABF_fw_species_threat_random_terr_ext), ],
  random = curve_ABF_fw_species_threat_random_terr_ext)

ABF_fw_birds_threat_sai <- calculate_sai(ABF_fw_species_birds_threat_sai_curves)[[1]]
ABF_fw_birds_threat_sai_plot <- plot_sai_curves(ABF_fw_species_birds_threat_sai_curves)

ggsave(ABF_fw_birds_threat_sai_plot,
       filename = "ABF_BIRDS_PLOT.pdf",
       device = "pdf", path = "./",
       width = 95, height = 95, units="mm",dpi = 300)  

# Amphibians as surrogate
ABF_fw_species_amphibians_threat_sai_curves <- combine_sai_curves(
  target = curve_ABF_fw_species_threat_fw_species[1:nrow(curve_ABF_fw_species_threat_random_terr_ext), ],
  surrogate = curve_ABF_fw_species_threat_amph[1:nrow(curve_ABF_fw_species_threat_random_terr_ext), ],
  random = curve_ABF_fw_species_threat_random_terr_ext)

ABF_fw_amphibians_threat_sai <- calculate_sai(ABF_fw_species_amphibians_threat_sai_curves)[[1]]
ABF_fw_amphibians_threat_sai_plot <- plot_sai_curves(ABF_fw_species_amphibians_threat_sai_curves)

ggsave(ABF_fw_amphibians_threat_sai_plot,
       filename = "ABF_AMPHIBIANS_PLOT.pdf",
       device = "pdf", path = "./",
       width = 95, height = 95, units="mm",dpi = 300)  

# Reptiles as surrogate
ABF_fw_species_reptiles_threat_sai_curves <- combine_sai_curves(
  target = curve_ABF_fw_species_threat_fw_species[1:nrow(curve_ABF_fw_species_threat_random_terr_ext), ],
  surrogate = curve_ABF_fw_species_threat_rept[1:nrow(curve_ABF_fw_species_threat_random_terr_ext), ],
  random = curve_ABF_fw_species_threat_random_terr_ext)

ABF_fw_reptiles_threat_sai <- calculate_sai(ABF_fw_species_reptiles_threat_sai_curves)[[1]]
ABF_fw_reptiles_threat_sai_plot <- plot_sai_curves(ABF_fw_species_reptiles_threat_sai_curves)

ggsave(ABF_fw_reptiles_threat_sai_plot,
       filename = "ABF_REPTILES_PLOT.pdf",
       device = "pdf", path = "./",
       width = 95, height = 95, units="mm",dpi = 300)  

# Mammals as surrogate
ABF_fw_species_mammals_threat_sai_curves <- combine_sai_curves(
  target = curve_ABF_fw_species_threat_fw_species[1:nrow(curve_ABF_fw_species_threat_random_terr_ext), ],
  surrogate = curve_ABF_fw_species_threat_mamm[1:nrow(curve_ABF_fw_species_threat_random_terr_ext), ],
  random = curve_ABF_fw_species_threat_random_terr_ext)

ABF_fw_mammals_threat_sai <- calculate_sai(ABF_fw_species_mammals_threat_sai_curves)[[1]]
ABF_fw_mammals_threat_sai_plot <- plot_sai_curves(ABF_fw_species_mammals_threat_sai_curves)

ggsave(ABF_fw_mammals_threat_sai_plot,
       filename = "ABF_MAMMALS_PLOT.pdf",
       device = "pdf", path = "./",
       width = 95, height = 95, units="mm",dpi = 300)  
