# This code will generate the curves for decapods taxa as target and combine (amphibians + birds + mammals + reptiles) as surrogate.

library(raster);library(dplyr);library(zonator)

setwd("") #Set proper working directory

#Load custom functions
source("functions_code/functions.R")

# 1. Load cell-by-species matrices
threatened_decapods_pa_matrix50 <- readRDS("data/decapods_t_matrix_zeros.rds")

# 2. Extract Zonation ranks ####

# CAZ algorithm
ranks_CAZ_decapods <- extract_zonation_ranks(taxon = "decapods", algorithm = "CAZ", out_dir = "output/zonation_runs/50")
ranks_CAZ_combined <- extract_zonation_ranks(taxon = "combined", algorithm = "CAZ", out_dir = "output/zonation_runs/50")

# ABF algorithm
ranks_ABF_decapods <- extract_zonation_ranks(taxon = "decapods", algorithm = "ABF", out_dir = "output/zonation_runs/50")
ranks_ABF_combined <- extract_zonation_ranks(taxon = "combined", algorithm = "ABF", out_dir = "output/zonation_runs/50") 

# 3. Generate ranks based on randomized cell sequences ####

ranks_random_decapods <- random_ranks(ranks_CAZ_decapods, runs = 100)
ranks_random_comb <- random_ranks(ranks_CAZ_combined, runs = 100) 

# 4. Calculate terrestrial extents over which to assess surrogacy

# Create vertebrate classes object
vert_classes <- c("decapods", "combined") 

# Extract Zonation rank maps

# Create a new folder with the following name xxxx_projected_062023

map_extent <- purrr::map(vert_classes, extract_zonation_maps, algorithm = "CAZ", out_dir = "output/zonation_runs", suffix = "062023", resolution = "50") %>% 
  purrr::set_names(vert_classes)

map_extent$decapods

# Extract combined terrestrial extent over which are any threatened species occur, regardless of class
terrestrial_extent <- map_extent[["combined"]] > 0 | map_extent[["decapods"]] > 0 

terrestrial_extent <- as.data.frame(terrestrial_extent, xy = TRUE) %>% 
  dplyr::mutate(cellID = 1:length(terrestrial_extent[])) %>% 
  dplyr::filter(complete.cases(.))

ranks_random_terrestrial <- random_ranks(terrestrial_extent$cellID, runs = 100)

# Decapods + Combined extent
ranks_random_decapods_comb <- random_ranks(as.data.frame(map_extent[["combined"]] > 0 | map_extent[["decapods"]] > 0, xy = TRUE) %>% 
                                            dplyr::mutate(cellID = 1:length(x[])) %>% dplyr::filter(complete.cases(.)) %>% 
                                            pull(cellID), runs = 100)

# 5. Generate Surrogacy Curves ####

# CAZ

# Decapods as target

# Optimal curve
curve_CAZ_decapods_threat_decapods <- get_sai_curves_target(
  target_matrix = threatened_decapods_pa_matrix50, 
  target_ranks = ranks_CAZ_decapods, 
  algorithm = "CAZ", 
  runs = 5)

saveRDS(curve_CAZ_decapods_threat_decapods, "output/zonation_runs/50/curve_CAZ_decapods_threat_decapods.rds")
curve_CAZ_decapods_threat_decapods <- readRDS("output/zonation_runs/50/curve_CAZ_decapods_threat_decapods.rds")

# Combine as surrogate
curve_CAZ_decapods_threat_comb <- get_sai_curves_surrogate(target_matrix = threatened_decapods_pa_matrix50, 
                                                          surrogate_ranks = ranks_CAZ_combined, 
                                                          algorithm = "CAZ",
                                                          runs = 5)

saveRDS(curve_CAZ_decapods_threat_comb, "output/zonation_runs/50/curve_CAZ_decapods_threat_comb.rds")
curve_CAZ_decapods_threat_comb <- readRDS("output/zonation_runs/50/curve_CAZ_decapods_threat_comb.rds")

# ABF

# Decapods as target

# Optimal curve
curve_ABF_decapods_threat_decapods <- get_sai_curves_target(
  target_matrix = threatened_decapods_pa_matrix50, 
  target_ranks = ranks_ABF_decapods, 
  algorithm = "ABF", 
  runs = 5)

saveRDS(curve_ABF_decapods_threat_decapods, "output/zonation_runs/50/curve_ABF_decapods_threat_decapods.rds")
curve_ABF_decapods_threat_decapods <- readRDS("output/zonation_runs/50/curve_ABF_decapods_threat_decapods.rds")

# Combine as surrogate
curve_ABF_decapods_threat_comb <- get_sai_curves_surrogate(target_matrix = threatened_decapods_pa_matrix50, 
                                                          surrogate_ranks = ranks_ABF_combined, 
                                                          algorithm = "ABF", 
                                                          runs = 5)

saveRDS(curve_ABF_decapods_threat_comb, 
        "output/zonation_runs/50/curve_ABF_decapods_threat_comb.rds")
curve_ABF_decapods_threat_comb <- readRDS("output/zonation_runs/50/curve_ABF_decapods_threat_comb.rds")

# 6. Generate random expectations ####

# Decapods as target

# CAZ

# Terrestrial extent
curve_CAZ_decapods_threat_random_terr_ext <- 
  get_sai_curves_random(target_matrix = threatened_decapods_pa_matrix50, 
                        random_ranks = ranks_random_terrestrial, 
                        reference_coordinates = ranks_CAZ_decapods$reference_coordinates, 
                        algorithm = "CAZ", 
                        runs = 10)

saveRDS(curve_CAZ_decapods_threat_random_terr_ext, 
        "output/zonation_runs/50/curve_CAZ_decapods_threat_random_terr_ext.rds")
curve_CAZ_decapods_threat_random_terr_ext <- readRDS(
  "output/zonation_runs/50/curve_CAZ_decapods_threat_random_terr_ext.rds")

# Combine as surrogate 
curve_CAZ_decapods_threat_random_comb <- get_sai_curves_random(target_matrix = threatened_decapods_pa_matrix50,
                                                              random_ranks = ranks_random_decapods_comb, 
                                                              reference_coordinates = ranks_CAZ_decapods$reference_coordinates, 
                                                              algorithm = "CAZ", 
                                                              runs = 10)

saveRDS(curve_CAZ_decapods_threat_random_comb, 
        "output/zonation_runs/50/curve_CAZ_decapods_threat_random_comb.rds")
curve_CAZ_decapods_threat_random_comb <- readRDS(
  "output/zonation_runs/50/curve_CAZ_decapods_threat_random_comb.rds")

# ABF
# Terrestrial extent
curve_ABF_decapods_threat_random_terr_ext <- get_sai_curves_random(target_matrix = threatened_decapods_pa_matrix50, 
                                                                  random_ranks = ranks_random_terrestrial, 
                                                                  reference_coordinates = ranks_ABF_decapods$reference_coordinates, 
                                                                  algorithm = "ABF", 
                                                                  runs = 10)

saveRDS(curve_ABF_decapods_threat_random_terr_ext , 
        "output/zonation_runs/50/curve_ABF_decapods_threat_random_terr_ext.rds")
curve_ABF_decapods_threat_random_terr_ext  <- readRDS(
  "output/zonation_runs/50/curve_ABF_decapods_threat_random_terr_ext.rds")

# Combine as surrogate
curve_ABF_decapods_threat_random_comb <- get_sai_curves_random(target_matrix = threatened_decapods_pa_matrix50, 
                                                              random_ranks = ranks_random_decapods_comb, 
                                                              reference_coordinates = ranks_ABF_decapods$reference_coordinates, 
                                                              algorithm = "ABF", runs = 10)

saveRDS(curve_ABF_decapods_threat_random_comb, "output/zonation_runs/50/curve_ABF_decapods_threat_random_comb.rds")
curve_ABF_decapods_threat_random_comb <- readRDS("output/zonation_runs/50/curve_ABF_decapods_threat_random_comb.rds")

# 7. Generate summary outputs
library(pracma);library(Rmisc);library(ggplot2)

# CAZ

# Decapods as target

# Combine as surrogate
CAZ_decapods_comb_threat_sai_curves <- combine_sai_curves(target = curve_CAZ_decapods_threat_decapods[1:nrow(curve_CAZ_decapods_threat_random_terr_ext), ],
                                                         surrogate = curve_CAZ_decapods_threat_comb[1:nrow(curve_CAZ_decapods_threat_random_terr_ext), ],
                                                         random = curve_CAZ_decapods_threat_random_terr_ext)

CAZ_decapods_comb_threat_sai <- calculate_sai(CAZ_decapods_comb_threat_sai_curves)[[1]]
CAZ_decapods_comb_threat_sai_plot <- plot_sai_curves(CAZ_decapods_comb_threat_sai_curves)

# ABF

# Decapods as target

# Combine as surrogate
ABF_decapods_comb_threat_sai_curves <- combine_sai_curves(target = curve_ABF_decapods_threat_decapods[1:nrow(curve_ABF_decapods_threat_random_terr_ext), ],
                                                         surrogate = curve_ABF_decapods_threat_comb[1:nrow(curve_ABF_decapods_threat_random_terr_ext), ],
                                                         random = curve_ABF_decapods_threat_random_terr_ext)

ABF_decapods_comb_threat_sai <- calculate_sai(ABF_decapods_comb_threat_sai_curves)[[1]]
ABF_decapods_comb_threat_sai_plot <- plot_sai_curves(ABF_decapods_comb_threat_sai_curves)

#Export figures

ggsave(ABF_decapods_comb_threat_sai_plot,
       filename = "ABF_PLOT_DECAPODS_COMBINED.pdf",
       device = "pdf", path = "./",
       width = 95, height = 95, units="mm",dpi = 300)  

ggsave(CAZ_decapods_comb_threat_sai_plot,
       filename = "CAZ_PLOT_DECAPODS_COMBINED.pdf",
       device = "pdf", path = "./",
       width = 95, height = 95, units="mm",dpi = 300)  
