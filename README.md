# FW-surrogacy

The code assess the extent to which the conservation efforts targeting threatened reptiles, birds, mammals, and amphibians—either individually or collectively—can act as surrogates for the conservation of threatened freshwater species, including combined odonates, decapods, and freshwater fish. This evaluation is accomplished by calculating a species accumulation index (SAI) of surrogate effectiveness, as defined by Cox et al. in 2022. The analysis is carried out within the context of two primary global conservation strategies: 1) Maximizing rarity-weighted richness, and 2) Maximizing inclusion of range-restricted species.

In addition, the code includes an evaluation of water stress as a potential surrogacy strategy for conserving threatened freshwater species. This assessment involves using the species accumulation index (SAI) to gauge the efficacy of water stress in identifying sites that most effectively represent threatened freshwater species. The evaluation is conducted within the framework of two conservation strategies: maximizing rarity-weighted richness and maximizing inclusion of range-restricted species.

To implement this, the code generates a complementarity-based ranking of conservation values for the specified targets using their respective algorithms, across the landscape of interest. The rank order is determined based on the baseline water stress layer sourced from the Aqueduct Water Risk Atlas. This layer measures the ratio of total water demand (including domestic, industrial, irrigation, and livestock consumptive and non-consumptive uses) to available renewable surface and groundwater supplies, as outlined by Garnet et al. in 2015.

Surrogacy analysis was conducted using R software version 4.3.0. The analyses were executed with the spatial conservation-planning software Zonation (v3.1) and the R package "zonator" (v.0.6.0). To perform the analyses, additional packages, including raster v3.6, dplyr v1.1.2, pracma v2.4.2, and Rmisc v1.5.1, are also required.

Instructions and additional comments to run the analyses:

- The R scripts for conducting the analyses are well-commented.
- Execute the codes sequentially in order by following the numerical order indicated in the R script names. 
- Spatial data files required for the analyses are stored in matrices and .rds format, to be loaded into R using the provided R scripts.
- The R scripts include steps and code for installing additional dependencies, such as zonator, which is designed to work in conjunction with the Zonation software.

Please download these files, and set your working directory appropriately in R. 

<img src="https://zenodo.org/badge/DOI/10.5281/zenodo.10286100.svg">
