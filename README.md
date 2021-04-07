# Materials for Modeling Travel and Importations on Bioko Island

This repository contains the code and documentation which support the data analysis and computer simulations for the publication "Quantifying malaria acquired during travel and its role in malaria elimination on Bioko Island" by Daniel T. Citron et al.

The supporting data which accompanies this code may be found [here](https://figshare.com/articles/dataset/Data_Supporting_Bioko_Island_Travel_Modeling/14380565).

* `aggregated_2015_2018_travel_data.csv` — contains population census data; trips recorded between 2015-2018 from the BIMEP MIS
* `travel_duration_data.csv` — contains data describing the duration of trips to different destinations from the 2018 BIMEP MIS
* `travel_dist_by_region.csv` — contains distances between each areaId and the centroid of each destination region
* `pfpr_draws.csv` — Draws from the joint posterior surface of malaria prevalence estimates (With permission, from [Guerra et al (2019)](https://www.nature.com/articles/s41467-019-10339-1)
* `negative_binomial_predictions_by_destination_region.csv` — Trip destination model results
* `trip_frequency_model_estimates.csv` — Frequency of travel model results

The R code inside the **scripts** folder is used to preprocess and analyze the data.

* `Trip_frequency_model.R` — Estimates the frequency with which residents leave home.
* `Trip_destination_choice_model.R` — Estimates the probability distribution of choosing a destination region, in a manner similar to a gravity model.
* `Trip_duration_model.R` — Fits an exponential model to the time spent away
* `region_to_areaId_mapping.R` — A script for disaggregating trips to regions into trips to specific map-areas within each region. We assume that the probability of visiting each map-area is proportional to population.

The R code insde the **simulation_configurations** folder serves as example code which maybe used to generate and plot the results from a single simulation run for each of the different scenarios discussed in the article.

* `baseline_simulation_configuration.R` — Baseline case, assuming that conditions do not change on Bioko Island
* `local_residual_configuration.R` — Estimating the local residual fraction, the fraction of cases that are attributable to on-island local transmission only, by reducing the force of infection on the mainland to zero.
* `travel_fraction_configuration.R` — Estimating the fraction of PfPR which is attributable to travel, by reducing the on-island local transmission to zero.
* `vaccination_configuration.R` — Estimating the long-term impact of a vaccine on malaria prevalence
* `travel_vaccination_configuration.R` — Estimating the impact of vaccinating and treating travelers to the mainland.

The software is called `macro.pfsi`, and is bundled as an R package in the directory **macro.pfsi** and as a compiled binary `macro.pfsi_1.0.tar.gz`. To install the `macro.pfsi` software, the easiest method is to use a terminal to navigate to the directory containing `macro.pfsi_1.0.tar.gz` and call `R CMD INSTALL macro.pfsi_1.0.tar.gz`. Alternatively, you can install it within an RStudio session using devtools `devtools::install_github(repo = "https://github.com/dtcitron/bioko_island_travel_materials",subdir = "macro.pfsi")` If you are on Windows you may need to first install [Rtools](https://cran.r-project.org/bin/windows/Rtools/), which provides support for packages with compiled code. More information is available at the [offical R documentation](https://cran.r-project.org/doc/manuals/r-devel/R-admin.html#Installing-packages). To install from source, you may open the .Rproj file within the **macro.pfsi** directory with RStudio and use the "Build and Reload" command. For more information, see the [RStudio documentation](https://support.rstudio.com/hc/en-us/articles/200486508-Building-Testing-and-Distributing-Packages).
