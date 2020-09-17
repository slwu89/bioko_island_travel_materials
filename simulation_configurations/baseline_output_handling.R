"qsub -l m_mem_free=8G -l fthread=1 -P proj_mmc -q all.q -N baseline_output_handling -o /ihme/malaria_modeling/dtcitron/SpatialUncertainty/ASTMH19/baseline/log_files -e /ihme/malaria_modeling/dtcitron/SpatialUncertainty/ASTMH19/baseline/log_files /ihme/singularity-images/rstudio/shells/execRscript.sh -s /ihme/malaria_modeling/dtcitron/SpatialUncertainty/ASTMH19/baseline_output_handling.R"


library(data.table)
library(ggplot2)
library(Matrix)

setwd("/ihme/malaria_modeling/dtcitron")
library(here, lib.loc = "/ihme/malaria_modeling/dtcitron/Rlibs")


pop.data <- fread(here("SpatialUncertainty/data_clean/aggregated_2015_2018_travel_data.csv"))
areaId.list <- sort(pop.data[year == 2018]$areaId)
pop.dt <- data.table(patch = c(0:(241-1)), areaId = areaId.list)
pop.dt <- merge(pop.dt, pop.data[year == 2018, .(areaId, pop)], by = "areaId")


pfpr.data <- fread(here("SpatialUncertainty/data_clean/pfpr_draws.csv"))
pfpr.data <- merge(pfpr.data, pop.data[year == 2018, .(areaId)], by = "areaId", all = FALSE)


ensemble.file.list <- list.files(path = here("SpatialUncertainty/ASTMH19/baseline/sim_output"),
                                 pattern = "pfsi_[[:digit:]]+.csv")

# load in the first one
df_curr <- fread(here("SpatialUncertainty/ASTMH19/baseline/sim_output",ensemble.file.list[1]))
df_curr <- merge(pop.dt, df_curr, by = "patch")
df_curr[, s := (S_resident_home + S_resident_away)/pop, by = c("time" , "patch" , "time")]
df_curr[, i := (I_resident_home + I_resident_away)/pop, by = c("time" , "patch" , "time")]
df_curr[, p := (P_resident_home + P_resident_away)/pop, by = c("time" , "patch" , "time")]
# copy it to create the corresponding matrix:
mat_curr <- as.matrix(df_curr)
# then we use that matrix to create holders to catch data for the means and standard deviations:
mat_mean <- mat_curr
mat_sd <- mat_mean
mat_sd[, 5:ncol(mat_curr)] <- mat_mean[, 5:ncol(mat_curr)]^2
# now we loop over the other files in the list of ensemble outputs
nrun = length(ensemble.file.list)
for (i in 2:nrun){
  df_curr <- fread(here("SpatialUncertainty/ASTMH19/baseline/sim_output",ensemble.file.list[i]))
  
  df_curr <- merge(pop.dt, df_curr, by = "patch")
  df_curr[, s := (S_resident_home + S_resident_away)/pop, by = c("time" , "patch" , "time")]
  df_curr[, i := (I_resident_home + I_resident_away)/pop, by = c("time" , "patch" , "time")]
  df_curr[, p := (P_resident_home + P_resident_away)/pop, by = c("time" , "patch" , "time")]
  
  mat_curr <- as.matrix(df_curr)
  mat_mean[, 5:ncol(mat_curr)] <- mat_mean[, 5:ncol(mat_curr)] + mat_curr[, 5:ncol(mat_curr)]
  mat_sd[, 5:ncol(mat_curr)] <-  mat_sd[, 5:ncol(mat_curr)] + mat_curr[, 5:ncol(mat_curr)]^2
}
mat_mean <- as.data.table(mat_mean)
mat_sd <- as.data.table(mat_sd)

mat_mean[, 5:ncol(mat_curr)] <- mat_mean[, 5:ncol(mat_curr)]/nrun
mat_sd[, 5:ncol(mat_curr)] <- mat_sd[, 5:ncol(mat_curr)]/nrun - mat_mean[, 5:ncol(mat_curr)]^2
mat_sd[mat_sd < 0] <- 0
mat_sd[, 5:ncol(mat_curr)] <- sqrt(mat_sd[, 5:ncol(mat_curr)])

fwrite(mat_mean, here("SpatialUncertainty/ASTMH19/baseline/baseline_pr_means_1000.csv"))
fwrite(mat_sd, here("SpatialUncertainty/ASTMH19/baseline/baseline_pr_sds_1000.csv"))