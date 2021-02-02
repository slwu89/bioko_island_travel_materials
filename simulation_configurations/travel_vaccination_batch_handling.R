"qsub -l m_mem_free=8G -l fthread=1 -P proj_mmc -q all.q -N vax_output_handling -o /ihme/malaria_modeling/dtcitron/SpatialUncertainty/ASTMH19/vaccine/log_files -e /ihme/malaria_modeling/dtcitron/SpatialUncertainty/ASTMH19/vaccine/log_files /ihme/singularity-images/rstudio/shells/execRscript.sh -s /ihme/malaria_modeling/dtcitron/SpatialUncertainty/ASTMH19/vaccine_output_handling.R"


library(data.table)
library(ggplot2)
library(Matrix)
library(here)

pop.data <- fread(here("data/clean/aggregated_2015_2018_travel_data.csv"))
areaId.list <- sort(pop.data[year == 2018]$areaId)
pop.dt <- data.table(patch = c(0:(241-1)), areaId = areaId.list)
pop.dt <- merge(pop.dt, pop.data[year == 2018, .(areaId, pop)], by = "areaId")


pfpr.data <- fread(here("data/clean/pfpr_draws.csv"))
pfpr.data <- merge(pfpr.data, pop.data[year == 2018, .(areaId)], by = "areaId", all = FALSE)



ensemble.file.list <- list.files(path = here("data/simulation_outputs/travel_vaccination/"),
                                 pattern = "travel_vaccination_pfsi_[[:digit:]]+.csv")

# load in the first one
df_curr <- fread(here("data/simulation_outputs/travel_vaccination",ensemble.file.list[1]))
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
  df_curr <- fread(here("data/simulation_outputs/travel_vaccination",ensemble.file.list[i]))
  
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

fwrite(mat_mean, here("data/simulation_outputs/travel_vaccination/travel_vaccine_pr_means_1000.csv"))
fwrite(mat_sd, here("data/simulation_outputs/travel_vaccination/travel_vaccine_pr_sds_1000.csv"))




h <- melt(mat_mean,
          id.vars = c("time", "areaId"), 
          measure.vars = c("s","i","p"),
          value.name = "fraction")
h.sd <- melt(mat_sd,
             id.vars = c("time", "areaId"), 
             measure.vars = c("s","i","p"),
             value.name = "fraction.sd")
trav.vacc.dat <- merge(h, h.sd, by = c("time", "areaId", "variable"))

trav.vacc.dat[, lower.bound := max(fraction - fraction.sd, 0), by = c("areaId", "variable", "time")]



trav.vacc.plot.335 <- ggplot(
  data = trav.vacc.dat[areaId == 335 & variable == "i" & time < 6*365]) + 
  geom_errorbar(mapping = aes(x = time, 
                              ymin = lower.bound,
                              ymax = fraction + fraction.sd), 
                alpha = .5, color = "red") + 
  scale_x_continuous(name = "Time (Years)", 
                     breaks = c(365, 365*2, 365*3, 365*4, 365*5), labels = c(1:5),
                     expand = c(0,0)) + 
  scale_y_continuous(name = "Prevalence", 
                     expand = c(0,0),
                     limits = c(0,.105),
                     breaks = c(0,.05,.1)) + #ylab("Prevalence") +
  geom_point(mapping = aes(x = time, y = fraction), color = "black", shape = 16, size = 1) +
  theme(panel.background = element_rect(fill = "white", colour = "white"),
        axis.line = element_line(color = "black", size = 1),
        axis.text = element_text(size = 18),
        axis.title = element_text(size = 22),
        panel.grid = element_line(color = NULL),
        legend.position = c(.8, .8))

trav.vacc.plot.335

## Vaccination counting:
ensemble.file.list <- list.files(path = here("data/simulation_outputs/travel_vaccination/"),
                                 pattern = "travel_vaccination_vaxx_[[:digit:]]+.csv")
df <- fread(here("data/simulation_outputs/travel_vaccination/",ensemble.file.list[1]))
data_current <- df[,.("vaxx_events" = sum(vaxx_events), run.number = 1), by = "time"]

nrun = length(ensemble.file.list)
for (i in 2:nrun){
  df <- fread(here("data/simulation_outputs/travel_vaccination/",ensemble.file.list[i]))
  data_current <- rbind(data_current, df[,.("vaxx_events" = sum(vaxx_events), run.number = i), by = "time"])
}

data_current[,sum(vaxx_events)/7, by = "run.number"]
