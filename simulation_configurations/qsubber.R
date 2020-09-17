##
#
# Qsubber for ASTMH 2019 jobs
# Daniel T Citron
# 10/29/19
#
##
# Submit a set of jobs, where each job corresponds to a single run of the simulation
# 10/29/19 - We will, at first, submit 500 jobs, and then tomorrow make sure we can analyze them

# Baseline

for (i in 1:1000){
  sys.sub <- paste0("qsub -l m_mem_free=2G -l fthread=1 -P proj_mmc -q all.q",
                    " -N baseline_",i,
                    " -o /ihme/malaria_modeling/dtcitron/SpatialUncertainty/ASTMH19/baseline/log_files", # where output will go
                    " -e /ihme/malaria_modeling/dtcitron/SpatialUncertainty/ASTMH19/baseline/log_files"  # where errors will go
  )
  
  
  # This is a shell script, says which R to use, and passes arguments to the R script
  shell <- "/ihme/singularity-images/rstudio/shells/execRscript.sh" 
  # R script to execute
  script <- paste("-s", "/ihme/malaria_modeling/dtcitron/SpatialUncertainty/ASTMH19/baseline_job.R")
  
  args = c(i) # this is a seed, passed to the initial conditions of the simulation
  
  # Here's the text of the qsub call
  print(paste(sys.sub, shell, script, "\\", args))
  
  # And we can call the qsub as if it were from the command line here:
  system(paste(sys.sub, shell, script, "\\", args))
}

# Vaccination
for (i in 1:1000){
  sys.sub <- paste0("qsub -l m_mem_free=16G -l fthread=4 -P proj_mmc -q all.q",
                    " -N vax_",i,
                    " -o /ihme/malaria_modeling/dtcitron/SpatialUncertainty/ASTMH19/vaccine/log_files", # where output will go
                    " -e /ihme/malaria_modeling/dtcitron/SpatialUncertainty/ASTMH19/vaccine/log_files"  # where errors will go
  )
  
  
  # This is a shell script, says which R to use, and passes arguments to the R script
  shell <- "/ihme/singularity-images/rstudio/shells/execRscript.sh" 
  # R script to execute
  script <- paste("-s", "/ihme/malaria_modeling/dtcitron/SpatialUncertainty/ASTMH19/vaccine_job.R")
  
  args = c(i) # this is a seed, passed to the initial conditions of the simulation
  
  # Here's the text of the qsub call
  print(paste(sys.sub, shell, script, "\\", args))
  
  # And we can call the qsub as if it were from the command line here:
  system(paste(sys.sub, shell, script, "\\", args))
}

# Travel Fraction
for (i in 101:1000){
  sys.sub <- paste0("qsub -l m_mem_free=2G -l fthread=1 -P proj_mmc -q all.q",
                    " -N tf_",i,
                    " -o /ihme/malaria_modeling/dtcitron/SpatialUncertainty/ASTMH19/travelfrac/log_files", # where output will go
                    " -e /ihme/malaria_modeling/dtcitron/SpatialUncertainty/ASTMH19/travelfrac/log_files"  # where errors will go
  )
  
  
  # This is a shell script, says which R to use, and passes arguments to the R script
  shell <- "/ihme/singularity-images/rstudio/shells/execRscript.sh" 
  # R script to execute
  script <- paste("-s", "/ihme/malaria_modeling/dtcitron/SpatialUncertainty/ASTMH19/travelfrac_job.R")
  
  args = c(i) # this is a seed, passed to the initial conditions of the simulation
  
  # Here's the text of the qsub call
  print(paste(sys.sub, shell, script, "\\", args))
  
  # And we can call the qsub as if it were from the command line here:
  system(paste(sys.sub, shell, script, "\\", args))
}


# Local Residual Transmission - zero off-island transmission
for (i in 965){
  sys.sub <- paste0("qsub -l m_mem_free=2G -l fthread=1 -P proj_mmc -q all.q",
                    " -N localresidual_",i,
                    " -o /ihme/malaria_modeling/dtcitron/SpatialUncertainty/ASTMH19/localresid/log_files", # where output will go
                    " -e /ihme/malaria_modeling/dtcitron/SpatialUncertainty/ASTMH19/localresid/log_files"  # where errors will go
  )
  
  # This is a shell script, says which R to use, and passes arguments to the R script
  shell <- "/ihme/singularity-images/rstudio/shells/execRscript.sh" 
  # R script to execute
  script <- paste("-s", "/ihme/malaria_modeling/dtcitron/SpatialUncertainty/ASTMH19/localresid_job.R")
  
  args = c(i) # this is a seed, passed to the initial conditions of the simulation
  
  # Here's the text of the qsub call
  print(paste(sys.sub, shell, script, "\\", args))
  
  # And we can call the qsub as if it were from the command line here:
  system(paste(sys.sub, shell, script, "\\", args))
}

# Local Residual Transmission - zero off-island transmission - "true" ICs
for (i in 1:1000){
  sys.sub <- paste0("qsub -l m_mem_free=2G -l fthread=1 -P proj_mmc -q all.q",
                    " -N lr_trueICs_",i,
                    " -o /ihme/malaria_modeling/dtcitron/SpatialUncertainty/ASTMH19/localresid_trueICs/log_files", # where output will go
                    " -e /ihme/malaria_modeling/dtcitron/SpatialUncertainty/ASTMH19/localresid_trueICs/log_files"  # where errors will go
  )
  
  # This is a shell script, says which R to use, and passes arguments to the R script
  shell <- "/ihme/singularity-images/rstudio/shells/execRscript.sh" 
  # R script to execute
  script <- paste("-s", "/ihme/malaria_modeling/dtcitron/SpatialUncertainty/ASTMH19/localresid_trueICs_job.R")
  
  args = c(i) # this is a seed, passed to the initial conditions of the simulation
  
  # Here's the text of the qsub call
  print(paste(sys.sub, shell, script, "\\", args))
  
  # And we can call the qsub as if it were from the command line here:
  system(paste(sys.sub, shell, script, "\\", args))
}


# Low Importations - with "true" ICs
for (i in c(770, 834)){ #1:1000){
  sys.sub <- paste0("qsub -l m_mem_free=2G -l fthread=1 -P proj_mmc -q all.q",
                    " -N lowimport_trueICs_",i,
                    " -o /ihme/malaria_modeling/dtcitron/SpatialUncertainty/ASTMH19/lowimport_trueICs/log_files", # where output will go
                    " -e /ihme/malaria_modeling/dtcitron/SpatialUncertainty/ASTMH19/lowimport_trueICs/log_files"  # where errors will go
  )
  
  # This is a shell script, says which R to use, and passes arguments to the R script
  shell <- "/ihme/singularity-images/rstudio/shells/execRscript.sh" 
  # R script to execute
  script <- paste("-s", "/ihme/malaria_modeling/dtcitron/SpatialUncertainty/ASTMH19/lowimport_trueICs_job.R")
  
  args = c(i) # this is a seed, passed to the initial conditions of the simulation
  
  # Here's the text of the qsub call
  print(paste(sys.sub, shell, script, "\\", args))
  
  # And we can call the qsub as if it were from the command line here:
  system(paste(sys.sub, shell, script, "\\", args))
}




# Low Importation rate - 50% of off-island transmission
for (i in 1:100){
  sys.sub <- paste0("qsub -l m_mem_free=2G -l fthread=1 -P proj_mmc -q all.q",
                    " -N lowimport_",i,
                    " -o /ihme/malaria_modeling/dtcitron/SpatialUncertainty/ASTMH19/lowimport/log_files", # where output will go
                    " -e /ihme/malaria_modeling/dtcitron/SpatialUncertainty/ASTMH19/lowimport/log_files"  # where errors will go
  )
  
  # This is a shell script, says which R to use, and passes arguments to the R script
  shell <- "/ihme/singularity-images/rstudio/shells/execRscript.sh" 
  # R script to execute
  script <- paste("-s", "/ihme/malaria_modeling/dtcitron/SpatialUncertainty/ASTMH19/lowimport_job.R")
  
  args = c(i) # this is a seed, passed to the initial conditions of the simulation
  
  # Here's the text of the qsub call
  print(paste(sys.sub, shell, script, "\\", args))
  
  # And we can call the qsub as if it were from the command line here:
  system(paste(sys.sub, shell, script, "\\", args))
}
