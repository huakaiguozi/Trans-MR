time_start <- Sys.time()
#---
#setwd("C:/Users/Administrator/Desktop/WY/R_codes/Trans_MR_project")
ver <- 'v3-1000'  #版本
date_mark <- format(Sys.time(), "%Y-%m-%d")
#===============================================================================结果存放位置
path_output_plots <- paste0('output/simulation_misspecified/plots','/',date_mark,'-',ver)
path_output_tables <- paste0('output/simulation_misspecified/tables','/',date_mark,'-',ver)

dir.create(path_output_plots,recursive = TRUE)
dir.create(path_output_tables,recursive = TRUE)

ite_times <- 1000
#===============================================================================执行顺序
source('src/requirements.R')
source('src/data_generating_functions_misspecified.R')
source('src/computing_functions_categorical.R')
source('src/computation_categorical3.R')

save.image("TLMR数据/simulation1_misspecified.RData")

source('src/drawing_categorical.R')

#-----------------------
time_end <- Sys.time()
time_consume <- difftime(time_end,time_start)
time_consume 
