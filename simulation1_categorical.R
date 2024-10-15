time_start <- Sys.time()
#---
setwd("C:/Users/Administrator/Desktop/WY/R_codes/Trans_MR_project")
ver <- 'v2'  #版本
date_mark <- format(Sys.time(), "%Y-%m-%d")
#===============================================================================结果存放位置
path_output_plots <- paste0('output/categorical/plots','/',date_mark,'-',ver)
path_output_tables <- paste0('output/categorical/tables','/',date_mark,'-',ver)

dir.create(path_output_plots,recursive = TRUE)
dir.create(path_output_tables,recursive = TRUE)
#===============================================================================执行顺序
source('src/requirements.R')
source('src/data_generating_functions_categorical.R')
source('src/computing_functions_categorical.R')
source('src/computation_categorical.R')
source('src/drawing_categorical.R')

#-----------------------
time_end <- Sys.time()
time_consume <- difftime(time_end,time_start)
time_consume 