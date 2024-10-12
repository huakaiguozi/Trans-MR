time_start <- Sys.time()
#---
setwd("C:/Users/Administrator/Desktop/WY/R_codes/Trans_MR_project")
ver <- 'v1'  #版本
date_mark <- format(Sys.time(), "%Y-%m-%d")
#===============================================================================结果存放位置
path_output_plots <- 'output/categorical/plots'
path_output_tables <- 'output/categorical/tables'
#===============================================================================执行顺序
source('src/requirements.R')
source('src/data_generating_functions_categorical.R')
source('src/computing_functions_categorical.R')
source('src/computation_categorical.R')
source('src/drawing.R')

#-----------------------
time_end <- Sys.time()
time_consume <- difftime(time_end,time_start)
time_consume 