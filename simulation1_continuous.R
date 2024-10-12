time_start <- Sys.time()
#---
setwd("C:/Users/Administrator/Desktop/WY/R_codes/Trans_MR_project")
ver <- 1  #版本
date_mark <- format(Sys.time(), "%Y-%m-%d")
#===============================================================================结果存放位置
path_output_plots <- paste0(getwd(),'/output/continuous/plots')
path_output_tables <- paste0(getwd(),'/output/continuous/tables')
#===============================================================================执行顺序
source('/src/requirements.R')
source('/src/data_generating_functions_continuous.R')
source('/src/computing_functions_continuous.R')
source('/src/computation_continuous.R')
source('/src/drawing.R')

#-----------------------
time_end <- Sys.time()
time_consume <- difftime(time_end,time_start)
time_consume 