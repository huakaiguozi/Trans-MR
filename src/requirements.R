#-----------------renv包管理
# renv.activate()
# renv::status()
# renv::hydrate() #扫描并添加本地的包到renv管理项目下的包
# renv::snapshot() 
# 每当您添加、更新或删除包时，运行 renv::snapshot() 以保持 renv.lock 文件的最新状态。
# 这确保了团队成员和部署环境能够准确地恢复项目依赖。
#-----------------
#detach('package:ggGenshin',unload = TRUE)
library("ggplot2")
library("fixest")
#install.packages("TwoSampleMR", repos = c("https://mrcieu.r-universe.dev", "https://cloud.r-project.org"))
library('TwoSampleMR')
library('mr.raps')
library('MendelianRandomization')
library("penalized")
#library(sisVIVE)
library('withr')
library('tidyr')
library('nloptr')
library('ggpubr')
library('glmnet')
library('MASS')
library('FNN')
library('dplyr')
library('ks') #核密度估计包
library('car')
library('ridge')
library('np') #不确定是否用到了
library('RColorBrewer')

#remotes::install_github("Schwarzeneggerjune/Genshinpalette",force = T) #原神调色板
library('Genshinpalette')
#devtools::install_github("RestlessTail/ggGenshin")
library('ggGenshin')
#并行部分
library('foreach')
library('parallel')
library('doParallel')