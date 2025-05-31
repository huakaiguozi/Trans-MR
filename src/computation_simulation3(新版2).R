#time_start <- Sys.time()
# 设置并行计算的核心数量
cores <- detectCores()
registerDoParallel(cores = cores-10)
#-----------------------------------------------设定
# n_0_list <- c(1000,2000,3000,4000,5000)  #样本量
# #n_0_list <- c(1000)  #样本量
# n_1_list <- c(400)
# gamma_list <- c(0,0.2,0.4,0.6,0.8,1.0) #gamma大小的调节参数
# heter_level_list <- c(0,0.1,0.2,0.3,0.4,0.5)
# p_Z_list <- c(3,5,10,15,20)
# inter_strength_list <- c(0,0.2,0.4,0.6,0.8,1.0)
# alpha_list <- c(0,0.2,0.4,0.6,0.8,1.0)
# eta_list <- c(0,0.2,0.4,0.6,0.8,1.0)
#---
# n_0_list <- c(2000,4000,6000,8000,10000,12000)  #样本量
n_0_list <- c(4000,8000,12000,16000,20000)  #样本量
#n_0_list2 <- c(500,800,1200,1600,2000,3000)
n_0_list2 <- c(500,1500,3000)
#n_0_list <- c(1000)  #样本量
n_1_list <- c(400)
gamma_list <- c(0,0.2,0.4,0.6,0.8,1) #gamma大小的调节参数
heter_level_list <- c(0,0.2,0.4,0.6,0.8,1)
#-----
# p_Z_list <- c(3,4,6,8,10,12)
p_Z_list <- c(3,4)
#-----
inter_strength_list <- c(0,0.2,0.4,0.6,0.8,1)
alpha_list <- c(0,0.2,0.4,0.6,0.8,1)
# beta_list <- c(1,2)
# beta_list <- list(c(1,1),c(10.6,14.8))
omegam_list <- c(1,2)
# #单独循环p_Z
# p_Z2_list <- c(3,6,9,12,15)

#似乎可以不遍历这个
eta_list <- c(1.0)
#beta_setting_list <- c(0,0.2,0.4,0.6,0.8,1.0)
# 暂时没用到的遍历参数：交互比例参数、水平多效性比例参数
ratio_inter_list <- c(1.0)
ratio_pleiotropy <- c(1.0)
#ratio_inter_list <- c(0,0.2,0.4,0.6,0.8,1.0)

# n_0_list <- c(1000)
# n_1_list <- c(400)
# gamma_list <- c(0.3)
# heter_level_list <- c(0.3)
# ratio_inter_list <- c(0.6)

#----Twosamplemr包方法

# 共11种方法
other_methods <- c()
my_methods <- c('TLMR-OR'
                ,'TLMR-IPW'
                ,'TLMR-AIPW')
method_names <- c(other_methods,my_methods)
len_methods <- length(method_names)
len_other_methods <- length(other_methods)
len_my_methods <- length(my_methods)
#path_file <- 'C:/Users/Administrator/Desktop/WY/R_codes/TSMR_simulation/241005'
#C:\Users\Administrator\Desktop\WY\R_codes\TSMR_simulation\241002
#-----------------------------------------------参数设定
#=======
# gamma_list <- c(0.2,0.4,0.6,0.8,1.0)
# heter_level_list <- c(0.2,0.4,0.6,0.8,1.0)
#=======
len_gamma <- length(gamma_list)
len_heter <- length(heter_level_list)
len_inter <- length(ratio_inter_list)
len_n0 <- length(n_0_list)
len_n02 <- length(n_0_list2)
len_n1 <- length(n_1_list)
len_pz <- length(p_Z_list)
len_inter_strength <- length(inter_strength_list)
len_alpha <- length(alpha_list)
len_beta <- length(omegam_list)
len_eta <- length(eta_list)

len_ratio_pleiotropy <- length(ratio_pleiotropy)
#参数空间
# beta=beta,eta_2=eta_2,gamma=gamma,ratio_pleiotropy=ratio_pleiotropy,heter_level=heter_level
pi_data1 <- matrix(round(with_seed(1, runif(90000,0,2)),2),nrow = 300)
pi_data2 <- matrix(round(with_seed(1, runif(90000,0,2)),2),nrow = 300)
alpha_data <- round(with_seed(1,runif(1000,0.2,0.8)),2)
gamma_data <- round(with_seed(1,runif(1000,0.1,0.2)),2)
eta1_data <- round(with_seed(1,runif(1000,-0.5,0.5)),2)
eta2_data <- round(with_seed(1,runif(1000,-0.5,0.5)),2)
omega_data <- round(with_seed(1, runif(90000,0,2)),2)
#--------------------------------名称排序
# c("IVW","Egger" ,"Weighted median" ,      "Weighted mode"  ,       "MR-RAPS"      ,
#   "Contamination Mixture" ,"MR-Robust"   ,          "MR-Lasso"  ,            "Trans-OR"       ,       "Trans-IPW"  ,
#   "Trans-AIPW")
# 
# c("IVW","Egger" ,"Weighted median" ,      "Weighted mode"  ,       "MR-RAPS"      ,
# "MR-conmix" ,"MR-Robust"   ,          "MR-Lasso"  ,            "TLMR-OR"       ,       "TLMR-IPW"  ,
# "TLMR-AIPW")
# 
# c("IVW","Egger" ,"Weighted median" ,      "Weighted mode"  ,       "MR-RAPS"      ,
#   "Contamination Mixture" ,"MR-conmix","MR-Robust"   ,          "MR-Lasso"  ,            "Trans-OR"       , "TLMR-OR",      "Trans-IPW"  ,"TLMR-IPW"  ,
#   "Trans-AIPW","TLMR-AIPW")
method_names_ordered <- factor(method_names,levels = method_names)
index_ordered <- factor(c('Bias','SD','MSE'),levels = c('Bias','SD','MSE'))
#--------------------------------结果保存
save_df <- data.frame()  #结果保存
col_names <- c('iteration time','omegam','p_Z', 'gamma','heterogeneity','n_0','n_1','inter_ratio','inter_strength','alpha_modifier','eta_modifier','method'
               ,paste0('omega',c(1:length(omega)))
               ,paste0('p_power',c(1:length(omega)))
               ,paste0('p_cover',c(1:length(omega))))

save_df2 <- data.frame()  #结果保存
col_names2 <- c('max iteration times','omegam','p_Z', 'gamma','heterogeneity','n_0','n_1','inter_ratio','inter_strength','alpha_modifier','eta_modifier','method','index',paste0('num',c(1:length(omega))))

save_power_list <- vector('list',length = len_n02)
save_power_list2 <- vector('list',length = len_n02)
#----------
zongxunhuancishu<-0
k_min <- 1

k_pz <- 1
p_Z_fix <- p_Z_list[k_pz]
k_gamma <- 1
k_gamma2 <-1 #length(gamma_list)
k_alpha <- len_alpha
k_alpha2 <- len_alpha
k_heter <- len_heter
k_n0 <- len_n0
k_n02 <- 1
k_inter <- len_inter
k_inter_strength <- len_inter_strength
k_eta <- len_eta
k_beta <- 1
k_beta2 <- 2
#k_ratio_pl <- len_ratio_pleiotropy
all_params_list <- vector('list') 
for(i_n0 in 1:len_n0){for (i_n1 in 1:len_n1){for(i_gamma in 1:len_gamma){for (i_heter in 1:len_heter){for(i_inter in 1:len_inter){for (i_pz in 1:len_pz) {for (i_inter_stren in 1:len_inter_strength) {for (i_alpha in 1:len_alpha) {for (i_eta in 1:len_eta) {for (i_beta in 1:len_beta) {for (i_n02 in 1:len_n02) {
  omegam <- omegam_list[i_beta]
  n_02 <- n_0_list2[i_n02]
  n_0 <- n_0_list[i_n0]
  
  heter_level <- heter_level_list[i_heter]
  ratio_inter <- ratio_inter_list[i_inter]
  gamma_multi <- gamma_list[i_gamma]
  inter_strength <- inter_strength_list[i_inter_stren]
  p_Z <- p_Z_list[i_pz]
  alpha_modifier <- alpha_list[i_alpha]
  eta_modifier <- eta_list[i_eta]
  param_list <- list(n_0=n_0
                     ,heter_level=heter_level
                     ,ratio_inter=ratio_inter
                     ,gamma_multi=gamma_multi
                     ,inter_strength=inter_strength
                     ,p_Z=p_Z
                     ,alpha_modifier=alpha_modifier
                     ,eta_modifier=eta_modifier
  )
  #gamma大小
  condition1 <- (heter_level == heter_level_list[k_heter] & n_0 == n_0_list[k_n0] & ratio_inter == ratio_inter_list[k_inter] & p_Z==p_Z_list[k_pz] & inter_strength==inter_strength_list[k_inter_strength] & alpha_modifier == alpha_list[k_alpha] & eta_modifier == eta_list[k_eta]& omegam == omegam_list[k_beta]& n_02 == n_0_list2[k_n02])
  # 异质性水平
  condition2 <- (gamma_multi == gamma_list[k_gamma] & n_0 == n_0_list[k_n0] & ratio_inter == ratio_inter_list[k_inter] & p_Z==p_Z_list[k_pz] & inter_strength==inter_strength_list[k_inter_strength] & alpha_modifier == alpha_list[k_alpha] & eta_modifier == eta_list[k_eta]& omegam == omegam_list[k_beta]& n_02 == n_0_list2[k_n02])
  # 样本量n0
  condition3 <- (heter_level == heter_level_list[k_heter] & gamma_multi == gamma_list[k_gamma]  & ratio_inter == ratio_inter_list[k_inter] & p_Z==p_Z_list[k_pz] & inter_strength==inter_strength_list[k_inter_strength] & alpha_modifier == alpha_list[k_alpha] & eta_modifier == eta_list[k_eta]& omegam == omegam_list[k_beta]& n_02 == n_0_list2[k_n02])
  # 修饰比率
  condition4 <- (heter_level == heter_level_list[k_heter] & gamma_multi == gamma_list[k_gamma] & n_0 == n_0_list[k_n0] & p_Z==p_Z_list[k_pz] & inter_strength==inter_strength_list[k_inter_strength] & alpha_modifier == alpha_list[k_alpha] & eta_modifier == eta_list[k_eta]& omegam == omegam_list[k_beta]& n_02 == n_0_list2[k_n02])
  # 工具个数：p_Z
  condition5 <- (heter_level == heter_level_list[k_heter] & gamma_multi == gamma_list[k_gamma] & n_0 == n_0_list[k_n0] & ratio_inter == ratio_inter_list[k_inter] & inter_strength==inter_strength_list[k_inter_strength] & alpha_modifier == alpha_list[k_alpha] & eta_modifier == eta_list[k_eta]& omegam == omegam_list[k_beta]& n_02 == n_0_list2[k_n02])
  # 修饰强度
  condition6 <- (heter_level == heter_level_list[k_heter] & gamma_multi == gamma_list[k_gamma2] & n_0 == n_0_list[k_n0] & ratio_inter == ratio_inter_list[k_inter] & p_Z==p_Z_list[k_pz] & alpha_modifier == alpha_list[k_alpha2] & eta_modifier == eta_list[k_eta]& omegam == omegam_list[k_beta]& n_02 == n_0_list2[k_n02])
  # eta
  condition7 <- (heter_level == heter_level_list[k_heter] & gamma_multi == gamma_list[k_gamma] & n_0 == n_0_list[k_n0] & ratio_inter == ratio_inter_list[k_inter] & p_Z==p_Z_list[k_pz] & inter_strength==inter_strength_list[k_inter_strength] & alpha_modifier == alpha_list[k_alpha] & omegam == omegam_list[k_beta]& n_02 == n_0_list2[k_n02])
  # alpha
  condition8 <- (heter_level == heter_level_list[k_heter] & gamma_multi == gamma_list[k_gamma] & n_0 == n_0_list[k_n0] & ratio_inter == ratio_inter_list[k_inter] & p_Z==p_Z_list[k_pz] & inter_strength==inter_strength_list[k_inter_strength] & eta_modifier == eta_list[k_eta]& omegam == omegam_list[k_beta] & n_02 == n_0_list2[k_n02] )
  # 另一个beta下遍历样本量
  condition9 <- (heter_level == heter_level_list[k_heter] & gamma_multi == gamma_list[k_gamma] & n_0 == n_0_list[k_n0] & ratio_inter == ratio_inter_list[k_inter] & p_Z==p_Z_list[k_pz] & inter_strength==inter_strength_list[k_inter_strength] & eta_modifier == eta_list[k_eta] & alpha_modifier == alpha_list[k_alpha] & omegam == omegam_list[k_beta2])
  # 遍历样本量、效应值
  condition10 <- (heter_level == heter_level_list[k_heter] & gamma_multi == gamma_list[k_gamma]  & ratio_inter == ratio_inter_list[k_inter] & p_Z==p_Z_list[k_pz] & inter_strength==inter_strength_list[k_inter_strength] & eta_modifier == eta_list[k_eta] & alpha_modifier == alpha_list[k_alpha]& n_02 == n_0_list2[k_n02])
  
  if(condition1|condition2|condition3|condition4|condition5|condition6|condition7|condition8|condition9|condition10){
    zongxunhuancishu<-zongxunhuancishu+1
    all_params_list[[zongxunhuancishu]] <- param_list
  }
  
}}}}}}}}}}}

#BL_times <- 0
xunhuanshijiweizhi <- 0
jishi1 <- Sys.time()
for(i_n0 in 1:len_n0){
  for (i_n1 in 1:len_n1){
    for(i_gamma in 1:len_gamma){
      for (i_heter in 1:len_heter){
        for(i_inter in 1:len_inter){
          for (i_pz in 1:len_pz) {
            for (i_inter_stren in 1:len_inter_strength) {
              for (i_alpha in 1:len_alpha) {
                for (i_eta in 1:len_eta) {
                  for (i_beta in 1:len_beta) {
                    for (i_n02 in 1:len_n02) {
                      ##---计数
                      #BL_times <- BL_times + 1
                      #print(paste('第',as.character(BL_times),'次'))
                      #---
                      # i_n0 <- 1
                      # i_n1 <- 1
                      # i_gamma <- 1
                      # i_heter <- 1
                      # i_inter <- 1
                      # i_pz <- 1
                      # i_inter_stren <- 1
                      # i_alpha <- 1
                      # i_eta <- 1
                      # ite_k <- 1
                      #===================================================================参数设置
                      beta <- 0.4
                      
                      n_0 <- n_0_list[i_n0]
                      n_02 <- n_0_list2[i_n02]
                      n_1 <- n_0
                      
                      heter_level <- heter_level_list[i_heter]
                      ratio_inter <- ratio_inter_list[i_inter]
                      gamma_multi <- gamma_list[i_gamma]
                      
                      inter_strength <- inter_strength_list[i_inter_stren]
                      alpha_modifier <- alpha_list[i_alpha]
                      eta_modifier <- eta_list[i_eta]
                      
                      #--------------------------
                      p_Z <- p_Z_list[i_pz]
                      # if(p_Z == p_Z_list[k_pz]){  #下面改成所有工具个数下都具有相同的强度
                      #   p_Z_max = p_Z
                      # }else{
                      #   p_Z_max = p_Z_list[len_pz]
                      # }
                      p_Z_max =p_Z
                      #-----------
                      # if(p_Z > p_Z_list[k_pz]){
                      #   n_mag <- p_Z/p_Z_list[k_pz]
                      #   alpha_modifier <- alpha_modifier * n_mag
                      #   inter_strength <- inter_strength * n_mag
                      # }
                      #U_md <- p_Z/p_Z_list[k_pz]
                      #-----------
                      p_C <- 2
                      omegam <- omegam_list[i_beta]
                      omega <- omega_data[1:p_C]
                      if(omegam == 2){
                        omega <- c(0.05,0.05)
                      }
                      len_omega <- length(omega)
                      eta_1 <- eta1_data[1:p_C] * eta_modifier
                      eta_2 <- eta2_data[1:p_C] * eta_modifier
                      
                      #-----
                      #alpha <- seq(2,4,2/(p_Z-1))*alpha_modifier
                      alpha <- alpha_data[1:p_Z_max] * alpha_modifier
                      pi <- pi_data1[1:p_Z_max,1:p_C] * inter_strength
                      pi2 <- pi_data2[1:p_Z_max,1:p_C] * inter_strength
                      
                      gamma <- gamma_data[1:p_Z_max] * gamma_multi
                      indice_gamma1 <- seq(1,p_Z_max,3)
                      num_gamma1 <- round(p_Z_max/3)
                      indice_gamma1 <- indice_gamma1[1:num_gamma1]
                      indice_gamma0 <- setdiff(1:p_Z_max,indice_gamma1)
                      
                      gamma[indice_gamma0] <- 0
                      # if(p_Z_max<10){
                      #   gamma[c(2:length(gamma))] <- 0
                      # }else{
                      #   gamma[c(1,2,4,5,7,8,10,11,13,14)] <- 0
                      # }
                      
                      
                      alpha_true <- alpha_data[1:p_Z] * alpha_modifier
                      pi_true <- pi_data1[1:p_Z,1:p_C] * inter_strength
                      pi2_true <- pi_data2[1:p_Z,1:p_C] * inter_strength
                      gamma_true <- gamma_data[1:p_Z] * gamma_multi
                      #----
                      param_X <- list(alpha=alpha,pi=pi,pi2=pi2,eta_1=eta_1
                                      #,U_md = U_md
                      )
                      param_Y <- list(beta=beta,eta_2=eta_2,gamma=gamma,ratio_pleiotropy=ratio_pleiotropy
                                      ,heter_level=heter_level,omega=omega
                                      #,U_md = U_md
                      )
                      #heritability <- 0.7
                      pval_i <- 0.05
                      #0.05/p_Z
                      #-------------------
                      condition1 <- (heter_level == heter_level_list[k_heter] & n_0 == n_0_list[k_n0] & ratio_inter == ratio_inter_list[k_inter] & p_Z==p_Z_list[k_pz] & inter_strength==inter_strength_list[k_inter_strength] & alpha_modifier == alpha_list[k_alpha] & eta_modifier == eta_list[k_eta] & omegam == omegam_list[k_beta]& n_02 == n_0_list2[k_n02])
                      condition2 <- (gamma_multi == gamma_list[k_gamma] & n_0 == n_0_list[k_n0] & ratio_inter == ratio_inter_list[k_inter] & p_Z==p_Z_list[k_pz] & inter_strength==inter_strength_list[k_inter_strength] & alpha_modifier == alpha_list[k_alpha] & eta_modifier == eta_list[k_eta] & omegam == omegam_list[k_beta]& n_02 == n_0_list2[k_n02])
                      condition3 <- (heter_level == heter_level_list[k_heter] & gamma_multi == gamma_list[k_gamma]  & ratio_inter == ratio_inter_list[k_inter] & p_Z==p_Z_list[k_pz] & inter_strength==inter_strength_list[k_inter_strength] & alpha_modifier == alpha_list[k_alpha] & eta_modifier == eta_list[k_eta] & omegam == omegam_list[k_beta]& n_02 == n_0_list2[k_n02])
                      condition4 <- (heter_level == heter_level_list[k_heter] & gamma_multi == gamma_list[k_gamma] & n_0 == n_0_list[k_n0] & p_Z==p_Z_list[k_pz] & inter_strength==inter_strength_list[k_inter_strength] & alpha_modifier == alpha_list[k_alpha] & eta_modifier == eta_list[k_eta] & omegam == omegam_list[k_beta]& n_02 == n_0_list2[k_n02])
                      condition5 <- (heter_level == heter_level_list[k_heter] & gamma_multi == gamma_list[k_gamma] & n_0 == n_0_list[k_n0] & ratio_inter == ratio_inter_list[k_inter] & inter_strength==inter_strength_list[k_inter_strength] & alpha_modifier == alpha_list[k_alpha] & eta_modifier == eta_list[k_eta]& omegam == omegam_list[k_beta]& n_02 == n_0_list2[k_n02])
                      condition6 <- (heter_level == heter_level_list[k_heter] & gamma_multi == gamma_list[k_gamma2] & n_0 == n_0_list[k_n0] & ratio_inter == ratio_inter_list[k_inter] & p_Z==p_Z_list[k_pz] & alpha_modifier == alpha_list[k_alpha2] & eta_modifier == eta_list[k_eta]& omegam == omegam_list[k_beta]& n_02 == n_0_list2[k_n02])
                      condition7 <- (heter_level == heter_level_list[k_heter] & gamma_multi == gamma_list[k_gamma] & n_0 == n_0_list[k_n0] & ratio_inter == ratio_inter_list[k_inter] & p_Z==p_Z_list[k_pz] & inter_strength==inter_strength_list[k_inter_strength] & alpha_modifier == alpha_list[k_alpha] & omegam == omegam_list[k_beta]& n_02 == n_0_list2[k_n02])
                      condition8 <- (heter_level == heter_level_list[k_heter] & gamma_multi == gamma_list[k_gamma] & n_0 == n_0_list[k_n0] & ratio_inter == ratio_inter_list[k_inter] & p_Z==p_Z_list[k_pz] & inter_strength==inter_strength_list[k_inter_strength] & eta_modifier == eta_list[k_eta] & omegam == omegam_list[k_beta]& n_02 == n_0_list2[k_n02])
                      condition9 <- (heter_level == heter_level_list[k_heter] & gamma_multi == gamma_list[k_gamma] & n_0 == n_0_list[k_n0] & ratio_inter == ratio_inter_list[k_inter] & p_Z==p_Z_list[k_pz] & inter_strength==inter_strength_list[k_inter_strength] & eta_modifier == eta_list[k_eta] & alpha_modifier == alpha_list[k_alpha] & omegam == omegam_list[k_beta2])
                      condition10 <- (heter_level == heter_level_list[k_heter] & gamma_multi == gamma_list[k_gamma]  & ratio_inter == ratio_inter_list[k_inter] & p_Z==p_Z_list[k_pz] & inter_strength==inter_strength_list[k_inter_strength] & eta_modifier == eta_list[k_eta] & alpha_modifier == alpha_list[k_alpha]& n_02 == n_0_list2[k_n02])
                      
                      if  ( !(condition1|condition2|condition3|condition4|condition5|condition6|condition7|condition8|condition9|condition10) ){
                        next
                      }else{
                        xunhuanshijiweizhi <- xunhuanshijiweizhi+1
                        jishi2 <- Sys.time()
                        xunhuanyongshi <- difftime(jishi2,jishi1,units = 'hours')
                        jishi1 <- jishi2
                        print(paste('当前循环位置：',(xunhuanshijiweizhi-1),'/',zongxunhuancishu,'，','本次循环用时：',as.character(round(xunhuanyongshi,5)),'小时，当前时间：',as.character(Sys.time())))
                      }
                      #----------------------
                      if(condition9){
                        n_0 <- n_02
                      }
                      #===========================================================================
                      ## 结果保存
                      param_df2 <- data.frame(matrix(rep(c(ite_times,omegam,p_Z,gamma_multi,heter_level,n_0,n_1,ratio_inter,inter_strength,alpha_modifier,eta_modifier),len_methods*3),nrow = len_methods*3,byrow = TRUE))
                      param_df2 <- cbind(param_df2,rep(method_names_ordered,3))
                      param_df2 <- cbind(param_df2,rep(index_ordered,each=len_methods))
                      
                      
                      param_df <- data.frame(matrix(rep(c(omegam,p_Z,gamma_multi,heter_level,n_0,n_1,ratio_inter,inter_strength,alpha_modifier,eta_modifier),len_methods*ite_times),nrow = len_methods*ite_times,byrow = TRUE))
                      param_df <- cbind(rep(c(1:ite_times),each=len_methods),param_df)
                      param_df <- cbind(param_df,rep(method_names_ordered,ite_times))
                      #----------------------------
                      beta_hat_ite <- foreach(ite_k=1:ite_times, .export =c('heter_level')
                                              , .packages = c('fixest','TwoSampleMR','glmnet','mr.raps','penalized','ridge','FNN','data.table'))%dopar%
                        {
                          #for(ite_k in 1:ite_times){
                          tryCatch({
                            #indicator_print <- 1
                            
                            seed <- ite_k * 100
                            set.seed(seed)
                            #-----------------------------------------------生成随机数
                            data_0 <- data_generate_Y(n=n_0,p_Z=p_Z,p_Z_max = p_Z_max
                                                      ,p_C=p_C,param_X=param_X,param_Y=param_Y)
                            Z_0 <- data_0$Z
                            C_0 <- data_0$C
                            X_0 <- data_0$X
                            Y_0 <- data_0$Y
                            #------------------------------------样本1
                            data_1 <- data_generate_X(n=n_1,p_Z=p_Z,p_Z_max = p_Z_max
                                                      ,p_C=p_C,param_X=param_X)
                            #alpha_true <- data_1$alpha
                            Z_1 <- data_1$Z
                            C_1 <- data_1$C
                            X_1 <- data_1$X
                            #----------------------------
                            Z_0_save <- Z_0
                            Z_1_save <- Z_1
                            
                            Z_0_sum <- matrix(apply(Z_0, 1, sum),ncol = 1)
                            Z_1_sum <- matrix(apply(Z_1, 1, sum),ncol = 1)
                            Z_0 <- Z_0_sum
                            Z_1 <- Z_1_sum
                            # Z_0 <- cbind(Z_0,Z_0_sum)
                            # Z_1 <- cbind(Z_1,Z_1_sum)
                            #----------------------------
                            conditions_met <- FALSE
                            while(!conditions_met){
                              #condition_list <- rep(FALSE,p_Z)
                              condition_list <- rep(FALSE,1)
                              # Z_0_sum <- matrix(apply(Z_0, 1, sum),ncol = 1)
                              # Z_1_sum <- matrix(apply(Z_1, 1, sum),ncol = 1)
                              # Z_0 <- cbind(Z_0,Z_0_sum)
                              # Z_1 <- cbind(Z_1,Z_1_sum)
                              for (i in 1:1 ){
                                P_0 <- cbind(Z_0[,i], C_0)
                                P_1 <- cbind(Z_1[,i], C_1)
                                P_0_strings <- apply(P_0, 1, paste, collapse = ",")
                                P_1_strings <- apply(P_1, 1, paste, collapse = ",")
                                
                                if( all(unique(P_0_strings) %in% unique(P_1_strings))) {
                                  condition_list[i] <- TRUE
                                }else{break}
                              }
                              if(all(condition_list)){
                                conditions_met <- TRUE
                              }else{
                                # 生成新的随机数据
                                data_0 <- data_generate_Y(n=n_0,p_Z=p_Z,p_Z_max = p_Z_max
                                                          ,p_C=p_C,param_X=param_X,param_Y=param_Y)
                                Z_0 <- data_0$Z
                                C_0 <- data_0$C
                                X_0 <- data_0$X
                                Y_0 <- data_0$Y
                                # 生成新的随机数据
                                data_1 <- data_generate_X(n=n_1,p_Z=p_Z,p_Z_max = p_Z_max
                                                          ,p_C=p_C,param_X=param_X)
                                #alpha_true <- data_1$alpha
                                Z_1 <- data_1$Z
                                C_1 <- data_1$C
                                X_1 <- data_1$X
                                #----------------------------
                                Z_0_save <- Z_0
                                Z_1_save <- Z_1
                                
                                Z_0_sum <- matrix(apply(Z_0, 1, sum),ncol = 1)
                                Z_1_sum <- matrix(apply(Z_1, 1, sum),ncol = 1)
                                Z_0 <- Z_0_sum
                                Z_1 <- Z_1_sum
                              }
                            }
                            #----------------------------
                            #cor(X_hat_matrix_IPW,vec(X_hat_matrix_IPW) * C_0)
                            # cor(X_hat_matrix_IPW, Z_0 * C_0)
                            # diag(solve(cor(model.matrix(model3)[, -1])))
                            # # 提取模型矩阵
                            # X <- model.matrix(model3)[, -1] # 去掉截距项
                            # kappa(X, exact = TRUE)
                            # 
                            # # 计算VIF
                            # vif_values <- diag(solve(cor(X)))  # 对角线元素即为VIF
                            # print(vif_values)
                            # 
                            # kappa(X, exact = TRUE)
                            
                            #===============================================================Trans-MR
                            
                            #===============================================================迁移估计量
                            #Trans_OR
                            X_hat_matrix_OR <- Trans_OLS2(Z_0=Z_0,Z_1=Z_1,C_0=C_0,C_1=C_1,X_1=X_1)$X_hat_matrix
                            #Trans_IPW
                            X_hat_matrix_IPW <- IPW2(Z_0=Z_0,Z_1=Z_1,C_0=C_0,C_1=C_1,X_1=X_1)$X_hat_matrix
                            #Trans_AIPW
                            X_hat_matrix_AIPW <- AIPW2(Z_0=Z_0,Z_1=Z_1,C_0=C_0,C_1=C_1,X_1=X_1
                            )$X_hat_matrix
                            #---------------------------------------------------------------Trans_OR0         21
                            #cor(X_hat_matrix_IPW,vec(X_hat_matrix_IPW) * C_0)
                            cor(X_hat_matrix_IPW, Z_0 * C_0)
                            diag(solve(cor(model.matrix(model3)[, -1])))
                            
                            O_x <- cbind(X_hat_matrix_OR,C_0,Z_0,as.vector(X_hat_matrix_OR) * C_0)
                            penalty_factor <- c(0,0,0,0,1,1)
                            #lambda_custum <- seq(0,0.1, length.out=101)
                            
                            cv_model <- cv.glmnet(O_x, Y_0, alpha = 0,nfolds = 5, lambda.min.ratio = 1e-5#, lambda = lambda_custum
                                                  , penalty.factor=penalty_factor
                                                  , intercept = T)
                            best_lambda <- cv_model$lambda.min
                            best_model <- glmnet(O_x, Y_0, alpha = 0
                                                 , penalty.factor=penalty_factor
                                                 , lambda = best_lambda,intercept = T)
                            
                            n_coef <- length(best_model$beta)
                            omega_hat2 <- best_model$beta[(n_coef-len_omega+1):n_coef]
                            
                            #sd_omega_hat2 <- summary(model2)$coefficients[,2][(n_coef-len_omega+1):n_coef]
                            
                            model2 <- lm(Y_0~X_hat_matrix_OR + C_0+ Z_0+ X_hat_matrix_OR * C_0 )  #有Z
                            cor(X_hat_matrix_OR, vec(X_hat_matrix_OR) * C_0)
                            diag(solve(cor(model.matrix(model2)[, -1])))
                            kappa(model.matrix(model2)[, -1], exact = TRUE)
                            
                            n_coef <- length(model2$coefficients)
                            omega_hat2 <- model2$coefficients[(n_coef-len_omega+1):n_coef]
                            sd_omega_hat2 <- summary(model2)$coefficients[,2][(n_coef-len_omega+1):n_coef]
                            #---------------------------------------------------------------Trans_IPW0          31
                            O_x <- cbind(X_hat_matrix_IPW,C_0,Z_0,as.vector(X_hat_matrix_IPW) * C_0)
                            cv_model <- cv.glmnet(O_x, Y_0, alpha = 0,nfolds = 5, lambda.min.ratio = 1e-5#, lambda = lambda_custum
                                                  , penalty.factor=penalty_factor
                                                  , intercept = T)
                            best_lambda <- cv_model$lambda.min
                            best_model <- glmnet(O_x, Y_0, alpha = 0
                                                 , penalty.factor=penalty_factor
                                                 , lambda = best_lambda,intercept = T)
                            n_coef <- length(best_model$beta)
                            omega_hat3 <- best_model$beta[(n_coef-len_omega+1):n_coef]
                            
                            #---
                            model3_1 <- lm(X_hat_matrix_IPW~Z_0*C_0)
                            coef_pi <- matrix(model3_1$coefficients[(length(model3_1$coefficients)-p_C+1):length(model3_1$coefficients)],ncol = p_C)
                            X_hat_matrix_IPW_2 <- X_hat_matrix_IPW - rowSums((Z_0 %*% coef_pi) * C_0)
                            #rowSums((Z_0 %*% coef_pi) * C_0)
                            
                            model3 <- lm(Y_0~X_hat_matrix_IPW + C_0+ Z_0+ X_hat_matrix_IPW * C_0)
                            cor(X_hat_matrix_IPW, vec(X_hat_matrix_IPW) * C_0)
                            diag(solve(cor(model.matrix(model3)[, -1])))
                            kappa(model.matrix(model3)[, -1], exact = TRUE)
                            
                            n_coef <- length(model3$coefficients)
                            omega_hat3 <- model3$coefficients[(n_coef-len_omega+1):n_coef]
                            sd_omega_hat3 <- summary(model3)$coefficients[,2][(n_coef-len_omega+1):n_coef]
                            #---------------------------------------------------------------Trans_AIPW0          41
                            O_x <- cbind(X_hat_matrix_AIPW,C_0,Z_0,as.vector(X_hat_matrix_AIPW) * C_0)
                            cv_model <- cv.glmnet(O_x, Y_0, alpha = 0,nfolds = 5, lambda.min.ratio = 1e-5#, lambda = lambda_custum
                                                  , penalty.factor=penalty_factor
                                                  , intercept = T)
                            best_lambda <- cv_model$lambda.min
                            best_model <- glmnet(O_x, Y_0, alpha = 0
                                                 , penalty.factor=penalty_factor
                                                 , lambda = best_lambda,intercept = T)
                            omega_hat4 <- best_model$beta[(n_coef-len_omega+1):n_coef]
                            
                            # model4 <- lm(Y_0~X_hat_matrix_AIPW + C_0+ Z_0+ X_hat_matrix_AIPW * C_0)
                            # omega_hat4 <- model4$coefficients[(n_coef-len_omega+1):n_coef]
                            # sd_omega_hat4 <- summary(model4)$coefficients[,2][(n_coef-len_omega+1):n_coef]
                            #----------------------------结果保存
                            #beta_hat_vec <- c(beta_hat_1, beta_hat_2, beta_hat_3,beta_hat_4)
                            omega_hat_list <- list(omega_hat2,omega_hat3,omega_hat4)
                            
                            p_power_OR <- 0
                            p_power_IPW <- 0
                            p_power_AIPW <- 0
                            p_power_list <- list(p_power_OR,p_power_IPW,p_power_AIPW)
                            
                            p_cover_OR <- 0
                            p_cover_IPW <- 0
                            p_cover_AIPW <- 0
                            p_cover_list <- list(p_cover_OR,p_cover_IPW,p_cover_AIPW)
                            
                            return_list <- list(omega=omega_hat_list,p=p_power_list,p_cover=p_cover_list
                                                ,seed=seed)
                            
                            #print(return_list)
                            return(return_list)  #还要返回p值
                            
                          }, error = function(e) {
                            #warnings(paste("Iteration", ite_k, "failed:", e$message))
                            return(NULL)
                            # message_list <- list(message = paste("Iteration", ite_k, "failed:", e$message)) 
                            # return(message_list)
                          })
                          
                        }
                      gc()
                      #-----------------由于内存问题，部分结果可能计算不出来
                      indice_errorbyRAM <- sapply(beta_hat_ite, is.null)
                      indice_noerrorbyRAM <- !indice_errorbyRAM
                      index_noerrorbyRAM <- c(1:length(indice_noerrorbyRAM))[indice_noerrorbyRAM]
                      num_effect_ites <- length(index_noerrorbyRAM)
                      
                      beta_hat_ite <- beta_hat_ite[indice_noerrorbyRAM]
                      #---
                      param_df <- param_df[param_df[[1]] %in% index_noerrorbyRAM, ]
                      param_df2[[1]] <- num_effect_ites
                      # 结果保存
                      #+++++++++++++++ 1
                      
                      # 加载必要的包
                      library(purrr)
                      
                      # 获取总模拟次数
                      ite_times_actual <- length(beta_hat_ite)
                      
                      # 初始化三个矩阵（行数 = 3*ite_times，列数 = 2）
                      omega_matrix <- matrix(nrow = 3 * ite_times_actual, ncol = 2)
                      p_power_matrix <- matrix(nrow = 3 * ite_times_actual, ncol = 2)
                      p_cover_matrix <- matrix(nrow = 3 * ite_times_actual, ncol = 2)
                      
                      # 设置列名（按参数命名）
                      colnames(omega_matrix) <- paste0('omega',c(1:p_C))
                      colnames(p_power_matrix) <- paste0('p_power',c(1:p_C))
                      colnames(p_cover_matrix) <- paste0('p_cover',c(1:p_C))
                      
                      # 遍历每个模拟结果并填充矩阵
                      for (i in 1:ite_times_actual) {
                        # 提取第 i 次模拟的子列表
                        sim <- beta_hat_ite[[i]]
                        
                        # 计算当前模拟在矩阵中的起始行位置
                        start_row <- (i - 1) * 3 + 1
                        
                        # 填充 OR 方法的结果（第 1 行）
                        omega_matrix[start_row, ] <- sim$omega[[1]]
                        p_power_matrix[start_row, ] <- sim$p[[1]]
                        p_cover_matrix[start_row, ] <- sim$p_cover[[1]]
                        
                        # 填充 IPW 方法的结果（第 2 行）
                        omega_matrix[start_row + 1, ] <- sim$omega[[2]]
                        p_power_matrix[start_row + 1, ] <- sim$p[[2]]
                        p_cover_matrix[start_row + 1, ] <- sim$p_cover[[2]]
                        
                        # 填充 AIPW 方法的结果（第 3 行）
                        omega_matrix[start_row + 2, ] <- sim$omega[[3]]
                        p_power_matrix[start_row + 2, ] <- sim$p[[3]]
                        p_cover_matrix[start_row + 2, ] <- sim$p_cover[[3]]
                      }
                      
                      save_df_ite <- cbind(param_df,omega_matrix,p_power_matrix,p_cover_matrix)
                      save_df_itep2<- cbind(param_df,omega_matrix,p_power_matrix,p_cover_matrix)
                      
                      save_df <- rbind(save_df,save_df_ite)
                      
                      #+++++++++++++++ 2
                      omega_matrix2 <- matrix(nrow = ite_times_actual, ncol = len_my_methods * p_C)
                      colnames(omega_matrix2) <- rep(paste0('omega',c(1:p_C)),len_my_methods)
                      for (i in 1:ite_times_actual) {
                        # 计算当前模拟的起始行和结束行
                        start_row <- (i - 1) * 3 + 1
                        end_row <- i * 3
                        sim_data <- omega_matrix[start_row:end_row, ]
                        omega_matrix2[i, ] <- c(sim_data[1, ], sim_data[2, ], sim_data[3, ])
                      }
                      
                      bias_vec <- as.vector(apply(omega_matrix2-rep(omega,3), 2, function(x) mean(x[abs(x) < 10], na.rm = TRUE)))
                      sd_vec <- as.vector(apply(omega_matrix2, 2, function(x) sd(x[abs(x) < 10], na.rm = TRUE)))
                      mse_vec <- bias_vec^2 + sd_vec^2
                      
                      bias_mat <- matrix(bias_vec,nrow = len_my_methods,byrow = T)
                      sd_mat <- matrix(sd_vec,nrow = len_my_methods,byrow = T)
                      mse_mat <- matrix(mse_vec,nrow = len_my_methods,byrow = T)
                      
                      save_df_ite2 <- cbind(param_df2,rbind(bias_mat,sd_mat,mse_mat))
                      save_df2 <- rbind(save_df2,save_df_ite2)
                      #+++++++++++++++ 3
                      if(condition9){
                        # 手算power
                        test_df <- save_df_ite
                        colnames(test_df) <- col_names
                        indice_power_list <- which(n_0 == n_0_list2)
                        save_power_list[[indice_power_list]] <- test_df
                        
                        # bootstrap算power
                        test_df2 <- save_df_itep2
                        colnames(test_df2) <- col_names
                        # power_vec <- c()
                        # cover_vec <- c()
                        for (m1 in 1:ncol(sd_mat)) {
                          sd_mid <- sd_mat[,m1]
                          sd_mid2 <- rep(sd_mid,ite_times_actual)
                          
                          omega_mid <- save_df_itep2[,12+m1]
                          
                          power_mid <- 2 * (1 - pnorm(abs(omega_mid / sd_mid2)))
                          cover_mid <- 2 * (1 - pnorm(abs((omega_mid-omega[m1]) / sd_mid2)))
                          
                          # power_vec <- c(power_vec,power_mid)
                          # cover_vec <- c(cover_vec,cover_mid)
                          test_df2[,14+m1] <- power_mid
                          test_df2[,16+m1] <- cover_mid
                        }
                        
                        save_power_list2[[indice_power_list]] <- test_df2
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
}
colnames(save_df) <- col_names
colnames(save_df2) <- col_names2
#-----------------------------------------------写文件
fwrite(save_df,file=paste(path_output_tables, '/data_all_',date_mark,'_',ver,'.csv',sep = ''))
fwrite(save_df2,file=paste(path_output_tables, '/statistics_all_',date_mark,'_',ver,'.csv',sep = ''))
#===============================================================================绘图1
# time_end <- Sys.time()
# time_consume <- difftime(time_end,time_start)
# time_consume 
stopImplicitCluster()
closeAllConnections()
