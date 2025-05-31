#time_start <- Sys.time()
# 设置并行计算的核心数量
#getDoParWorkers()

cores <- detectCores()
registerDoParallel(cores = cores-4)
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
n_1_list <- c(400)
gamma_list <- c(0,0.2,0.4,0.6,0.8,1) #gamma大小的调节参数
heter_level_list <- c(0,0.2,0.4,0.6,0.8,1)
#-----
p_Z_list <- c(3,2,4,6,8,10)
#-----
inter_strength_list <- c(0,0.2,0.4,0.6,0.8,1)
alpha_list <- c(0,0.2,0.4,0.6,0.8,1)
beta_list <- c(0.05,0.1,0.2,0.3,0.4)
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
other_methods <- c("IVW","Egger","Weighted median","Weighted mode","MR-RAPS","MR-Conmix", "MR-Robust", "MR-Lasso")
my_methods <- c('TLMR-OR'
                ,'TLMR-IPW'
                ,'TLMR-AIPW')
method_names <- c(other_methods,my_methods)
my_methods2 <- c('TLMR-OR')
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
len_beta <- length(beta_list)
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
#--------------------------------名称排序
method_names_ordered <- factor(method_names,levels = method_names)
index_ordered <- factor(c('Bias','SD','MSE'),levels = c('Bias','SD','MSE'))
#--------------------------------结果保存
save_df <- data.frame()  #结果保存
col_names <- c('iteration time','beta','p_Z', 'gamma','heterogeneity','n_0','n_1','inter_ratio','inter_strength','alpha_modifier','eta_modifier','method','beta_hat','p_power','p_cover')

save_df2 <- data.frame()  #结果保存
col_names2 <- c('max iteration times','beta','p_Z', 'gamma','heterogeneity','n_0','n_1','inter_ratio','inter_strength','alpha_modifier','eta_modifier','method','index','num')

save_power_list <- vector('list',length = len_n0)
save_power_list2 <- vector('list',length = len_n0)
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
k_beta <- 5
k_beta2 <- 1
#k_ratio_pl <- len_ratio_pleiotropy
all_params_list <- vector('list') 
for(i_n0 in 1:len_n0){for (i_n1 in 1:len_n1){for(i_gamma in 1:len_gamma){for (i_heter in 1:len_heter){for(i_inter in 1:len_inter){for (i_pz in 1:len_pz) {for (i_inter_stren in 1:len_inter_strength) {for (i_alpha in 1:len_alpha) {for (i_eta in 1:len_eta) {for (i_beta in 1:len_beta) {for (i_n02 in 1:len_n02) {
  beta <- beta_list[i_beta]
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
  condition1 <- (heter_level == heter_level_list[k_heter] & n_0 == n_0_list[k_n0] & ratio_inter == ratio_inter_list[k_inter] & p_Z==p_Z_list[k_pz] & inter_strength==inter_strength_list[k_inter_strength] & alpha_modifier == alpha_list[k_alpha] & eta_modifier == eta_list[k_eta]& beta == beta_list[k_beta]& n_02 == n_0_list2[k_n02])
  # 异质性水平
  condition2 <- (gamma_multi == gamma_list[k_gamma] & n_0 == n_0_list[k_n0] & ratio_inter == ratio_inter_list[k_inter] & p_Z==p_Z_list[k_pz] & inter_strength==inter_strength_list[k_inter_strength] & alpha_modifier == alpha_list[k_alpha] & eta_modifier == eta_list[k_eta]& beta == beta_list[k_beta]& n_02 == n_0_list2[k_n02])
  # 样本量n0
  condition3 <- (heter_level == heter_level_list[k_heter] & gamma_multi == gamma_list[k_gamma]  & ratio_inter == ratio_inter_list[k_inter] & p_Z==p_Z_list[k_pz] & inter_strength==inter_strength_list[k_inter_strength] & alpha_modifier == alpha_list[k_alpha] & eta_modifier == eta_list[k_eta]& beta == beta_list[k_beta]& n_02 == n_0_list2[k_n02])
  # 修饰比率
  condition4 <- (heter_level == heter_level_list[k_heter] & gamma_multi == gamma_list[k_gamma] & n_0 == n_0_list[k_n0] & p_Z==p_Z_list[k_pz] & inter_strength==inter_strength_list[k_inter_strength] & alpha_modifier == alpha_list[k_alpha] & eta_modifier == eta_list[k_eta]& beta == beta_list[k_beta]& n_02 == n_0_list2[k_n02])
  # 工具个数：p_Z
  condition5 <- (heter_level == heter_level_list[k_heter] & gamma_multi == gamma_list[k_gamma] & n_0 == n_0_list[k_n0] & ratio_inter == ratio_inter_list[k_inter] & inter_strength==inter_strength_list[k_inter_strength] & alpha_modifier == alpha_list[k_alpha] & eta_modifier == eta_list[k_eta]& beta == beta_list[k_beta]& n_02 == n_0_list2[k_n02])
  # 修饰强度
  condition6 <- (heter_level == heter_level_list[k_heter] & gamma_multi == gamma_list[k_gamma2] & n_0 == n_0_list[k_n0] & ratio_inter == ratio_inter_list[k_inter] & p_Z==p_Z_list[k_pz] & alpha_modifier == alpha_list[k_alpha2] & eta_modifier == eta_list[k_eta]& beta == beta_list[k_beta]& n_02 == n_0_list2[k_n02])
  # eta
  condition7 <- (heter_level == heter_level_list[k_heter] & gamma_multi == gamma_list[k_gamma] & n_0 == n_0_list[k_n0] & ratio_inter == ratio_inter_list[k_inter] & p_Z==p_Z_list[k_pz] & inter_strength==inter_strength_list[k_inter_strength] & alpha_modifier == alpha_list[k_alpha] & beta == beta_list[k_beta]& n_02 == n_0_list2[k_n02])
  # alpha
  condition8 <- (heter_level == heter_level_list[k_heter] & gamma_multi == gamma_list[k_gamma] & n_0 == n_0_list[k_n0] & ratio_inter == ratio_inter_list[k_inter] & p_Z==p_Z_list[k_pz] & inter_strength==inter_strength_list[k_inter_strength] & eta_modifier == eta_list[k_eta]& beta == beta_list[k_beta] & n_02 == n_0_list2[k_n02] )
  # 另一个beta下遍历样本量
  condition9 <- (heter_level == heter_level_list[k_heter] & gamma_multi == gamma_list[k_gamma] & n_0 == n_0_list[k_n0] & ratio_inter == ratio_inter_list[k_inter] & p_Z==p_Z_list[k_pz] & inter_strength==inter_strength_list[k_inter_strength] & eta_modifier == eta_list[k_eta] & alpha_modifier == alpha_list[k_alpha] & beta == beta_list[k_beta2])
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
                      beta <- beta_list[i_beta]
                      
                      n_0 <- n_0_list[i_n0]
                      n_02 <- n_0_list2[i_n02]
                      n_1 <- n_0
                      
                      n_divides <- floor(log2(n_0/400))
                      
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
                                      ,heter_level=heter_level
                                      #,U_md = U_md
                      )
                      #heritability <- 0.7
                      pval_i <- 0.05
                      #0.05/p_Z
                      #-------------------
                      condition1 <- (heter_level == heter_level_list[k_heter] & n_0 == n_0_list[k_n0] & ratio_inter == ratio_inter_list[k_inter] & p_Z==p_Z_list[k_pz] & inter_strength==inter_strength_list[k_inter_strength] & alpha_modifier == alpha_list[k_alpha] & eta_modifier == eta_list[k_eta] & beta == beta_list[k_beta]& n_02 == n_0_list2[k_n02])
                      condition2 <- (gamma_multi == gamma_list[k_gamma] & n_0 == n_0_list[k_n0] & ratio_inter == ratio_inter_list[k_inter] & p_Z==p_Z_list[k_pz] & inter_strength==inter_strength_list[k_inter_strength] & alpha_modifier == alpha_list[k_alpha] & eta_modifier == eta_list[k_eta] & beta == beta_list[k_beta]& n_02 == n_0_list2[k_n02])
                      condition3 <- (heter_level == heter_level_list[k_heter] & gamma_multi == gamma_list[k_gamma]  & ratio_inter == ratio_inter_list[k_inter] & p_Z==p_Z_list[k_pz] & inter_strength==inter_strength_list[k_inter_strength] & alpha_modifier == alpha_list[k_alpha] & eta_modifier == eta_list[k_eta] & beta == beta_list[k_beta]& n_02 == n_0_list2[k_n02])
                      condition4 <- (heter_level == heter_level_list[k_heter] & gamma_multi == gamma_list[k_gamma] & n_0 == n_0_list[k_n0] & p_Z==p_Z_list[k_pz] & inter_strength==inter_strength_list[k_inter_strength] & alpha_modifier == alpha_list[k_alpha] & eta_modifier == eta_list[k_eta] & beta == beta_list[k_beta]& n_02 == n_0_list2[k_n02])
                      condition5 <- (heter_level == heter_level_list[k_heter] & gamma_multi == gamma_list[k_gamma] & n_0 == n_0_list[k_n0] & ratio_inter == ratio_inter_list[k_inter] & inter_strength==inter_strength_list[k_inter_strength] & alpha_modifier == alpha_list[k_alpha] & eta_modifier == eta_list[k_eta]& beta == beta_list[k_beta]& n_02 == n_0_list2[k_n02])
                      condition6 <- (heter_level == heter_level_list[k_heter] & gamma_multi == gamma_list[k_gamma2] & n_0 == n_0_list[k_n0] & ratio_inter == ratio_inter_list[k_inter] & p_Z==p_Z_list[k_pz] & alpha_modifier == alpha_list[k_alpha2] & eta_modifier == eta_list[k_eta]& beta == beta_list[k_beta]& n_02 == n_0_list2[k_n02])
                      condition7 <- (heter_level == heter_level_list[k_heter] & gamma_multi == gamma_list[k_gamma] & n_0 == n_0_list[k_n0] & ratio_inter == ratio_inter_list[k_inter] & p_Z==p_Z_list[k_pz] & inter_strength==inter_strength_list[k_inter_strength] & alpha_modifier == alpha_list[k_alpha] & beta == beta_list[k_beta]& n_02 == n_0_list2[k_n02])
                      condition8 <- (heter_level == heter_level_list[k_heter] & gamma_multi == gamma_list[k_gamma] & n_0 == n_0_list[k_n0] & ratio_inter == ratio_inter_list[k_inter] & p_Z==p_Z_list[k_pz] & inter_strength==inter_strength_list[k_inter_strength] & eta_modifier == eta_list[k_eta] & beta == beta_list[k_beta]& n_02 == n_0_list2[k_n02])
                      condition9 <- (heter_level == heter_level_list[k_heter] & gamma_multi == gamma_list[k_gamma] & n_0 == n_0_list[k_n0] & ratio_inter == ratio_inter_list[k_inter] & p_Z==p_Z_list[k_pz] & inter_strength==inter_strength_list[k_inter_strength] & eta_modifier == eta_list[k_eta] & alpha_modifier == alpha_list[k_alpha] & beta == beta_list[k_beta2])
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
                      param_df2 <- data.frame(matrix(rep(c(ite_times,beta,p_Z,gamma_multi,heter_level,n_0,n_1,ratio_inter,inter_strength,alpha_modifier,eta_modifier),len_methods*3),nrow = len_methods*3,byrow = TRUE))
                      param_df2 <- cbind(param_df2,rep(method_names_ordered,3))
                      param_df2 <- cbind(param_df2,rep(index_ordered,each=len_methods))
                      
                      
                      param_df <- data.frame(matrix(rep(c(beta,p_Z,gamma_multi,heter_level,n_0,n_1,ratio_inter,inter_strength,alpha_modifier,eta_modifier),len_methods*ite_times),nrow = len_methods*ite_times,byrow = TRUE))
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
                            #----------------
                            quantile_vec_all <- seq(0,1,1/n_divides)
                            quantile_vec <- quantile_vec_all[-c(1,length(quantile_vec_all))]
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
                            conditions_met <- FALSE
                            while(!conditions_met){
                              #------------------------------
                              C_all <- rbind(C_0,C_1)
                              
                              C_0_cl <- C_0
                              C_1_cl <- C_1
                              for (i in c(1:dim(C_0)[2])) {
                                if( length(unique(C_0[,i])) > 10 ){
                                  #quant_vec <- c(0.3333,0.6667)
                                  C_qt_vector <- quantile(C_all[,i],quantile_vec)
                                  C_qt_label <- cut( C_0[,i],c(-Inf,C_qt_vector,Inf),labels=c(1:n_divides) )
                                  C_0_cl[,i] <- C_qt_label
                                  C_qt_label <- cut( C_1[,i],c(-Inf,C_qt_vector,Inf),labels=c(1:n_divides) )
                                  C_1_cl[,i] <- C_qt_label
                                }
                              }
                              # 替换成分类变量；原始数据存储至origin
                              C_0_origin <- C_0
                              C_1_origin <- C_1
                              C_0 <- C_0_cl
                              C_1 <- C_1_cl
                              #------------------------------
                              
                              condition_list <- rep(FALSE,p_Z)
                              for (i in 1:p_Z ){
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
                                #seed <- seed + 10 
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
                              }
                            }
                            #----------------------------
                            #Twosamle_package2(p_Z,p_C,data_0,data_1)
                            #===============================================================Twosample-MR包     1
                            beta_hat_1_list <- Twosamle_package(p_Z=p_Z,p_C=p_C,data_0=data_0,data_1=data_1,beta=beta)
                            beta_hat_1 <- beta_hat_1_list$beta
                            #===============================================================Trans-MR
                            #----------------------------看是否存在修饰作用的代码
                            modifier_matrix <- matrix(nrow = p_Z,ncol = p_C)
                            for(k in 1:p_Z){
                              for (j in 1:p_C) {
                                data_chow_test <- data.frame(Z=Z_1[,k],X=X_1,C=C_1[,j])
                                data_chow_test <- data_chow_test[order(data_chow_test$C), ]
                                C_group_vec <- table(data_chow_test$C)
                                
                                end_indice <- cumsum(C_group_vec)
                                start_indice <- c(1,head(end_indice,-1)+1)
                                
                                beta_group_vec <- c()
                                se_group_vec <- c()
                                for(i in 1:length(C_group_vec)){
                                  data_chow_test_sub <- data_chow_test[(start_indice[i]:end_indice[i]),]
                                  model_group_C <- lm(X~Z,data = data_chow_test_sub)
                                  res_model_group <- summary(model_group_C)$coefficients
                                  beta_group_vec <- c(beta_group_vec,res_model_group[2,1])
                                  se_group_vec <- c(se_group_vec,res_model_group[2,2])
                                }
                                res_meta_analysis <- metafor::rma.uni(yi=beta_group_vec,sei = se_group_vec,method = 'DL')
                                Q_statistics_group <- res_meta_analysis$QE
                                p_Q_statistics_group <- res_meta_analysis$QEp
                                if(p_Q_statistics_group <(0.005)){
                                  modifier_matrix[k,j] <- 1
                                }else{
                                  modifier_matrix[k,j] <- 0
                                }
                              }
                            }
                            #----------------------------
                            b_Y_list2 <- c()
                            sd_Y_list2 <- c()
                            b_Y_list3 <- c()
                            sd_Y_list3 <- c()
                            b_Y_list4 <- c()
                            sd_Y_list4 <- c()
                            
                            indice_modify <- which(rowSums(modifier_matrix) > 0)
                            indice_notmodify <- setdiff(c(1:p_Z),which(rowSums(modifier_matrix) > 0))
                            p_Zm <- length(indice_modify)
                            p_Znm <- length(indice_notmodify)
                            
                            
                            if(p_Znm>0){
                              
                              data3_0 <- data_0
                              data3_1 <- data_1
                              
                              data3_0$Z <- data3_0$Z[,indice_notmodify,drop=FALSE]
                              data3_0$beta <- data3_0$beta[indice_notmodify,drop=FALSE]
                              data3_0$sd <- data3_0$sd[indice_notmodify,drop=FALSE]
                              data3_0$eaf <- data3_0$eaf[indice_notmodify,drop=FALSE]
                              data3_0$R2 <- data3_0$R2[indice_notmodify,drop=FALSE]
                              data3_0$F_stats <- data3_0$F_stats[indice_notmodify,drop=FALSE]
                              
                              data3_1$Z <- data3_1$Z[,indice_notmodify,drop=FALSE]
                              data3_1$beta <- data3_1$beta[indice_notmodify,drop=FALSE]
                              data3_1$sd <- data3_1$sd[indice_notmodify,drop=FALSE]
                              data3_1$eaf <- data3_1$eaf[indice_notmodify,drop=FALSE]
                              data3_1$R2 <- data3_1$R2[indice_notmodify,drop=FALSE]
                              data3_1$F_stats <- data3_1$F_stats[indice_notmodify,drop=FALSE]
                              
                              beta_hat_res <- Twosamle_package3(p_Z=length(indice_notmodify),p_C=p_C,data_0=data3_0,data_1=data3_1)
                              beta_hat_degrade <- beta_hat_res$beta
                              se_degrade <- beta_hat_res$se
                              #----------------------------------
                              if(p_Zm == 0){
                                beta_hat_2 <- beta_hat_degrade
                                beta_hat_3 <- beta_hat_degrade
                                beta_hat_4 <- beta_hat_degrade
                                b_Y_list2 <- c(b_Y_list2,beta_hat_degrade)
                                b_Y_list3 <- c(b_Y_list3,beta_hat_degrade)
                                b_Y_list4 <- c(b_Y_list4,beta_hat_degrade)
                                sd_Y_list2 <- c(sd_Y_list2,se_degrade)
                                sd_Y_list3 <- c(sd_Y_list3,se_degrade)
                                sd_Y_list4 <- c(sd_Y_list4,se_degrade)
                                sd_Y_list2_inv <- 1/sd_Y_list2
                                sd_Y_list3_inv <- 1/sd_Y_list3
                                sd_Y_list4_inv <- 1/sd_Y_list4
                              }else{
                                Z_0m <- Z_0[,indice_modify,drop=FALSE]
                                Z_1m <- Z_1[,indice_modify,drop=FALSE]
                                # if(sum(indice_modify) == 1){
                                #   Z_0m <- matrix(Z_0m,ncol = 1)
                                #   Z_1m <- matrix(Z_1m,ncol = 1)
                                # }
                                #===============================================================迁移估计量
                                #Trans_OR
                                X_hat_matrix_OR <- Trans_OLS2(Z_0=Z_0m,Z_1=Z_1m,C_0=C_0_origin,C_1=C_1_origin,X_1=X_1)$X_hat_matrix
                                #Trans_IPW
                                X_hat_matrix_IPW <- IPW2(Z_0=Z_0m,Z_1=Z_1m,C_0=C_0,C_1=C_1,X_1=X_1)$X_hat_matrix
                                #Trans_AIPW
                                X_hat_matrix_AIPW <- AIPW2(Z_0=Z_0m,Z_1=Z_1m,C_0=C_0,C_1=C_1,X_1=X_1
                                                           ,C_0_origin=C_0_origin,C_1_origin=C_1_origin)$X_hat_matrix
                                #-----------------------------------------
                                
                                #---------------------------------------------------------------Trans_OR0         21
                                if(p_Zm >= 1){
                                  for (i in 1:p_Zm) {  
                                    model2 <- lm(Y_0~X_hat_matrix_OR[,i] + C_0_origin + Z_0m[,i])  #有Z
                                    summary_model2 <- summary(model2)
                                    pval_model2 <- summary_model2$coefficients[,4]
                                    Z_pval <- pval_model2[length(pval_model2)]
                                    if(Z_pval < pval_i){
                                      b_Y_list2[i] <- model2$coefficients[2]
                                      sd_Y_list2[i] <- summary(model2)$coefficients[,'Std. Error'][2]
                                    }else{
                                      model2_1 <- lm(Y_0~X_hat_matrix_OR[,i] + C_0_origin)  #无Z
                                      b_Y_list2[i] <- model2_1$coefficients[2]
                                      sd_Y_list2[i] <- summary(model2_1)$coefficients[,'Std. Error'][2]
                                    }
                                  }
                                }
                                b_Y_list2 <- c(b_Y_list2, beta_hat_degrade)
                                sd_Y_list2 <- c(sd_Y_list2, se_degrade)
                                
                                sd_Y_list2_inv <- 1/sd_Y_list2 
                                beta_hat_2 <- sum(sd_Y_list2_inv^2 * b_Y_list2)/sum(sd_Y_list2_inv^2)
                                #---------------------------------------------------------------Trans_IPW0          31
                                if(p_Zm>=1){
                                  
                                  for (i in 1:p_Zm) {
                                    
                                    model3 <- lm(Y_0~X_hat_matrix_IPW[,i] + C_0_origin + Z_0m[,i])  #有Z
                                    summary_model3 <- summary(model3)
                                    pval_model3 <- summary_model3$coefficients[,4]
                                    Z_pval <- pval_model3[length(pval_model3)]
                                    if(Z_pval < pval_i){
                                      b_Y_list3[i] <- model3$coefficients[2]
                                      sd_Y_list3[i] <- summary(model3)$coefficients[,'Std. Error'][2]
                                    }else{
                                      model3_1 <- lm(Y_0~X_hat_matrix_IPW[,i] + C_0_origin + Z_0[,i])  #无Z
                                      b_Y_list3[i] <- model3_1$coefficients[2]
                                      sd_Y_list3[i] <- summary(model3_1)$coefficients[,'Std. Error'][2]
                                    }
                                  }
                                }
                                
                                b_Y_list3 <- c(b_Y_list3, beta_hat_degrade)
                                sd_Y_list3 <- c(sd_Y_list3, se_degrade)
                                
                                sd_Y_list3_inv <- 1/sd_Y_list3 
                                beta_hat_3 <- sum(sd_Y_list3_inv^2 * b_Y_list3)/sum(sd_Y_list3_inv^2)
                                #---------------------------------------------------------------Trans_AIPW0          41
                                if(p_Zm>=1){
                                  
                                  for (i in 1:p_Zm) {
                                    model4 <- lm(Y_0~X_hat_matrix_AIPW[,i] + C_0_origin + Z_0m[,i])  #有Z
                                    summary_model4 <- summary(model4)
                                    pval_model4 <- summary_model4$coefficients[,4]
                                    Z_pval <- pval_model4[length(pval_model4)]
                                    if(Z_pval<pval_i){
                                      b_Y_list4[i] <- model4$coefficients[2]
                                      sd_Y_list4[i] <- summary(model4)$coefficients[,'Std. Error'][2]
                                    }else{
                                      model4_1 <- lm(Y_0~X_hat_matrix_AIPW[,i] + C_0_origin + Z_0[,i])  #无Z
                                      b_Y_list4[i] <- model4_1$coefficients[2]
                                      sd_Y_list4[i] <- summary(model4_1)$coefficients[,'Std. Error'][2]
                                    }
                                  }
                                }
                                b_Y_list4 <- c(b_Y_list4, beta_hat_degrade)
                                sd_Y_list4 <- c(sd_Y_list4, se_degrade)
                                
                                sd_Y_list4_inv <- 1/sd_Y_list4
                                beta_hat_4 <- sum(sd_Y_list4_inv^2 * b_Y_list4)/sum(sd_Y_list4_inv^2)
                              }
                            }else{
                              #===============================================================迁移估计量
                              #Trans_OR
                              X_hat_matrix_OR <- Trans_OLS2(Z_0=Z_0,Z_1=Z_1,C_0=C_0_origin,C_1=C_1_origin,X_1=X_1)$X_hat_matrix
                              #Trans_IPW
                              X_hat_matrix_IPW <- IPW2(Z_0=Z_0,Z_1=Z_1,C_0=C_0,C_1=C_1,X_1=X_1)$X_hat_matrix
                              #Trans_AIPW
                              X_hat_matrix_AIPW <- AIPW2(Z_0=Z_0,Z_1=Z_1,C_0=C_0,C_1=C_1,X_1=X_1
                                                         ,C_0_origin=C_0_origin,C_1_origin=C_1_origin)$X_hat_matrix
                              #---------------------------------------------------------------Trans_OR0         21
                              for (i in 1:p_Z) {  
                                model2 <- lm(Y_0~X_hat_matrix_OR[,i] + C_0_origin + Z_0[,i])  #有Z
                                summary_model2 <- summary(model2)
                                pval_model2 <- summary_model2$coefficients[,4]
                                Z_pval <- pval_model2[length(pval_model2)]
                                #print(as.numeric(Z_pval))
                                #if(Z_pval < pval_i){
                                if(Z_pval < pval_i){
                                  b_Y_list2[i] <- model2$coefficients[2]
                                  sd_Y_list2[i] <- summary(model2)$coefficients[,'Std. Error'][2]
                                }else{
                                  model2_1 <- lm(Y_0~X_hat_matrix_OR[,i] + C_0_origin)  #无Z
                                  b_Y_list2[i] <- model2_1$coefficients[2]
                                  sd_Y_list2[i] <- summary(model2_1)$coefficients[,'Std. Error'][2]
                                }
                                
                              }
                              # cor(X_hat_matrix_OR[,1],Z_0[,1])
                              sd_Y_list2_inv <- 1/sd_Y_list2 
                              beta_hat_2 <- sum(sd_Y_list2_inv^2 * b_Y_list2)/sum(sd_Y_list2_inv^2)
                              #---------------------------------------------------------------Trans_IPW0          31
                              for (i in 1:p_Z) {
                                
                                model3 <- lm(Y_0~X_hat_matrix_IPW[,i] + C_0_origin + Z_0[,i])  #有Z
                                summary_model3 <- summary(model3)
                                pval_model3 <- summary_model3$coefficients[,4]
                                Z_pval <- pval_model3[length(pval_model3)]
                                #print(as.numeric(Z_pval))
                                if(Z_pval < pval_i){
                                  b_Y_list3[i] <- model3$coefficients[2]
                                  sd_Y_list3[i] <- summary(model3)$coefficients[,'Std. Error'][2]
                                }else{
                                  model3_1 <- lm(Y_0~X_hat_matrix_IPW[,i] + C_0_origin)  #无Z
                                  b_Y_list3[i] <- model3_1$coefficients[2]
                                  sd_Y_list3[i] <- summary(model3_1)$coefficients[,'Std. Error'][2]
                                }
                                
                                
                              }
                              sd_Y_list3_inv <- 1/sd_Y_list3 
                              beta_hat_3 <- sum(sd_Y_list3_inv^2 * b_Y_list3)/sum(sd_Y_list3_inv^2)
                              #---------------------------------------------------------------Trans_AIPW0          41
                              for (i in 1:p_Z) {
                                model4 <- lm(Y_0~X_hat_matrix_AIPW[,i] + C_0_origin + Z_0[,i])  #有Z
                                summary_model4 <- summary(model4)
                                pval_model4 <- summary_model4$coefficients[,4]
                                Z_pval <- pval_model4[length(pval_model4)]
                                if(Z_pval < pval_i){
                                  b_Y_list4[i] <- model4$coefficients[2]
                                  sd_Y_list4[i] <- summary(model4)$coefficients[,'Std. Error'][2]
                                }else{
                                  model4_1 <- lm(Y_0~X_hat_matrix_AIPW[,i] + C_0_origin)  #无Z
                                  b_Y_list4[i] <- model4_1$coefficients[2]
                                  sd_Y_list4[i] <- summary(model4_1)$coefficients[,'Std. Error'][2]
                                }
                              }
                              sd_Y_list4_inv <- 1/sd_Y_list4
                              beta_hat_4 <- sum(sd_Y_list4_inv^2 * b_Y_list4)/sum(sd_Y_list4_inv^2)
                            }
                            
                            #----------------------------结果保存
                            beta_hat_vec <- c(beta_hat_1, beta_hat_2, beta_hat_3,beta_hat_4)
                            #--------------------------------------------------IVW方差计算
                            sd_OR <- sqrt(1/sum(sd_Y_list2_inv^2))
                            sd_IPW <- sqrt(1/sum(sd_Y_list3_inv^2))
                            sd_AIPW <- sqrt(1/sum(sd_Y_list4_inv^2))
                            p_power_OR <- 2 * (1 - pnorm(abs(beta_hat_2 / sd_OR)))
                            p_power_IPW <- 2 * (1 - pnorm(abs(beta_hat_3 / sd_IPW)))
                            p_power_AIPW <- 2 * (1 - pnorm(abs(beta_hat_4 / sd_AIPW)))
                            p_power_mm <- c(p_power_OR,p_power_IPW,p_power_AIPW)
                            
                            p_cover_OR <- 2 * (1 - pnorm(abs((beta_hat_2-beta) / sd_OR)))
                            p_cover_IPW <- 2 * (1 - pnorm(abs((beta_hat_3-beta) / sd_IPW)))
                            p_cover_AIPW <- 2 * (1 - pnorm(abs((beta_hat_4-beta) / sd_AIPW)))
                            p_cover_mm <- c(p_cover_OR,p_cover_IPW,p_cover_AIPW)
                            #--------------------------------------------------
                            p_vec <- c(beta_hat_1_list$p,p_power_mm)
                            p_cover_vec <- c(beta_hat_1_list$p_cover,p_cover_mm) 
                            
                            p_vec2 <- c(beta_hat_1_list$p,rep(NA,len_methods-8))
                            p_cover_vec2 <- c(beta_hat_1_list$p_cover,rep(NA,len_methods-8)) 
                            
                            return_list <- list(beta=beta_hat_vec,p=p_vec,p_cover=p_cover_vec
                                                ,p2=p_vec2,p_cover2=p_cover_vec2,seed=seed)
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
                      ite_return <- do.call(rbind,beta_hat_ite)
                      ite_df2 <- do.call(rbind,ite_return[,1])  #beta值矩阵
                      #ite_p <- do.call(rbind,ite_return[,2])   #p值矩阵
                      #ite_p_cover <- do.call(rbind,ite_return[,3])  #p值矩阵
                      
                      beta_hat_ite_vec <- as.vector(unlist(ite_return[,1]))   # (len_methods个beta_hat) * (循环次数)  个beta_hat值
                      p_ite_vec <- as.vector(unlist(ite_return[,2]))    # (len_methods个beta_hat) * (循环次数)  个p值
                      p_cover_ite_vec <- as.vector(unlist(ite_return[,3]))   # (len_methods个beta_hat) * (循环次数)  个p值
                      
                      p_ite_vec2 <- as.vector(unlist(ite_return[,4]))    # (len_methods个beta_hat) * (循环次数)  个p值
                      p_cover_ite_vec2 <- as.vector(unlist(ite_return[,5]))   # (len_methods个beta_hat) * (循环次数)  个p值
                      
                      
                      save_df_ite <- cbind(param_df,beta_hat_ite_vec,p_ite_vec,p_cover_ite_vec)
                      save_df_itep2<- cbind(param_df,beta_hat_ite_vec,p_ite_vec2,p_cover_ite_vec2)
                      
                      save_df <- rbind(save_df,save_df_ite)
                      
                      #+++++++++++++++ 2
                      
                      bias_vec <- as.vector(apply(ite_df2-beta, 2, function(x) mean(x[abs(x) < 10], na.rm = TRUE)))
                      sd_vec <- as.vector(apply(ite_df2, 2, function(x) sd(x[abs(x) < 10], na.rm = TRUE)))
                      mse_vec <- bias_vec^2 + sd_vec^2
                      save_df_ite2 <- cbind(param_df2,c(bias_vec,sd_vec,mse_vec))
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
                        sd_vec_my_methods <- sd_vec[(len_other_methods+1):len_methods]
                        test_df2 <- test_df2 %>%
                          mutate(
                            # 如果method在my_methods中，计算p值
                            p_power = ifelse(method %in% my_methods,
                                             2 * (1 - pnorm(abs(beta_hat / sd_vec_my_methods[match(method, my_methods)]))),
                                             p_power),
                            # 如果method在my_methods中，计算p_power
                            p_cover = ifelse(method %in% my_methods,
                                             2 * (1 - pnorm(abs((beta_hat - beta) / sd_vec_my_methods[match(method, my_methods)]))),
                                             p_cover)
                          )
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