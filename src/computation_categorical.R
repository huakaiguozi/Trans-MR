#time_start <- Sys.time()
# 设置并行计算的核心数量
cores <- detectCores()
registerDoParallel(cores = cores-1)
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
n_0_list <- c(2000,4000,6000,8000,10000,12000)  #样本量
#n_0_list <- c(1000)  #样本量
n_1_list <- c(400)
gamma_list <- c(0,0.2,0.4,0.6,0.8,1) #gamma大小的调节参数
heter_level_list <- c(0,0.2,0.4,0.6,0.8,1)
#-----
p_Z_list <- c(3,4,6,8,10,12,15)
#-----
inter_strength_list <- c(0,0.2,0.4,0.6,0.8,1)
alpha_list <- c(0,0.2,0.4,0.6,0.8,1)

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

ite_times <- 5   #预计10.00结束（21.55开始）
date_mark <- format(Sys.time(), "%Y-%m-%d")
#----Twosamplemr包方法

# 共11种方法
other_methods <- c("IVW","Egger","Weighted median","Weighted mode","MR-RAPS","Contamination Mixture", "MR-Robust", "MR-Lasso")
my_methods <- c('Trans-OR'
                ,'Trans-IPW'
                ,'Trans-AIPW')
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
len_n1 <- length(n_1_list)
len_pz <- length(p_Z_list)
len_inter_strength <- length(inter_strength_list)
len_alpha <- length(alpha_list)
len_eta <- length(eta_list)

len_ratio_pleiotropy <- length(ratio_pleiotropy)
#参数空间
# beta=beta,eta_2=eta_2,gamma=gamma,ratio_pleiotropy=ratio_pleiotropy,heter_level=heter_level
pi_data1 <- matrix(round(with_seed(1, runif(90000,5,15)),1),nrow = 300)
pi_data2 <- matrix(round(with_seed(1, runif(90000,6,16)),1),nrow = 300)
alpha_data <- round(with_seed(1,runif(1000,2,4)),1)
gamma_data <- round(with_seed(1,runif(1000,0.1,0.3)),2)

#--------------------------------名称排序
method_names_ordered <- factor(method_names,levels = method_names)
index_ordered <- factor(c('Bias','SD','MSE'),levels = c('Bias','SD','MSE'))
#--------------------------------结果保存
save_df <- data.frame()  #结果保存
col_names <- c('iteration time','p_Z', 'gamma','heterogeneity','n_0','n_1','inter_ratio','inter_strength','alpha_modifier','eta_modifier','method','beta_hat','p_power','p_cover')

save_df2 <- data.frame()  #结果保存
col_names2 <- c('max iteration times','p_Z', 'gamma','heterogeneity','n_0','n_1','inter_ratio','inter_strength','alpha_modifier','eta_modifier','method','index','num')

save_power_list <- vector('list',length = len_n0)
#----------
zongxunhuancishu<-0
k_min <- 1

k_pz <- 1
p_Z_fix <- p_Z_list[k_pz]
k_gamma <- 1
k_alpha <- 1
k_alpha2 <- 2
k_heter <- len_heter
k_n0 <- 6
k_inter <- len_inter
k_inter_strength <- len_inter_strength
k_eta <- len_eta

#k_ratio_pl <- len_ratio_pleiotropy
all_params_list <- vector('list') 
for(i_n0 in 1:len_n0){for (i_n1 in 1:len_n1){for(i_gamma in 1:len_gamma){for (i_heter in 1:len_heter){for(i_inter in 1:len_inter){for (i_pz in 1:len_pz) {for (i_inter_stren in 1:len_inter_strength) {for (i_alpha in 1:len_alpha) {for (i_eta in 1:len_eta) {
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
  condition1 <- (heter_level == heter_level_list[k_heter] & n_0 == n_0_list[k_n0] & ratio_inter == ratio_inter_list[k_inter] & p_Z==p_Z_list[k_pz] & inter_strength==inter_strength_list[k_inter_strength] & alpha_modifier == alpha_list[k_alpha] & eta_modifier == eta_list[k_eta])
  # 异质性水平
  condition2 <- (gamma_multi == gamma_list[k_gamma] & n_0 == n_0_list[k_n0] & ratio_inter == ratio_inter_list[k_inter] & p_Z==p_Z_list[k_pz] & inter_strength==inter_strength_list[k_inter_strength] & alpha_modifier == alpha_list[k_alpha] & eta_modifier == eta_list[k_eta])
  # 样本量n0
  condition3 <- (heter_level == heter_level_list[k_heter] & gamma_multi == gamma_list[k_gamma]  & ratio_inter == ratio_inter_list[k_inter] & p_Z==p_Z_list[k_pz] & inter_strength==inter_strength_list[k_inter_strength] & alpha_modifier == alpha_list[k_alpha] & eta_modifier == eta_list[k_eta])
  # 修饰比率
  condition4 <- (heter_level == heter_level_list[k_heter] & gamma_multi == gamma_list[k_gamma] & n_0 == n_0_list[k_n0] & p_Z==p_Z_list[k_pz] & inter_strength==inter_strength_list[k_inter_strength] & alpha_modifier == alpha_list[k_alpha] & eta_modifier == eta_list[k_eta])
  # 工具个数：p_Z
  condition5 <- (heter_level == heter_level_list[k_heter] & gamma_multi == gamma_list[k_gamma] & n_0 == n_0_list[k_n0] & ratio_inter == ratio_inter_list[k_inter] & inter_strength==inter_strength_list[k_inter_strength] & alpha_modifier == alpha_list[k_alpha] & eta_modifier == eta_list[k_eta])
  # 修饰强度
  condition6 <- (heter_level == heter_level_list[k_heter] & gamma_multi == gamma_list[k_gamma] & n_0 == n_0_list[k_n0] & ratio_inter == ratio_inter_list[k_inter] & p_Z==p_Z_list[k_pz] & alpha_modifier == alpha_list[k_alpha2] & eta_modifier == eta_list[k_eta])
  # eta
  condition7 <- (heter_level == heter_level_list[k_heter] & gamma_multi == gamma_list[k_gamma] & n_0 == n_0_list[k_n0] & ratio_inter == ratio_inter_list[k_inter] & p_Z==p_Z_list[k_pz] & inter_strength==inter_strength_list[k_inter_strength] & alpha_modifier == alpha_list[k_alpha] )
  # alpha
  condition8 <- (heter_level == heter_level_list[k_heter] & gamma_multi == gamma_list[k_gamma] & n_0 == n_0_list[k_n0] & ratio_inter == ratio_inter_list[k_inter] & p_Z==p_Z_list[k_pz] & inter_strength==inter_strength_list[k_inter_strength] & eta_modifier == eta_list[k_eta])
  
  if(condition1|condition2|condition3|condition4|condition5|condition6|condition7|condition8){
    zongxunhuancishu<-zongxunhuancishu+1
    all_params_list[[zongxunhuancishu]] <- param_list
  }
  
}}}}}}}}}

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
                  n_0 <- n_0_list[i_n0]
                  n_1 <- n_0
                  
                  heter_level <- heter_level_list[i_heter]
                  ratio_inter <- ratio_inter_list[i_inter]
                  gamma_multi <- gamma_list[i_gamma]
                  
                  inter_strength <- inter_strength_list[i_inter_stren]
                  alpha_modifier <- alpha_list[i_alpha]
                  eta_modifier <- eta_list[i_eta]
                  #--------------------------
                  p_Z <- p_Z_list[i_pz]
                  if(p_Z == p_Z_list[k_pz]){
                    p_Z_max = p_Z
                  }else{
                    p_Z_max = p_Z_list[len_pz]
                  }
                  
                  p_C <- 2
                  
                  beta <- 1
                  eta_1 <- rep(0.15*eta_modifier, p_C)
                  eta_2 <- rep(0.18*eta_modifier, p_C)
                  
                  #-----
                  #alpha <- seq(2,4,2/(p_Z-1))*alpha_modifier
                  alpha <- alpha_data[1:p_Z_max] * alpha_modifier
                  pi <- pi_data1[1:p_Z_max,1:p_C] * inter_strength
                  pi2 <- pi_data2[1:p_Z_max,1:p_C] * inter_strength
                  
                  #gamma <- seq(0.3,0.1,-0.2/(p_Z-1)) * gamma_multi
                  gamma <- gamma_data[1:p_Z_max] * gamma_multi
                  #rep(0.2,p_Z) * gamma_multi
                  
                  alpha_true <- alpha_data[1:p_Z] * alpha_modifier
                  pi_true <- pi_data1[1:p_Z,1:p_C] * inter_strength
                  pi2_true <- pi_data2[1:p_Z,1:p_C] * inter_strength
                  gamma_true <- gamma_data[1:p_Z] * gamma_multi
                  #----
                  param_X <- list(alpha=alpha,pi=pi,pi2=pi2,eta_1=eta_1)
                  param_Y <- list(beta=beta,eta_2=eta_2,gamma=gamma,ratio_pleiotropy=ratio_pleiotropy,heter_level=heter_level)
                  #heritability <- 0.7
                  #-------------------
                  condition1 <- (heter_level == heter_level_list[k_heter] & n_0 == n_0_list[k_n0] & ratio_inter == ratio_inter_list[k_inter] & p_Z==p_Z_list[k_pz] & inter_strength==inter_strength_list[k_inter_strength] & alpha_modifier == alpha_list[k_alpha] & eta_modifier == eta_list[k_eta])
                  condition2 <- (gamma_multi == gamma_list[k_gamma] & n_0 == n_0_list[k_n0] & ratio_inter == ratio_inter_list[k_inter] & p_Z==p_Z_list[k_pz] & inter_strength==inter_strength_list[k_inter_strength] & alpha_modifier == alpha_list[k_alpha] & eta_modifier == eta_list[k_eta])
                  condition3 <- (heter_level == heter_level_list[k_heter] & gamma_multi == gamma_list[k_gamma]  & ratio_inter == ratio_inter_list[k_inter] & p_Z==p_Z_list[k_pz] & inter_strength==inter_strength_list[k_inter_strength] & alpha_modifier == alpha_list[k_alpha] & eta_modifier == eta_list[k_eta])
                  condition4 <- (heter_level == heter_level_list[k_heter] & gamma_multi == gamma_list[k_gamma] & n_0 == n_0_list[k_n0] & p_Z==p_Z_list[k_pz] & inter_strength==inter_strength_list[k_inter_strength] & alpha_modifier == alpha_list[k_alpha] & eta_modifier == eta_list[k_eta])
                  condition5 <- (heter_level == heter_level_list[k_heter] & gamma_multi == gamma_list[k_gamma] & n_0 == n_0_list[k_n0] & ratio_inter == ratio_inter_list[k_inter] & inter_strength==inter_strength_list[k_inter_strength] & alpha_modifier == alpha_list[k_alpha] & eta_modifier == eta_list[k_eta])
                  condition6 <- (heter_level == heter_level_list[k_heter] & gamma_multi == gamma_list[k_gamma] & n_0 == n_0_list[k_n0] & ratio_inter == ratio_inter_list[k_inter] & p_Z==p_Z_list[k_pz] & alpha_modifier == alpha_list[k_alpha2] & eta_modifier == eta_list[k_eta])
                  condition7 <- (heter_level == heter_level_list[k_heter] & gamma_multi == gamma_list[k_gamma] & n_0 == n_0_list[k_n0] & ratio_inter == ratio_inter_list[k_inter] & p_Z==p_Z_list[k_pz] & inter_strength==inter_strength_list[k_inter_strength] & alpha_modifier == alpha_list[k_alpha] )
                  condition8 <- (heter_level == heter_level_list[k_heter] & gamma_multi == gamma_list[k_gamma] & n_0 == n_0_list[k_n0] & ratio_inter == ratio_inter_list[k_inter] & p_Z==p_Z_list[k_pz] & inter_strength==inter_strength_list[k_inter_strength] & eta_modifier == eta_list[k_eta])
                  
                  if  ( !(condition1|condition2|condition3|condition4|condition5|condition6|condition7|condition8) ){
                    next
                  }else{
                    xunhuanshijiweizhi <- xunhuanshijiweizhi+1
                    jishi2 <- Sys.time()
                    xunhuanyongshi <- difftime(jishi2,jishi1,units = 'hours')
                    jishi1 <- jishi2
                    print(paste('当前循环位置：',(xunhuanshijiweizhi-1),'/',zongxunhuancishu,'，','本次循环用时：',as.character(round(xunhuanyongshi,5)),'小时，当前时间：',as.character(Sys.time())))
                  }
                  #===========================================================================
                  ## 结果保存
                  param_df2 <- data.frame(matrix(rep(c(ite_times,p_Z,gamma_multi,heter_level,n_0,n_1,ratio_inter,inter_strength,alpha_modifier,eta_modifier),len_methods*3),nrow = len_methods*3,byrow = TRUE))
                  param_df2 <- cbind(param_df2,rep(method_names_ordered,3))
                  param_df2 <- cbind(param_df2,rep(index_ordered,each=len_methods))
                  
                  
                  param_df <- data.frame(matrix(rep(c(p_Z,gamma_multi,heter_level,n_0,n_1,ratio_inter,inter_strength,alpha_modifier,eta_modifier),len_methods*ite_times),nrow = len_methods*ite_times,byrow = TRUE))
                  param_df <- cbind(rep(c(1:ite_times),each=len_methods),param_df)
                  param_df <- cbind(param_df,rep(method_names_ordered,ite_times))
                  #----------------------------
                  beta_hat_ite <- foreach(ite_k=1:ite_times, .export =c('heter_level', 'data_generate_Y', 'data_generate_X', 'Twosamle_package', 'Trans_OLS2', 'IPW2', 'AIPW2')
                                          , .packages = c('fixest','TwoSampleMR','glmnet','mr.raps','penalized','ridge'))%dopar%
                    {
                      #for(ite_k in 1:ite_times){
                      tryCatch({
                        indicator_print <- 0
                        
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
                        conditions_met <- FALSE
                        while(!conditions_met){
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
                        #C <- rbind(C_0,C_1)
                        
                        #Twosamle_package2(p_Z,p_C,data_0,data_1)
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
                            if(p_Q_statistics_group <0.05){
                              modifier_matrix[k,j] <- 1
                            }else{
                              modifier_matrix[k,j] <- 0
                            }
                          }
                        }
                        # chow_res <- strucchange::Fstats(X~Z.1+Z.2+Z.3,data = data_chow_test,from = 5996)
                        # sctest(chow_res)
                        
                        
                        #----------------------------
                        #----------------------------
                        b_Y_list2 <- rep(0,p_Z)
                        sd_Y_list2 <- rep(0,p_Z)
                        b_Y_list3 <- rep(0,p_Z)
                        sd_Y_list3 <- rep(0,p_Z)
                        b_Y_list4 <- rep(0,p_Z)
                        sd_Y_list4 <- rep(0,p_Z)
                        #===============================================================Twosample-MR包     1
                        #-------------------
                        beta_hat_1_list <- Twosamle_package(p_Z=p_Z,p_C=p_C,data_0=data_0,data_1=data_1)
                        beta_hat_1 <- beta_hat_1_list$b
                        
                        #===============================================================Trans_OR
                        X_hat_matrix_OR <- Trans_OLS2(Z_0=Z_0,Z_1=Z_1,C_0=C_0,C_1=C_1,X_1=X_1)$X_hat_matrix
                        #===============================================================Trans_OR0         21
                        for (i in 1:p_Z) {  
                          
                          model2 <- lm(Y_0~X_hat_matrix_OR[,i] + C_0+Z_0[,i])  #有Z
                          summary_model2 <- summary(model2)
                          pval_model2 <- summary_model2$coefficients[,4]
                          Z_pval <- pval_model2[length(pval_model2)]
                          
                          pearson_r <- abs(cor(X_hat_matrix_OR[,i],Z_0[,i]))
                          
                          # vif_df2 <- car::vif(model2)
                          # vif_2_vec <- vif_df2[,1]
                          if(indicator_print == 1){
                            print(paste0(Z_pval,', ',pearson_r))
                          }
                          #print(vif_2_vec)
                          
                          if( (Z_pval < 0.05) & (pearson_r < 0.7) ){    #有水平多效性：岭回归
                            model_ridge2 <- ridge::linearRidge(Y_0~X_hat_matrix_OR[,i] + C_0+Z_0[,i])
                            summary_ridge2 <- summary(model_ridge2)
                            summary_chosen2 <- summary_ridge2$summaries[[summary_ridge2$chosen.nPCs]]   # 提取选择的 nPCs 对应的 summary
                            coefficients_chosen2 <- summary_chosen2$coefficients   # 从 chosen_summary 提取 coefficients
                            coef_chosen2 <- coefficients_chosen2[,1]
                            se_chosen2 <- coefficients_chosen2[,3][-1]/model_ridge2$scales
                            
                            b_Y_list2[i] <-coef_chosen2[2]
                            sd_Y_list2[i] <- se_chosen2[1]
                          }else{   #无水平多效性：线性回归
                            model2_1 <- lm(Y_0~X_hat_matrix_OR[,i] + C_0)  #无Z
                            b_Y_list2[i] <- model2_1$coefficients[2]
                            sd_Y_list2[i] <- summary(model2_1)$coefficients[,'Std. Error'][2]
                          }
                          
                        }
                        # cor(X_hat_matrix_OR[,1],Z_0[,1])
                        sd_Y_list2_inv <- 1/sd_Y_list2 
                        beta_hat_2 <- sum(sd_Y_list2_inv^2 * b_Y_list2)/sum(sd_Y_list2_inv^2)
                        #===============================================================Trans_IPW
                        X_hat_matrix_IPW <- IPW2(Z_0=Z_0,Z_1=Z_1,C_0=C_0,C_1=C_1,X_1=X_1)$X_hat_matrix
                        #===============================================================Trans_IPW0          31
                        for (i in 1:p_Z) {
                          
                          model3 <- lm(Y_0~X_hat_matrix_IPW[,i] + C_0+Z_0[,i])  #有Z
                          summary_model3 <- summary(model3)
                          pval_model3 <- summary_model3$coefficients[,4]
                          Z_pval <- pval_model3[length(pval_model3)]
                          
                          pearson_r <- cor(X_hat_matrix_IPW[,i],Z_0[,i])
                          if(indicator_print == 1){
                            print(paste0(Z_pval,', ',pearson_r))
                          }
                          if( (Z_pval < 0.05) & (pearson_r < 0.7) ){    #有水平多效性：岭回归
                            model_ridge3 <- ridge::linearRidge(Y_0~X_hat_matrix_IPW[,i] + C_0+Z_0[,i])
                            summary_ridge3 <- summary(model_ridge3)
                            summary_chosen3 <- summary_ridge3$summaries[[summary_ridge3$chosen.nPCs]]   # 提取选择的 nPCs 对应的 summary
                            coefficients_chosen3 <- summary_chosen3$coefficients   # 从 chosen_summary 提取 coefficients
                            coef_chosen3 <- coefficients_chosen3[,1]
                            se_chosen3 <- coefficients_chosen3[,3][-1]/model_ridge3$scales
                            
                            b_Y_list3[i] <-coef_chosen3[2]
                            sd_Y_list3[i] <- se_chosen3[1]
                          }else{   #无水平多效性：线性回归
                            model3_1 <- lm(Y_0~X_hat_matrix_IPW[,i] + C_0)  #无Z
                            b_Y_list3[i] <- model3_1$coefficients[2]
                            sd_Y_list3[i] <- summary(model3_1)$coefficients[,'Std. Error'][2]
                          }
                          
                        }
                        sd_Y_list3_inv <- 1/sd_Y_list3 
                        beta_hat_3 <- sum(sd_Y_list3_inv^2 * b_Y_list3)/sum(sd_Y_list3_inv^2)
                        #===============================================================Trans_AIPW
                        X_hat_matrix_AIPW <- AIPW2(Z_0=Z_0,Z_1=Z_1,C_0=C_0,C_1=C_1,X_1=X_1)$X_hat_matrix
                        #===============================================================Trans_AIPW0          41
                        for (i in 1:p_Z) {
                          model4 <- lm(Y_0~X_hat_matrix_AIPW[,i] + C_0+Z_0[,i])  #有Z
                          summary_model4 <- summary(model4)
                          pval_model4 <- summary_model4$coefficients[,4]
                          Z_pval <- pval_model4[length(pval_model4)]
                          
                          pearson_r <- cor(X_hat_matrix_AIPW[,i],Z_0[,i])
                          if(indicator_print == 1){
                            print(paste0(Z_pval,', ',pearson_r))
                          }
                          if( (Z_pval < 0.05) & (pearson_r < 0.7) ){    #有水平多效性：岭回归
                            model_ridge4 <- ridge::linearRidge(Y_0~X_hat_matrix_AIPW[,i] + C_0+Z_0[,i])
                            summary_ridge4 <- summary(model_ridge4)
                            summary_chosen4 <- summary_ridge4$summaries[[summary_ridge4$chosen.nPCs]]   # 提取选择的 nPCs 对应的 summary
                            coefficients_chosen4 <- summary_chosen4$coefficients   # 从 chosen_summary 提取 coefficients
                            coef_chosen4 <- coefficients_chosen4[,1]
                            se_chosen4 <- coefficients_chosen4[,3][-1]/model_ridge4$scales
                            
                            b_Y_list4[i] <-coef_chosen4[2]
                            sd_Y_list4[i] <- se_chosen4[1]
                          }else{   #无水平多效性：线性回归
                            model4_1 <- lm(Y_0~X_hat_matrix_AIPW[,i] + C_0)  #无Z
                            b_Y_list4[i] <- model4_1$coefficients[2]
                            sd_Y_list4[i] <- summary(model4_1)$coefficients[,'Std. Error'][2]
                          }
                        }
                        sd_Y_list4_inv <- 1/sd_Y_list4
                        beta_hat_4 <- sum(sd_Y_list4_inv^2 * b_Y_list4)/sum(sd_Y_list4_inv^2)
                        #----------------------------结果保存
                        beta_hat_vec <- c(beta_hat_1, beta_hat_2, beta_hat_3,beta_hat_4)
                        p_vec <- c(beta_hat_1_list$p,rep(NA,len_methods-8))
                        p_cover_vec <- c(beta_hat_1_list$p_cover,rep(NA,len_methods-8)) 
                        
                        return_list <- list(beta=beta_hat_vec,p=p_vec,p_cover=p_cover_vec,seed=seed)
                        #print(return_list)
                        return(return_list)  #还要返回p值
                        
                      }, error = function(e) {
                        message(paste("Iteration", ite_k, "failed:", e$message))
                        return(NULL)
                      })
                      
                    }
                  
                  # 结果保存
                  #+++++++++++++++ 1
                  ite_return <- do.call(rbind,beta_hat_ite)
                  ite_df2 <- do.call(rbind,ite_return[,1])  #beta值矩阵
                  #ite_p <- do.call(rbind,ite_return[,2])   #p值矩阵
                  #ite_p_cover <- do.call(rbind,ite_return[,3])  #p值矩阵
                  
                  beta_hat_ite_vec <- as.vector(unlist(ite_return[,1]))   # (len_methods个beta_hat) * (循环次数)  个beta_hat值
                  p_ite_vec <- as.vector(unlist(ite_return[,2]))    # (len_methods个beta_hat) * (循环次数)  个p值
                  p_cover_ite_vec <- as.vector(unlist(ite_return[,3]))   # (len_methods个beta_hat) * (循环次数)  个p值
                  
                  save_df_ite <- cbind(param_df,beta_hat_ite_vec,p_ite_vec,p_cover_ite_vec)
                  save_df <- rbind(save_df,save_df_ite)
                  
                  #+++++++++++++++ 2
                  
                  bias_vec <- as.vector(apply(ite_df2-beta, 2, function(x) mean(x[abs(x) < 10], na.rm = TRUE)))
                  sd_vec <- as.vector(apply(ite_df2, 2, function(x) sd(x[abs(x) < 10], na.rm = TRUE)))
                  mse_vec <- bias_vec^2 + sd_vec^2
                  save_df_ite2 <- cbind(param_df2,c(bias_vec,sd_vec,mse_vec))
                  save_df2 <- rbind(save_df2,save_df_ite2)
                  #+++++++++++++++ 3
                  if(condition3){
                    test_df <- save_df_ite
                    colnames(test_df) <- col_names
                    sd_vec_my_methods <- sd_vec[(len_other_methods+1):len_methods]
                    test_df <- test_df %>%
                      mutate(
                        # 如果method在my_methods中，计算p值
                        p_power = ifelse(method %in% my_methods, 
                                         2 * (1 - pnorm(abs(beta_hat / sd_vec_my_methods[match(method, my_methods)]))),
                                         p_power),
                        # 如果method在my_methods中，计算p_power
                        p_cover = ifelse(method %in% my_methods, 
                                         2 * (1 - pnorm(abs((beta_hat - 1) / sd_vec_my_methods[match(method, my_methods)]))),
                                         p_cover)
                      )
                    
                    indice_power_list <- which(n_0 == n_0_list)
                    save_power_list[[indice_power_list]] <- test_df
                  }
                  # sd_vec_my_methods <- sd_vec[(len_other_methods+1):len_methods]
                  # test_df <- test_df %>%
                  #   mutate(
                  #     # 如果method在my_methods中，计算p值
                  #     p_power = ifelse(method %in% my_methods, 
                  #                  2 * (1 - pnorm(abs(beta_hat / sd_vec_my_methods[match(method, my_methods)]))),
                  #                  p_power),
                  #     # 如果method在my_methods中，计算p_power
                  #     p_cover = ifelse(method %in% my_methods, 
                  #                      2 * (1 - pnorm(abs((beta_hat - 1) / sd_vec_my_methods[match(method, my_methods)]))),
                  #                      p_cover)
                  #   )
                  
                  
                  #2 * (1 - pnorm(abs(res_list$beta/res_list$se)))
                  
                  #+++++++++++++++
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
write.csv(save_df,file=paste(path_output_tables, '/data_all_',date_mark,'_',ver,'.csv',sep = ''))
write.csv(save_df2,file=paste(path_output_tables, '/statistics_all_',date_mark,'_',ver,'.csv',sep = ''))
#===============================================================================绘图1
# time_end <- Sys.time()
# time_consume <- difftime(time_end,time_start)
# time_consume 