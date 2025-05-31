# 迭代留一法主函数（针对元分析）
iterative_loo_meta <- function(data, yi_col, vi_col, max_iter = 5, i2_threshold = 25) {
  # 参数说明:
  # data: 输入数据集（需包含效应量(yi)和方差(vi)列）
  # yi_col: 效应量列名（字符类型）
  # vi_col: 方差列名（字符类型）
  # max_iter: 最大迭代次数（默认5）
  # i2_threshold: I²异质性容忍阈值（默认25%，即低于25%时停止）
  
  current_data <- data              # 当前数据集（动态更新）
  removed_studies <- character(0)   # 记录被剔除的研究
  results_history <- list()         # 记录每次迭代的元分析结果
  
  for (iter in 1:max_iter) {
    # ---------- 步骤1: 执行留一交叉验证 ----------
    loo_results <- leave1out(
      rma(
        yi = get(yi_col), 
        vi = get(vi_col), 
        data = current_data, 
        method = "REML"  # 使用限制性最大似然估计
      )
    )
    
    # ---------- 步骤2: 计算每个研究剔除后的I²值 ----------
    i2_values <- loo_results$I2
    max_i2_idx <- which.min(i2_values)  # 找到剔除后I²最大的研究
    
    # ---------- 步骤4: 剔除对异质性贡献最大的研究 ----------
    removed_studies <- c(removed_studies, rownames(current_data)[max_i2_idx])
    results_history[[iter]] <- list(
      model = loo_results,
      removed = rownames(current_data)[max_i2_idx]
    )
    
    current_data <- current_data[-max_i2_idx, ]
    # ---------- 步骤3: 记录结果并检查终止条件 ----------
    current_model <- #leave1out(
      rma(
        yi = get(yi_col), 
        vi = get(vi_col), 
        data = current_data, 
        method = "REML"
      )
    #)
    
    # 如果当前I²已低于阈值，提前终止
    if (all(current_model$I2 < i2_threshold)) break
    
  }
  
  # ---------- 返回结果 ----------
  return_list <-   list(
    final_model = current_model,
    removed_studies = removed_studies,
    history = results_history
  )
  
  return(return_list)
}
#--------------------
#time_start <- Sys.time()
# 设置并行计算的核心数量
cores <- detectCores()
registerDoParallel(cores = cores-4)

#参数空间
pi_data1 <- matrix(round(with_seed(1, runif(90000,0,2)),2),nrow = 300)
pi_data2 <- matrix(round(with_seed(1, runif(90000,0,2)),2),nrow = 300)
alpha_data <- round(with_seed(1,runif(1000,0.2,0.8)),2)
gamma_data <- -round(with_seed(1,runif(1000,0.2,0.4)),2)
eta1_data <- round(with_seed(1,runif(1000,-0.5,0.5)),2)
eta2_data <- round(with_seed(1,runif(1000,-0.5,0.5)),2)
#===================================================================参数设置
beta <- 0.4

n_0 <- 20000
n_02 <- 20000
n_1 <- n_0

heter_level <- 1
ratio_inter <- 1
gamma_multi <- 1

inter_strength <- 1
alpha_modifier <- 1
eta_modifier <- 1

ratio_pleiotropy <- 1

Intercept <- -2
#--------------------------
p_Z <- 3
p_Z_max =3
#-----------
p_C <- 2

eta_1 <- eta1_data[1:p_C] * eta_modifier
eta_2 <- eta2_data[1:p_C] * eta_modifier

#-----
#alpha <- seq(2,4,2/(p_Z-1))*alpha_modifier
alpha <- alpha_data[1:p_Z_max] * alpha_modifier
pi <- pi_data1[1:p_Z_max,1:p_C] * inter_strength
#pi <- c(0.53,1.35)
pi2 <- pi_data2[1:p_Z_max,1:p_C] * inter_strength
#omega <- omega_data[1:p_C]

gamma <- gamma_data[1:p_Z_max] * gamma_multi
#gamma <- 0.13
gamma <- gamma_data[1:p_Z_max] * gamma_multi #/ p_Z_max
indice_gamma1 <- seq(1,p_Z_max,3)
num_gamma1 <- round(p_Z_max/3)
indice_gamma1 <- indice_gamma1[1:num_gamma1]
indice_gamma0 <- setdiff(1:p_Z_max,indice_gamma1)

gamma[indice_gamma0] <- 0


alpha_true <- alpha_data[1:p_Z] * alpha_modifier
pi_true <- pi_data1[1:p_Z,1:p_C] * inter_strength
pi2_true <- pi_data2[1:p_Z,1:p_C] * inter_strength
gamma_true <- gamma_data[1:p_Z] * gamma_multi



#===============================================================================
gamma_choose <- seq(-1,1,0.1)
plot_list <- list()
for (g in 1:length(gamma_choose)) {
  registerDoParallel(cores = cores-4)
  gamma <- c(gamma_choose[g],0,0)
  
  #----
  param_X <- list(alpha=alpha,pi=pi,pi2=pi2,eta_1=eta_1
                  #,U_md = U_md
  )
  # param_Y <- list(beta=beta,eta_2=eta_2,gamma=gamma,ratio_pleiotropy=ratio_pleiotropy
  #                 ,heter_level=heter_level,omega=omega,Intercept=Intercept
  #                 #,U_md = U_md
  # )
  param_Y <- list(beta=beta,eta_2=eta_2,gamma=gamma,ratio_pleiotropy=ratio_pleiotropy
                  ,heter_level=heter_level,Intercept=Intercept
                  #,U_md = U_md
  )
  #heritability <- 0.7
  pval_i <- 0.05
  #0.05/p_Z
  
  #ite_times <- 1000
  #===========================================================================
  ## 结果保存
  other_methods <- c("IVW","Egger","Weighted median","Weighted mode","MR-RAPS","MR-Conmix", "MR-Robust", "MR-Lasso")
  my_methods <- c('TLMR-OR'
                  ,'TLMR-IPW'
                  ,'TLMR-AIPW')
  
  method_names <- c(other_methods,my_methods)
  
  methods_ordered <- factor(method_names,levels = method_names)
  len_methods <- length(methods_ordered)
  index_ordered <- factor(c('Bias','SD','MSE'),levels = c('Bias','SD','MSE'))
  
  param_df2 <- data.frame(matrix(rep(c(ite_times,beta,p_Z,gamma_multi,heter_level,n_0,n_1,ratio_inter,inter_strength,alpha_modifier,eta_modifier),len_methods*3),nrow = len_methods*3,byrow = TRUE))
  param_df2 <- cbind(param_df2,rep(methods_ordered,3))
  param_df2 <- cbind(param_df2,rep(index_ordered,each=len_methods))
  
  
  param_df <- data.frame(matrix(rep(c(beta,p_Z,gamma_multi,heter_level,n_0,n_1,ratio_inter,inter_strength,alpha_modifier,eta_modifier),len_methods*ite_times),nrow = len_methods*ite_times,byrow = TRUE))
  param_df <- cbind(rep(c(1:ite_times),each=len_methods),param_df)
  param_df <- cbind(param_df,rep(methods_ordered,ite_times))
  
  
  col_names <- c('iteration time','beta','p_Z', 'gamma','heterogeneity','n_0','n_1','inter_ratio','inter_strength','alpha_modifier','eta_modifier','method','beta_hat','p_power','p_cover')
  
  col_names2 <- c('max iteration times','beta','p_Z', 'gamma','heterogeneity','n_0','n_1','inter_ratio','inter_strength','alpha_modifier','eta_modifier','method','index','num')
  
  
  # > beta_hat_vec
  # [1] 0.3907408 0.3696579 0.3281686 0.3100789 0.3454174 0.5404447 0.3884075 0.3907408
  # [9] 0.4635852 0.4530766 0.4530766
  #===========================================================================
  beta_hat_ite <- foreach(ite_k=1:ite_times, .export =c('heter_level')
                          , .packages = c('fixest','TwoSampleMR','glmnet','mr.raps','penalized','ridge','FNN','data.table','dplyr','metafor'))%dopar%
    {
      #for(ite_k in 1:ite_times){
      tryCatch({
        #indicator_print <- 1
        
        seed <- ite_k * 100
        set.seed(seed)
        #-----------------------------------------------生成随机数
        data_0 <- data_generate_Y_new2(n=n_0,p_Z=p_Z,p_Z_max = p_Z_max
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
            data_0 <- data_generate_Y_new2(n=n_0,p_Z=p_Z,p_Z_max = p_Z_max
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
        
        #===============================================================Twosample-MR包     1
        beta_hat_1_list <- Twosamle_package(p_Z=p_Z,p_C=p_C,data_0=data_0,data_1=data_1,beta=beta)
        beta_hat_1 <- beta_hat_1_list$beta
        #===============================================================Trans-MR
        #----------------------------看是否存在修饰作用的代码
        modifier_matrix <- matrix(nrow = p_Z,ncol = p_C)
        p_modifier_matrix <- matrix(nrow = p_Z,ncol = p_C)
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
            res_meta_analysis <- metafor::rma.uni(yi=beta_group_vec,sei = se_group_vec,method = 'REML') # method = 'DL'
            Q_statistics_group <- res_meta_analysis$QE
            p_Q_statistics_group <- res_meta_analysis$QEp
            p_modifier_matrix[k,j] <- p_Q_statistics_group
            if(p_Q_statistics_group <(0.05)){
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
          row_degrade <- c(beta_hat_degrade,se_degrade)
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
            X_hat_matrix_OR <- Trans_OLS2(Z_0=Z_0m,Z_1=Z_1m,C_0=C_0,C_1=C_1,X_1=X_1)$X_hat_matrix
            #Trans_IPW
            X_hat_matrix_IPW <- IPW2(Z_0=Z_0m,Z_1=Z_1m,C_0=C_0,C_1=C_1,X_1=X_1)$X_hat_matrix
            #Trans_AIPW
            X_hat_matrix_AIPW <- AIPW2(Z_0=Z_0m,Z_1=Z_1m,C_0=C_0,C_1=C_1,X_1=X_1
            )$X_hat_matrix
            #-----------------------------------------
            
            #---------------------------------------------------------------Trans_OR0         21
            #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            if(p_Zm >= 1){
              
              for (i in 1:p_Zm) {
                
                model2 <- glm(Y_0~X_hat_matrix_OR[,i] + C_0, family = binomial(link = "log"))  #有Z
                b_Y_list2[i] <- model2$coefficients[2]
                sd_Y_list2[i] <- summary(model2)$coefficients[,'Std. Error'][2]
                
              }
              
              data_2 <- data.frame(b_Y = b_Y_list2,
                                   v_Y = sd_Y_list2^2)
              I2 <- rma(
                yi = b_Y,#yi_col, 
                vi = v_Y,#vi_col, 
                data = data_2, 
                method = "REML"  # 使用限制性最大似然估计
              )$I2
              
              if(I2 >= 25 ){
                ite_loo <- iterative_loo_meta(data= data_2, yi_col = 'b_Y', vi_col = 'v_Y')
                indice_removed_data <- as.numeric(ite_loo$removed_studies)
                indice_remainde_data <- setdiff(1:p_Zm,indice_removed_data) #看情况可以删
                
                removed_data <- data_2[indice_removed_data,]
                remained_data <- data_2[-indice_removed_data,]
                
                
                for(m in indice_removed_data){
                  model2 <- glm(Y_0~X_hat_matrix_OR[,m] + C_0 + Z_0m[,m], family = binomial(link = "log"))  #有Z
                  
                  data_2[m,1] <- model2$coefficients[2]
                  data_2[m,2] <- summary(model2)$coefficients[,'Std. Error'][2]
                  
                }
                
              }
            }
            data_2 <- rbind(data_2,row_degrade)
            
            sd_Y_list2_inv <- sqrt(1/data_2[,2])
            b_Y_list2 <- data_2[,1]
            beta_hat_2 <- sum(sd_Y_list2_inv^2 * b_Y_list2)/sum(sd_Y_list2_inv^2)
            #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            #---------------------------------------------------------------Trans_IPW0          31
            #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            if(p_Zm>=1){
              
              for (i in 1:p_Zm) {
                
                model3 <- glm(Y_0~X_hat_matrix_IPW[,i] + C_0, family = binomial(link = "log"))  #有Z
                b_Y_list3[i] <- model3$coefficients[2]
                sd_Y_list3[i] <- summary(model3)$coefficients[,'Std. Error'][2]
                
              }
              # sd_Y_list3_inv <- 1/sd_Y_list3
              # beta_hat_3 <- sum(sd_Y_list3_inv^2 * b_Y_list3)/sum(sd_Y_list3_inv^2)
              data_3 <- data.frame(b_Y = b_Y_list3,
                                   v_Y = sd_Y_list3^2)
              I2 <- rma(
                yi = b_Y,#yi_col, 
                vi = v_Y,#vi_col, 
                data = data_3, 
                method = "REML"  # 使用限制性最大似然估计
              )$I2
              
              if(I2 >= 25 ){
                ite_loo <- iterative_loo_meta(data= data_3, yi_col = 'b_Y', vi_col = 'v_Y')
                indice_removed_data <- as.numeric(ite_loo$removed_studies)
                indice_remainde_data <- setdiff(1:p_Zm,indice_removed_data) #看情况可以删
                
                removed_data <- data_3[indice_removed_data,]
                remained_data <- data_3[-indice_removed_data,]
                
                
                for(m in indice_removed_data){
                  model3 <- glm(Y_0~X_hat_matrix_IPW[,m] + C_0 + Z_0m[,m], family = binomial(link = "log"))  #有Z
                  
                  data_3[m,1] <- model3$coefficients[2]
                  data_3[m,2] <- summary(model3)$coefficients[,'Std. Error'][2]
                  
                }
                
              }
            }
            data_3 <- rbind(data_3,row_degrade)
            
            sd_Y_list3_inv <- sqrt(1/data_3[,2])
            b_Y_list3 <- data_3[,1]
            beta_hat_3 <- sum(sd_Y_list3_inv^2 * b_Y_list3)/sum(sd_Y_list3_inv^2)
            #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            #---------------------------------------------------------------Trans_AIPW0          41
            #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            if(p_Zm>=1){
              for (i in 1:p_Zm) {
                
                model4 <- glm(Y_0~X_hat_matrix_AIPW[,i] + C_0, family = binomial(link = "log"))  #有Z
                b_Y_list4[i] <- model4$coefficients[2]
                sd_Y_list4[i] <- summary(model4)$coefficients[,'Std. Error'][2]
                
              }
              
              data_4 <- data.frame(b_Y = b_Y_list4,
                                   v_Y = sd_Y_list4^2)
              I2 <- rma(
                yi = b_Y,#yi_col, 
                vi = v_Y,#vi_col, 
                data = data_4, 
                method = "REML"  # 使用限制性最大似然估计
              )$I2
              
              if(I2 >= 25 ){
                ite_loo <- iterative_loo_meta(data= data_4, yi_col = 'b_Y', vi_col = 'v_Y')
                indice_removed_data <- as.numeric(ite_loo$removed_studies)
                indice_remainde_data <- setdiff(1:p_Zm,indice_removed_data) #看情况可以删
                
                removed_data <- data_4[indice_removed_data,]
                remained_data <- data_4[-indice_removed_data,]
                
                
                for(m in indice_removed_data){
                  model4 <- glm(Y_0~X_hat_matrix_AIPW[,m] + C_0 + Z_0m[,m], family = binomial(link = "log"))  #有Z
                  
                  data_4[m,1] <- model4$coefficients[2]
                  data_4[m,2] <- summary(model4)$coefficients[,'Std. Error'][2]
                  
                }
                
              }
            }
            data_4 <- rbind(data_4,row_degrade)
            
            sd_Y_list4_inv <- sqrt(1/data_4[,2])
            b_Y_list4 <- data_4[,1]
            beta_hat_4 <- sum(sd_Y_list4_inv^2 * b_Y_list4)/sum(sd_Y_list4_inv^2)
            #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          }
        }else{
          #===============================================================迁移估计量
          #Trans_OR
          X_hat_matrix_OR <- Trans_OLS2(Z_0=Z_0,Z_1=Z_1,C_0=C_0,C_1=C_1,X_1=X_1)$X_hat_matrix
          #Trans_IPW
          X_hat_matrix_IPW <- IPW2(Z_0=Z_0,Z_1=Z_1,C_0=C_0,C_1=C_1,X_1=X_1)$X_hat_matrix
          #Trans_AIPW
          X_hat_matrix_AIPW <- AIPW2(Z_0=Z_0,Z_1=Z_1,C_0=C_0,C_1=C_1,X_1=X_1
          )$X_hat_matrix
          #---------------------------------------------------------------Trans_OR0         21
          #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          for (i in 1:p_Z) {
            
            model2 <- glm(Y_0~X_hat_matrix_OR[,i] + C_0, family = binomial(link = "log"))  #有Z
            b_Y_list2[i] <- model2$coefficients[2]
            sd_Y_list2[i] <- summary(model2)$coefficients[,'Std. Error'][2]
            
          }
          
          data_2 <- data.frame(b_Y = b_Y_list2,
                               v_Y = sd_Y_list2^2)
          I2 <- rma(
            yi = b_Y,#yi_col, 
            vi = v_Y,#vi_col, 
            data = data_2, 
            method = "REML"  # 使用限制性最大似然估计
          )$I2
          
          if(I2 >= 25 ){
            ite_loo <- iterative_loo_meta(data= data_2, yi_col = 'b_Y', vi_col = 'v_Y')
            indice_removed_data <- as.numeric(ite_loo$removed_studies)
            indice_remainde_data <- setdiff(1:p_Z,indice_removed_data) #看情况可以删
            
            removed_data <- data_2[indice_removed_data,]
            remained_data <- data_2[-indice_removed_data,]
            
            
            for(m in indice_removed_data){
              model2 <- glm(Y_0~X_hat_matrix_OR[,m] + C_0 + Z_0[,m], family = binomial(link = "log"))  #有Z
              
              data_2[m,1] <- model2$coefficients[2]
              data_2[m,2] <- summary(model2)$coefficients[,'Std. Error'][2]
              
            }
            
          }
          
          sd_Y_list2_inv <- sqrt(1/data_2[,2])
          b_Y_list2 <- data_2[,1]
          beta_hat_2 <- sum(sd_Y_list2_inv^2 * b_Y_list2)/sum(sd_Y_list2_inv^2)
          #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          #---------------------------------------------------------------Trans_IPW0          31
          #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          for (i in 1:p_Z) {
            
            model3 <- glm(Y_0~X_hat_matrix_IPW[,i] + C_0, family = binomial(link = "log"))  #有Z
            b_Y_list3[i] <- model3$coefficients[2]
            sd_Y_list3[i] <- summary(model3)$coefficients[,'Std. Error'][2]
            
          }
          # sd_Y_list3_inv <- 1/sd_Y_list3
          # beta_hat_3 <- sum(sd_Y_list3_inv^2 * b_Y_list3)/sum(sd_Y_list3_inv^2)
          data_3 <- data.frame(b_Y = b_Y_list3,
                               v_Y = sd_Y_list3^2)
          I2 <- rma(
            yi = b_Y,#yi_col, 
            vi = v_Y,#vi_col, 
            data = data_3, 
            method = "REML"  # 使用限制性最大似然估计
          )$I2
          
          if(I2 >= 25 ){
            ite_loo <- iterative_loo_meta(data= data_3, yi_col = 'b_Y', vi_col = 'v_Y')
            indice_removed_data <- as.numeric(ite_loo$removed_studies)
            indice_remainde_data <- setdiff(1:p_Z,indice_removed_data) #看情况可以删
            
            removed_data <- data_3[indice_removed_data,]
            remained_data <- data_3[-indice_removed_data,]
            
            
            for(m in indice_removed_data){
              model3 <- glm(Y_0~X_hat_matrix_IPW[,m] + C_0 + Z_0[,m], family = binomial(link = "log"))  #有Z

              data_3[m,1] <- model3$coefficients[2]
              data_3[m,2] <- summary(model3)$coefficients[,'Std. Error'][2]
              
            }
            
          }
          
          sd_Y_list3_inv <- sqrt(1/data_3[,2])
          b_Y_list3 <- data_3[,1]
          beta_hat_3 <- sum(sd_Y_list3_inv^2 * b_Y_list3)/sum(sd_Y_list3_inv^2)
          #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          #---------------------------------------------------------------Trans_AIPW0          41
          #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          for (i in 1:p_Z) {
            
            model4 <- glm(Y_0~X_hat_matrix_AIPW[,i] + C_0, family = binomial(link = "log"))  #有Z
            b_Y_list4[i] <- model4$coefficients[2]
            sd_Y_list4[i] <- summary(model4)$coefficients[,'Std. Error'][2]
            
          }
          
          data_4 <- data.frame(b_Y = b_Y_list4,
                               v_Y = sd_Y_list4^2)
          I2 <- rma(
            yi = b_Y,#yi_col, 
            vi = v_Y,#vi_col, 
            data = data_4, 
            method = "REML"  # 使用限制性最大似然估计
          )$I2
          
          if(I2 >= 25 ){
            ite_loo <- iterative_loo_meta(data= data_4, yi_col = 'b_Y', vi_col = 'v_Y')
            indice_removed_data <- as.numeric(ite_loo$removed_studies)
            indice_remainde_data <- setdiff(1:p_Z,indice_removed_data) #看情况可以删
            
            removed_data <- data_4[indice_removed_data,]
            remained_data <- data_4[-indice_removed_data,]
            
            
            for(m in indice_removed_data){
              model4 <- glm(Y_0~X_hat_matrix_AIPW[,m] + C_0 + Z_0[,m], family = binomial(link = "log"))  #有Z
              
              data_4[m,1] <- model4$coefficients[2]
              data_4[m,2] <- summary(model4)$coefficients[,'Std. Error'][2]
              
            }
            
          }
          
          sd_Y_list4_inv <- sqrt(1/data_4[,2])
          b_Y_list4 <- data_4[,1]
          beta_hat_4 <- sum(sd_Y_list4_inv^2 * b_Y_list4)/sum(sd_Y_list4_inv^2)
          #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
       
        }
        
        #----------------------------结果保存
        beta_hat_vec <- c(beta_hat_1, beta_hat_2, beta_hat_3,beta_hat_4)
        #--------------------------------------------------IVW方差计算
        sd_OR <- sqrt((1/beta_hat_2^2) * 1/sum(sd_Y_list2_inv^2))
        sd_IPW <- sqrt((1/beta_hat_3^2) * 1/sum(sd_Y_list3_inv^2))
        sd_AIPW <- sqrt((1/beta_hat_4^2) * 1/sum(sd_Y_list4_inv^2))
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
        # #=========================================================================
        # #Trans_OR
        # X_hat_matrix_OR <- Trans_OLS2(Z_0=Z_0,Z_1=Z_1,C_0=C_0,C_1=C_1,X_1=X_1)$X_hat_matrix
        # #Trans_IPW
        # X_hat_matrix_IPW <- IPW2(Z_0=Z_0,Z_1=Z_1,C_0=C_0,C_1=C_1,X_1=X_1)$X_hat_matrix
        # #Trans_AIPW
        # X_hat_matrix_AIPW <- AIPW2(Z_0=Z_0,Z_1=Z_1,C_0=C_0,C_1=C_1,X_1=X_1
        # )$X_hat_matrix
        # 
        # # model2_1 <- glm(Y_0~X_hat_matrix_OR + C_0, family = binomial(link = "log"))
        # # model2_2 <- glm(Y_0~X_hat_matrix_IPW + C_0, family = binomial(link = "log"))
        # # model2_3 <- glm(Y_0~X_hat_matrix_AIPW + C_0, family = binomial(link = "log"))
        # 
        # model2_1 <- glm(Y_0~X_hat_matrix_OR + C_0+ Z_0, family = binomial(link = "log"))
        # model2_2 <- glm(Y_0~X_hat_matrix_IPW + C_0+ Z_0, family = binomial(link = "log"))
        # model2_3 <- glm(Y_0~X_hat_matrix_AIPW + C_0+ Z_0, family = binomial(link = "log"))
        # 
        # # model2_1 <- glm(Y_0~X_hat_matrix_OR + C_0, family = binomial(link = "log"))
        # # model2_2 <- glm(Y_0~X_hat_matrix_IPW + C_0, family = binomial(link = "log"))
        # # model2_3 <- glm(Y_0~X_hat_matrix_AIPW + C_0, family = binomial(link = "log"))
        # 
        # # model2_1
        # # model2_2
        # 
        # beta_1 <- model2_1$coefficients[2]
        # beta_2 <- model2_2$coefficients[2]
        # beta_3 <- model2_3$coefficients[2]
        # 
        # #omega_hat <- c(omega_1,omega_2,omega_3)
        # 
        # return_list <- list(beta1=beta_1,beta2=beta_2,beta3=beta_3)
        # #=========================================================================
      }, error = function(e) {
        #warnings(paste("Iteration", ite_k, "failed:", e$message))
        return(NULL)
        # message_list <- list(message = paste("Iteration", ite_k, "failed:", e$message)) 
        # return(message_list)
      })
      
    }
  
  gc()
  #---------------------------------------------------------
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
  
  colnames(save_df_ite) <- col_names
  colnames(save_df_itep2) <- col_names
  #save_df <- rbind(save_df,save_df_ite)
  
  #+++++++++++++++ 2
  
  
  color_setting <- c(
    # 前8色（他人方法）冷调中性色组
    "#5F7A8AB3",  # 蓝灰（色相205°, 明度50%, 饱和度30%）
    "#A9A9A9B3",  # 深灰（中性基准色）
    "#5B6B95B3",  # 蓝紫灰（色相235°, 明度52%, 透明度70%）
    "#BC8F8FB3",  # 玫瑰棕（低纯度粉灰）
    "#7FAC8FB3",  # 绿灰（色相150°, 明度58%）
    "#C2A87FB3",  # 卡其灰（明度60%, 降饱和度处理）
    "#9CB8D6B3",  # 淡钢蓝（饱和度25%）
    "#B07D4FB3",   # 暖灰棕（色相35°, 明度55%）
    
    # 后3色（自研方法）高对比强调色
    "#FF6B6B",    # 珊瑚红（色相0°, 饱和度93%, 明度85%）
    "#4ECDC4",    # 青蓝（色相175°, 饱和度80%, 明度78%）
    #"#FFC843"     # 琥珀黄（色相45°, 饱和度90%, 明度85%）
    "#E3D15CB3"
  )
  
  # color_setting <- c( "#708090B3",   "#A9A9A9B3", "#6B6799B3", "#BC8F8FB3",  "#8FBC8FB3",  "#D2B48CB3",  "#B0C4DEB3","#CD853FB3","#FF6B6B", "#4ECDC4", "#FFD93D")
  # color_setting <- c("#708090","#A9A9A9","#6B6799","#BC8F8F","#8FBC8F","#D2B48C","#B0C4DE","#CD853F","#FF6B6B","#4ECDC4","#FFD93D")
  # color_setting <- c('#3A8870', '#A66D45', '#6B6799', '#B23D6D', '#5C8F3A', '#CC9B32', '#7D5E2E', '#808080', '#FF6B4A', '#6E8FCC', '#5CBD7D')
  # color_setting <- c("#1B9E77", "#D95F02" ,"#7570B3", "#E7298A" ,"#66A61E" ,"#E6AB02" ,"#A6761D" ,"#666666","#66C2A5", "#FC8D62", "#8DA0CB")
  zihao1 <- 15 #x轴标题字号
  zihao2 <- 12 #x轴刻度字号
  zihao3 <- 15 #图例文本字号
  kedu_up <- ceiling(quantile(save_df_ite$beta_hat,0.995)*10)/10
  kedu_down <- floor(quantile(save_df_ite$beta_hat,0.005)*10)/10
  kedu_fen <- 0.1
  breaks_set <- round(seq(kedu_down,kedu_up,kedu_fen),1)
  
  hline1 <- beta[1]
  #hline2 <- beta[2]
  
  P1<-ggplot(data=save_df_ite,aes( y = beta_hat)) +
    geom_boxplot(aes(fill = method),outlier.alpha=0.5) +
    scale_fill_manual(values=color_setting)+
    geom_hline(aes(yintercept = hline1), linetype = "dashed", color = "darkred") +
    theme_bw() +
    scale_y_continuous(limits = c(kedu_down, kedu_up), 
                       breaks = breaks_set)+
    #xlab(expression(beta[1])) +
    ylab("Estimation") +
    theme(axis.text.x = element_text(angle = 0,vjust = 0.5,hjust = 0.5))+ 
    theme(strip.text.x = element_text(size = 35),
          strip.text.y = element_text(size = 32)
    ) +
    theme(plot.title = element_text(size = 25,hjust = 0.5,vjust = 0.5),
          axis.text = element_text(size = zihao2),
          axis.title = element_text(size = zihao1),
          legend.text = element_text(size = zihao3),
          legend.title = element_text(size = zihao3))+  # 调整图例文本大小
    theme(axis.text.x = element_blank(),  # 隐藏 x 轴坐标标签
          axis.ticks.x = element_blank())+ # 隐藏 x 轴刻度线
    guides(fill=guide_legend(title="Method"))
  
  plot_list[[g]] <- P1
  stopImplicitCluster()
  closeAllConnections()
}


# 自动计算行列数
calc_layout <- function(n) {
  rows <- floor(sqrt(n))  # 计算行数
  cols <- ceiling(n / rows)  # 计算列数
  return(c(rows, cols))
}

layout <- calc_layout(length(plot_list))  # 获取布局
rows <- layout[1]
cols <- layout[2]

p_tests <- ggarrange(plotlist = plot_list, nrow = rows, ncol = cols,
                     common.legend = TRUE,
                     legend="right"
                     )
print(p_tests)

p_tests <- ggarrange(pyt1,pyt2,pyt3,pyt4,pyt5,
                     ncol=1,nrow=5,
                     #labels = c('(A)','(B)','(C)','(D)','(E)','(F)'),
                     common.legend = TRUE,
                     legend="right")


pyt_combine <- ggarrange(pyt1,pyt2,pyt3,pyt4,pyt5,
                         ncol=1,nrow=5,
                         #labels = c('(A)','(B)','(C)','(D)','(E)','(F)'),
                         common.legend = TRUE,
                         legend="right")

ggsave(paste('beta_',date_mark,'_',ver,'.png',sep = ''),plot = P1, bg="white", path = path_output_plots,width = 12,height = 9,dpi = 600)
ggsave(paste('beta_',date_mark,'_',ver,'.pdf',sep = ''),plot = P1, bg="white", path = path_output_plots,width = 12,height = 9,dpi = 600)

#===============================================================================绘图1
# time_end <- Sys.time()
# time_consume <- difftime(time_end,time_start)
# time_consume 
stopImplicitCluster()
closeAllConnections()
