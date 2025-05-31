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
                        , .packages = c('fixest','TwoSampleMR','glmnet','mr.raps','penalized','ridge','FNN','data.table','dplyr'))%dopar%
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
          if(p_Zm >= 1){
            for (i in 1:p_Zm) {  
              model2 <- glm(Y_0~X_hat_matrix_OR[,i] + C_0 + Z_0m[,i], family = binomial(link = "log"))  #有Z
              summary_model2 <- summary(model2)
              pval_model2 <- summary_model2$coefficients[,4]
              Z_pval <- pval_model2[length(pval_model2)]
              if(Z_pval < pval_i){
                b_Y_list2[i] <- model2$coefficients[2]
                sd_Y_list2[i] <- summary_model2$coefficients[,'Std. Error'][2]
              }else{
              model2_1 <- glm(Y_0~X_hat_matrix_OR[,i] + C_0, family = binomial(link = "log"))  #无Z
              b_Y_list2[i] <- model2_1$coefficients[2]
              sd_Y_list2[i] <- summary(model2_1)$coefficients[,'Std. Error'][2]
              }
            }
          }
          b_Y_list2 <- c(b_Y_list2, beta_hat_degrade)
          sd_Y_list2 <- c(sd_Y_list2, se_degrade)
          
          sd_Y_list2_inv <- 1/sd_Y_list2
          beta_hat_2 <- sum(sd_Y_list2_inv^2 * b_Y_list2)/sum(sd_Y_list2_inv^2)
          # sd_Y_list2_inv <- 1/(exp(b_Y_list2) * sd_Y_list2) 
          # beta_hat_2 <- log(sum(sd_Y_list2_inv^2 * exp(b_Y_list2))/sum(sd_Y_list2_inv^2))
          #---------------------------------------------------------------Trans_IPW0          31
          if(p_Zm>=1){
            
            for (i in 1:p_Zm) {
              
              model3 <- glm(Y_0~X_hat_matrix_IPW[,i] + C_0 + Z_0m[,i], family = binomial(link = "log"))  #有Z
              summary_model3 <- summary(model3)
              pval_model3 <- summary_model3$coefficients[,4]
              Z_pval <- pval_model3[length(pval_model3)]
              print(as.numeric(Z_pval))
              if(Z_pval < pval_i){
                b_Y_list3[i] <- model3$coefficients[2]
                sd_Y_list3[i] <- summary(model3)$coefficients[,'Std. Error'][2]
              }else{
              model3_1 <- glm(Y_0~X_hat_matrix_IPW[,i] + C_0, family = binomial(link = "log"))  #无Z
              b_Y_list3[i] <- model3_1$coefficients[2]
              sd_Y_list3[i] <- summary(model3_1)$coefficients[,'Std. Error'][2]
              }
            }
          }
          
          b_Y_list3 <- c(b_Y_list3, beta_hat_degrade)
          sd_Y_list3 <- c(sd_Y_list3, se_degrade)
          
          sd_Y_list3_inv <- 1/sd_Y_list3
          beta_hat_3 <- sum(sd_Y_list3_inv^2 * b_Y_list3)/sum(sd_Y_list3_inv^2)
          # sd_Y_list3_inv <- 1/(exp(b_Y_list3) * sd_Y_list3) 
          # beta_hat_3 <- log(sum(sd_Y_list3_inv^2 * exp(b_Y_list3))/sum(sd_Y_list3_inv^2))
          #---------------------------------------------------------------Trans_AIPW0          41
          if(p_Zm>=1){
            
            for (i in 1:p_Zm) {
              model4 <- glm(Y_0~X_hat_matrix_AIPW[,i] + C_0 + Z_0m[,i], family = binomial(link = "log"))  #有Z
              summary_model4 <- summary(model4)
              pval_model4 <- summary_model4$coefficients[,4]
              Z_pval <- pval_model4[length(pval_model4)]
              if(Z_pval<pval_i){
                b_Y_list4[i] <- model4$coefficients[2]
                sd_Y_list4[i] <- summary(model4)$coefficients[,'Std. Error'][2]
              }else{
              model4_1 <- glm(Y_0~X_hat_matrix_AIPW[,i] + C_0, family = binomial(link = "log"))  #无Z
              b_Y_list4[i] <- model4_1$coefficients[2]
              sd_Y_list4[i] <- summary(model4_1)$coefficients[,'Std. Error'][2]
              }
            }
          }
          b_Y_list4 <- c(b_Y_list4, beta_hat_degrade)
          sd_Y_list4 <- c(sd_Y_list4, se_degrade)
          
          sd_Y_list4_inv <- 1/sd_Y_list4
          beta_hat_4 <- sum(sd_Y_list4_inv^2 * b_Y_list4)/sum(sd_Y_list4_inv^2)
          # sd_Y_list4_inv <- 1/(exp(b_Y_list4) * sd_Y_list4) 
          # beta_hat_4 <- log(sum(sd_Y_list4_inv^2 * exp(b_Y_list4))/sum(sd_Y_list4_inv^2))
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
        for (i in 1:p_Z) {  
          model2 <- glm(Y_0~X_hat_matrix_OR[,i] + C_0 + Z_0[,i], family = binomial(link = "log"))  #有Z
          summary_model2 <- summary(model2)
          pval_model2 <- summary_model2$coefficients[,4]
          Z_pval <- pval_model2[length(pval_model2)]
          print(as.numeric(Z_pval))
          #if(Z_pval < pval_i){
          if(Z_pval < pval_i){
            b_Y_list2[i] <- model2$coefficients[2]
            sd_Y_list2[i] <- summary(model2)$coefficients[,'Std. Error'][2]
            print(model2$coefficients[2])
          }else{
          model2_1 <- glm(Y_0~X_hat_matrix_OR[,i] + C_0, family = binomial(link = "log"))  #无Z
          b_Y_list2[i] <- model2_1$coefficients[2]
          sd_Y_list2[i] <- summary(model2_1)$coefficients[,'Std. Error'][2]
          print(model2_1$coefficients[2])
          }
        }
        # cor(X_hat_matrix_OR[,1],Z_0[,1])
        sd_Y_list2_inv <- 1/sd_Y_list2
        beta_hat_2 <- sum(sd_Y_list2_inv^2 * b_Y_list2)/sum(sd_Y_list2_inv^2)
        # sd_Y_list2_inv <- 1/(exp(b_Y_list2) * sd_Y_list2) 
        # beta_hat_2 <- log(sum(sd_Y_list2_inv^2 * exp(b_Y_list2))/sum(sd_Y_list2_inv^2))
        #---------------------------------------------------------------Trans_IPW0          31
        for (i in 1:p_Z) {
          
          # model_r3_1 <- lm(X_hat_matrix_IPW[,i]~Z_0[,i])
          # X_hat3_r <- model_r3_1$residuals #X_hat_matrix_IPW[,i] - summary(model_r3_1)$coefficients [2,1] * Z_0[,i]
          # model_r3_2 <- glm(Y_0~X_hat3_r + C_0 + Z_0[,i], family = binomial(link = "log"))
          # #----------------------
          # O_x <- cbind(X_hat_matrix_IPW[,i], C_0, Z_0[,i])
          # cv_model <- cv.glmnet(
          #   x = O_x,
          #   y = Y_0,
          #   family = binomial(link = "log"),
          #   alpha = 1,  # Lasso 正则化
          #   #standardize = TRUE  # 建议标准化以提升稳定性
          #   nfolds = 5
          # )
          # 
          # best_lambda <- cv_model$lambda.min
          # 
          # best_model <- glmnet(
          #   x = O_x,
          #   y = Y_0,
          #   family = binomial(link = "log"),
          #   alpha = 1,  # Lasso 正则化
          #   lambda = best_lambda,
          #   standardize = TRUE  # 建议标准化以提升稳定性
          #   #nfolds = 5
          # )
          
          #----------------------
          model3 <- glm(Y_0~X_hat_matrix_IPW[,i] + C_0 + Z_0[,i], family = binomial(link = "log"))  #有Z
          summary_model3 <- summary(model3)
          pval_model3 <- summary_model3$coefficients[,4]
          Z_pval <- pval_model3[length(pval_model3)]
          print(as.numeric(Z_pval))
          if(Z_pval < pval_i){
            b_Y_list3[i] <- model3$coefficients[2]
            sd_Y_list3[i] <- summary(model3)$coefficients[,'Std. Error'][2]
            print(model3$coefficients[2])
          }else{
          model3_1 <- glm(Y_0~X_hat_matrix_IPW[,i] + C_0, family = binomial(link = "log"))  #无Z
          b_Y_list3[i] <- model3_1$coefficients[2]
          sd_Y_list3[i] <- summary(model3_1)$coefficients[,'Std. Error'][2]
          print(model3_1$coefficients[2])
          }
          
        }
        sd_Y_list3_inv <- 1/sd_Y_list3
        beta_hat_3 <- sum(sd_Y_list3_inv^2 * b_Y_list3)/sum(sd_Y_list3_inv^2)
        # sd_Y_list3_inv <- 1/(exp(b_Y_list3) * sd_Y_list3) 
        # beta_hat_3 <- log(sum(sd_Y_list3_inv^2 * exp(b_Y_list3))/sum(sd_Y_list3_inv^2))
        #---------------------------------------------------------------Trans_AIPW0          41
        for (i in 1:p_Z) {
          model4 <- glm(Y_0~X_hat_matrix_AIPW[,i] + C_0 + Z_0[,i], family = binomial(link = "log"))  #有Z
          summary_model4 <- summary(model4)
          pval_model4 <- summary_model4$coefficients[,4]
          Z_pval <- pval_model4[length(pval_model4)]
          if(Z_pval < pval_i){
            b_Y_list4[i] <- model4$coefficients[2]
            sd_Y_list4[i] <- summary(model4)$coefficients[,'Std. Error'][2]
          }else{
          model4_1 <- glm(Y_0~X_hat_matrix_AIPW[,i] + C_0, family = binomial(link = "log"))  #无Z
          b_Y_list4[i] <- model4_1$coefficients[2]
          sd_Y_list4[i] <- summary(model4_1)$coefficients[,'Std. Error'][2]
          }
        }
        sd_Y_list4_inv <- 1/sd_Y_list4
        beta_hat_4 <- sum(sd_Y_list4_inv^2 * b_Y_list4)/sum(sd_Y_list4_inv^2)
        # sd_Y_list4_inv <- 1/(exp(b_Y_list4) * sd_Y_list4) 
        # beta_hat_4 <- log(sum(sd_Y_list4_inv^2 * exp(b_Y_list4))/sum(sd_Y_list4_inv^2))
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

bias_vec <- as.vector(apply(ite_df2-beta, 2, function(x) mean(x[abs(x) < 10], na.rm = TRUE)))
sd_vec <- as.vector(apply(ite_df2, 2, function(x) sd(x[abs(x) < 10], na.rm = TRUE)))
mse_vec <- bias_vec^2 + sd_vec^2
save_df_ite2 <- cbind(param_df2,c(bias_vec,sd_vec,mse_vec))
#save_df2 <- rbind(save_df2,save_df_ite2)

colnames(save_df_ite2) <- col_names2
#+++++++++++++++ 3

# 手算power
test_df <- save_df_ite
colnames(test_df) <- col_names
indice_power_list <- which(n_0 == n_0_list2)
#save_power_list[[indice_power_list]] <- test_df
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
#save_power_list2[[indice_power_list]] <- test_df2


#---------------------------------------------------------
# hist(X_0, breaks = 20, probability = TRUE, main = "Histogram with Normal Curve")
# curve(dnorm(x, mean = mean(x), sd = sd(x)), col = "red", lwd = 2, add = TRUE)

# 构建空列表用于存储各个 omega 值
beta1_list <- list()
beta2_list <- list()
beta3_list <- list()

# 遍历 beta_hat_ite，提取 omega1, omega2, omega3
for (i in 1:length(beta_hat_ite)) {
  beta1_list[[i]] <- as.numeric(beta_hat_ite[[i]]$beta1)
  beta2_list[[i]] <- as.numeric(beta_hat_ite[[i]]$beta2)
  beta3_list[[i]] <- as.numeric(beta_hat_ite[[i]]$beta3)
}

# 将列表拼接成矩阵，每行为一次取值
beta1_matrix <- do.call(rbind, beta1_list)
beta2_matrix <- do.call(rbind, beta2_list)
beta3_matrix <- do.call(rbind, beta3_list)

apply(sweep(beta2_matrix, 2, beta, FUN = "-"),2,mean)
#gamma 0.5   -0.04468899 -0.14231651
#      0.25  -0.02296020 -0.07177885
summary(beta1_matrix)
summary(beta2_matrix)
summary(beta3_matrix)


bias_1 <- apply(sweep(beta1_matrix, 2, beta, FUN = "-"),2,mean)
bias_2 <- apply(sweep(beta2_matrix, 2, beta, FUN = "-"),2,mean)
bias_3 <- apply(sweep(beta3_matrix, 2, beta, FUN = "-"),2,mean)

sd_1 <- apply(beta1_matrix, 2, sd)
sd_2 <- apply(beta2_matrix, 2, sd)
sd_3 <- apply(beta3_matrix, 2, sd)

mse_1 <- bias_1^2 + sd_1^2
mse_2 <- bias_2^2 + sd_2^2
mse_3 <- bias_3^2 + sd_3^2



dt2 <- c()
name_vec1 <- rep(paste0('beta'),3)
name_vec2 <- c('bias','sd','mse')
num_vec1 <- c(bias_1[1],sd_1[1],mse_1[1])
num_vec2 <- c(bias_2[1],sd_2[1],mse_2[1])
num_vec3 <- c(bias_3[1],sd_3[1],mse_3[1])
dt2 <- rbind(name_vec1,name_vec2,num_vec1,num_vec2,num_vec3)

rnames <- c('param','index','TLMR-OR','TLMR-IPW','TLMR-AIPW')
dt2 <- cbind(rnames,dt2)
fwrite(dt2,file=paste(path_output_tables, '/mse',date_mark,'_',ver,'.csv',sep = ''))
#-----------------------------------------------------------------------

#methods_ordered <- factor(c('TLMR-OR','TLMR-IPW','TLMR-AIPW'),levels = c('TLMR-OR','TLMR-IPW','TLMR-AIPW'))

# beta.1 <- c(beta1_matrix[,1],beta2_matrix[,1],beta3_matrix[,1])
# methods_vec <- c(rep(methods_ordered[1],dim(beta1_matrix)[1]) 
#                  ,rep(methods_ordered[2],dim(beta2_matrix)[1])
#                  ,rep(methods_ordered[3],dim(beta3_matrix)[1]))
# 
# #beta.2 <- c(beta1_matrix[,2],beta2_matrix[,2],beta3_matrix[,2])
# 
# # methods_vec <- c(rep(methods_ordered[1],dim(beta1_matrix)[1]) 
# #                  ,rep(methods_ordered[2],dim(beta2_matrix)[1])
# #                  ,rep(methods_ordered[3],dim(beta3_matrix)[1]))
# 
# #beta_names <- factor(c('',''),levels = c(''))
# 
# 
# data_beta.1 <- data.frame(beta=beta.1,method=methods_vec) 
# #data_beta.2 <- data.frame(beta=beta.2,method=methods_vec) 

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

P1

ggsave(paste('beta_',date_mark,'_',ver,'.png',sep = ''),plot = P1, bg="white", path = path_output_plots,width = 12,height = 9,dpi = 600)
ggsave(paste('beta_',date_mark,'_',ver,'.pdf',sep = ''),plot = P1, bg="white", path = path_output_plots,width = 12,height = 9,dpi = 600)

#===============================================================================绘图1
# time_end <- Sys.time()
# time_consume <- difftime(time_end,time_start)
# time_consume 
stopImplicitCluster()
closeAllConnections()
