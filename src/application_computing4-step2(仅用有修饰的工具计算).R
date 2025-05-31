Z_0 <- Z_0[,indice_modify]
Z_1 <- Z_1[,indice_modify]

p_Z <- ncol(Z_1)

Z_1 <- apply(Z_1, 2, function(x) as.numeric(as.character(x)))

res_1 <- summary_compute(Z=Z_1,X=X_1,C=AS_1)
beta_1 <- res_1$beta
sd_1 <- res_1$sd
p_1 <- res_1$p
eaf_1 <- res_1$eaf
R2_1 <- res_1$R2
F_stats_1 <- res_1$F_stats

data_1 <- list(Z=Z_1,C=C_1,X=X_1,beta=beta_1,sd=sd_1,p=p_1,eaf=eaf_1,R2=R2_1,F_stats=F_stats_1)

res_0 <- summary_compute(Z=Z_0,X=Y_0,C=AS_0)
beta_0 <- res_0$beta
sd_0 <- res_0$sd
p_0 <- res_0$p
eaf_0 <- res_0$eaf
R2_0 <- res_0$R2
F_stats_0 <- res_0$F_stats

Z_0 <- apply(Z_0, 2, function(x) as.numeric(as.character(x)))

data_0 <- list(Z=Z_1,C=C_1,X=X_1, Y=Y_0,beta=beta_0,sd=sd_0,p=p_0,eaf=eaf_0,R2=R2_0,F_stats=F_stats_0)

beta_hat_1_list <- Twosamle_package(p_Z=p_Z,p_C=p_C,data_0=data_0,data_1=data_1,beta=1)

beta_hat_1 <- beta_hat_1_list$beta

b_Y_list2 <- c()
sd_Y_list2 <- c()
b_Y_list3 <- c()
sd_Y_list3 <- c()
b_Y_list4 <- c()
sd_Y_list4 <- c()

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
  model2 <- lm(Y_0~X_hat_matrix_OR[,i] + C_0 + Z_0[,i])  #有Z
  summary_model2 <- summary(model2)
  pval_model2 <- summary_model2$coefficients[,4]
  Z_pval <- pval_model2[length(pval_model2)]
  print(as.numeric(Z_pval))
  #if(Z_pval < pval_i){
  if(Z_pval < pval_i){
    b_Y_list2[i] <- model2$coefficients[2]
    sd_Y_list2[i] <- summary(model2)$coefficients[,'Std. Error'][2]
  }else{
    model2_1 <- lm(Y_0~X_hat_matrix_OR[,i] + C_0)  #无Z
    b_Y_list2[i] <- model2_1$coefficients[2]
    sd_Y_list2[i] <- summary(model2_1)$coefficients[,'Std. Error'][2]
  }
  
}
# cor(X_hat_matrix_OR[,1],Z_0[,1])
sd_Y_list2_inv <- 1/sd_Y_list2 
beta_hat_2 <- sum(sd_Y_list2_inv^2 * b_Y_list2)/sum(sd_Y_list2_inv^2)
#---------------------------------------------------------------Trans_IPW0          31
for (i in 1:p_Z) {
  
  model3 <- lm(Y_0~X_hat_matrix_IPW[,i] + C_0 + Z_0[,i])  #有Z
  summary_model3 <- summary(model3)
  pval_model3 <- summary_model3$coefficients[,4]
  Z_pval <- pval_model3[length(pval_model3)]
  print(as.numeric(Z_pval))
  if(Z_pval < pval_i){
    b_Y_list3[i] <- model3$coefficients[2]
    sd_Y_list3[i] <- summary(model3)$coefficients[,'Std. Error'][2]
  }else{
    model3_1 <- lm(Y_0~X_hat_matrix_IPW[,i] + C_0)  #无Z
    b_Y_list3[i] <- model3_1$coefficients[2]
    sd_Y_list3[i] <- summary(model3_1)$coefficients[,'Std. Error'][2]
  }
  
  
}
sd_Y_list3_inv <- 1/sd_Y_list3 
beta_hat_3 <- sum(sd_Y_list3_inv^2 * b_Y_list3)/sum(sd_Y_list3_inv^2)
#---------------------------------------------------------------Trans_AIPW0          41
for (i in 1:p_Z) {
  model4 <- lm(Y_0~X_hat_matrix_AIPW[,i] + C_0 + Z_0[,i])  #有Z
  summary_model4 <- summary(model4)
  pval_model4 <- summary_model4$coefficients[,4]
  Z_pval <- pval_model4[length(pval_model4)]
  if(Z_pval < pval_i){
    b_Y_list4[i] <- model4$coefficients[2]
    sd_Y_list4[i] <- summary(model4)$coefficients[,'Std. Error'][2]
  }else{
    model4_1 <- lm(Y_0~X_hat_matrix_AIPW[,i] + C_0)  #无Z
    b_Y_list4[i] <- model4_1$coefficients[2]
    sd_Y_list4[i] <- summary(model4_1)$coefficients[,'Std. Error'][2]
  }
}
sd_Y_list4_inv <- 1/sd_Y_list4
beta_hat_4 <- sum(sd_Y_list4_inv^2 * b_Y_list4)/sum(sd_Y_list4_inv^2)

#beta_hat_vec <- c(beta_hat_1, beta_hat_2, beta_hat_3,beta_hat_4)
beta_hat_vec <- c(beta_hat_1, beta_hat_2)

