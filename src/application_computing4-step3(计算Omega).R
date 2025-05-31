#Trans_OR
X_hat_matrix_OR <- Trans_OLS2(Z_0=Z_0,Z_1=Z_1,C_0=C_0,C_1=C_1,X_1=X_1)$X_hat_matrix
#Trans_IPW
X_hat_matrix_IPW <- IPW2(Z_0=Z_0,Z_1=Z_1,C_0=C_0,C_1=C_1,X_1=X_1)$X_hat_matrix
#Trans_AIPW
X_hat_matrix_AIPW <- AIPW2(Z_0=Z_0,Z_1=Z_1,C_0=C_0,C_1=C_1,X_1=X_1
)$X_hat_matrix

omega_list_1 <- matrix(nrow = p_Zm,ncol = p_C)
omega_list_2 <- matrix(nrow = p_Zm,ncol = p_C)
omega_list_3 <- matrix(nrow = p_Zm,ncol = p_C)

cov_tensor_1 <- array(dim=c(p_C,p_C,p_Zm))
cov_tensor_2 <- array(dim=c(p_C,p_C,p_Zm))
cov_tensor_3 <- array(dim=c(p_C,p_C,p_Zm))

for (q in 1:p_Zm) {
  
  model2_1 <- lm(Y_0~X_hat_matrix_OR[,q] + C_0_origin + Z_0m[,q]+ X_hat_matrix_OR[,q] * C_0_origin )  #有Z
  # model2_2 <- lm(Y_0~X_hat_matrix_IPW[,q] + C_0+ Z_0m[,q]+ X_hat_matrix_IPW[,q] * C_0)
  # model2_3 <- lm(Y_0~X_hat_matrix_AIPW[,q] + C_0+ Z_0m[,q]+ X_hat_matrix_AIPW[,q] * C_0)
  # 
  omega_1 <- model2_1$coefficients[c((length(model2_1$coefficients)-p_C+1):length(model2_1$coefficients))]
  # omega_2 <- model2_2$coefficients[c((length(model2_1$coefficients)-p_C+1):length(model2_1$coefficients))]
  # omega_3 <- model2_3$coefficients[c((length(model2_1$coefficients)-p_C+1):length(model2_1$coefficients))]
  # 
  omega_list_1[q,] <- omega_1
  # omega_list_2[q,] <- omega_2
  # omega_list_3[q,] <- omega_3
  
  cov_tensor_1[,,q] <- vcov(model2_1)[c((length(model2_1$coefficients)-p_C+1):length(model2_1$coefficients)),c((length(model2_1$coefficients)-p_C+1):length(model2_1$coefficients))]
  # cov_tensor_2[,,q] <- vcov(model2_2)[c((length(model2_1$coefficients)-p_C+1):length(model2_1$coefficients)),c((length(model2_1$coefficients)-p_C+1):length(model2_1$coefficients))]
  # cov_tensor_3[,,q] <- vcov(model2_3)[c((length(model2_1$coefficients)-p_C+1):length(model2_1$coefficients)),c((length(model2_1$coefficients)-p_C+1):length(model2_1$coefficients))]
  
}

cov_tensor_1_inv <- array(dim=c(p_C,p_C,p_Zm))
for (w in 1:dim(omega_list_1)[1]) {
  
  cov1_inv <- solve(cov_tensor_1[,,w])
  cov_tensor_1_inv[,,w] <- cov1_inv
  
}

cov_inv_sum_1_inv <-solve( apply(cov_tensor_1_inv, c(1, 2), sum)  )

weight_tensor_1 <- array(dim=c(p_C,p_C,p_Zm))
for (e in 1:dim(omega_list_1)[1]) {
  
  weight_tensor_1[,,e] <- cov_inv_sum_1_inv %*% cov_tensor_1_inv[,,e]
  
}

omega_1_slices <- matrix(nrow = p_Zm,ncol = p_C)
for (r in 1:dim(omega_list_1)[1]) {
 
  omega_1_slice <- weight_tensor_1[,,r] %*% omega_list_1[r,]
  omega_1_slices[r,] <- omega_1_slice
  
}

omega_1_ivw <- apply(omega_1_slices,2,sum) 

Wald_statistic_1 <- t(omega_1_ivw) %*% solve(cov_inv_sum_1_inv) %*% omega_1_ivw

p_value_1 <- 1 - pchisq(Wald_statistic_1, df = 3)


#-----------------
var_1 <- array(dim = 9)
omega_11 <- omega_list_1[,3]

for (t in 1:9) {
  var_1[t] <- cov_tensor_1[3,3,t]
}

theta_hat <- omega_11
se_theta <- sqrt(var_1)

# 计算权重矩阵（逐元素计算权重）
w <- 1 / (se_theta^2)  # 计算每个参数的权重

# 计算 IVW 估计值（逐维计算）
theta_ivw <- sum(w * theta_hat) / sum(w)  

# 计算标准误（逐维计算）
se_ivw <- sqrt(1 / sum(w))  

# 计算 Z 统计量
z_score <- theta_ivw / se_ivw  

# 计算 P 值（双侧检验）
p_value <- 2 * (1 - pnorm(abs(z_score))) 


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
