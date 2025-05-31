residual_1 <- sum(Y_0-X_0 * (-0.024))/dim(Y_0)[1]
residual_1


residual_2 <- sum(Y_0-X_0 * (-0.005))/dim(Y_0)[1]
residual_2



effects_1 <- c(-0.01,-0.005,-0.01,-0.012,-0.015,-0.009,-0.014,-0.014,-0.024)


effects_1 <- beta_hat_vec
effects_1

residuals_1 <- c()
cor_1 <- c()
R2_1 <- c()

for(i in 1:length(effects_1)){
  effect <- effects_1[i]
  Y_0_hat <- X_0 * effect
  
  residual_i <- sum(abs(Y_0-Y_0_hat))/dim(Y_0)[1]
  residuals_1 <- c(residuals_1,residual_i)
  
  cor_i <- cor(Y_0, Y_0_hat)
  # print(head(X_0 * effect))
  cor_1 <- c(cor_1, cor_i)
  
  RSS <- sum((Y_0 - Y_0_hat)^2)
  TSS <- sum((Y_0 - mean(Y_0))^2)
  R2_i <- 1- (RSS/TSS)
  R2_1 <- c(R2_1,R2_i)
  
}


residuals_1
cor_1
R2_1
