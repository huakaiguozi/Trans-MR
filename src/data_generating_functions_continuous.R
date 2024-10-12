data_generate_X <- function(n,p_Z,p_Z_max,p_C,param_X){
  #set.seed(seed)
  #---------param_X
  #h <- param_X$heritability
  alpha <- param_X$alpha
  #---
  pi <- param_X$pi
  pi2 <- param_X$pi2
  #---
  eta_1 <- param_X$eta_1
  U <- rnorm(n,0,1)
  epsi_x <- rnorm(n,0,1)
  #---------
  if(p_C == 1){
    p_C_list <- c(1,0,0)
  }else if(p_C == 2){
    p_C_list <- c(1,1,0)
  }else if(p_C >= 3){
    p_C_list <- c(1,1,(p_C-2))
  }
  #---------
  snp_pr_vector <- seq(0.12,0.4,(0.4-0.12)/(p_Z_max-1))
  Z <- matrix(rbinom(n = n*p_Z_max, size = 2, prob = rep(snp_pr_vector,each=n)), nrow = n)
  # #---------分类+连续协变量产生
  # prob_1 <- logistic_function(0) 
  # C_1 <- matrix(rbinom(n = p_C_list[1] * n, size = 1, prob = prob_1), nrow = n)
  # C_2 <- abs(matrix(rnorm(n = p_C_list[2] * n, mean = 2, sd = 2), nrow = n))
  # C_3 <- abs(matrix(rnorm(n = p_C_list[3] * n, mean = 1, sd = 1), nrow = n))
  #---------连续协变量产生
  C_1 <- abs(matrix(rnorm(n = p_C_list[1] * n, mean = 2, sd = 1), nrow = n))
  C_2 <- abs(matrix(rnorm(n = p_C_list[2] * n, mean = 3, sd = 2), nrow = n))
  C_3 <- abs(matrix(rnorm(n = p_C_list[3] * n, mean = 1.5, sd = 1), nrow = n))
  #---------------分类协变量产生
  # prob_1 <- logistic_function(0) 
  # prob_2 <- logistic_function(0.5) 
  # prob_3 <- logistic_function(1) 
  # C_1 <- matrix(rbinom(n = p_C_list[1] * n, size = 1, prob = prob_1), nrow = n)
  # C_2 <- matrix(rbinom(n = p_C_list[2] * n, size = 2, prob = prob_2), nrow = n)
  # C_3 <- matrix(rbinom(n = p_C_list[3] * n, size = 3, prob = prob_3), nrow = n)
  #---------------
  C <- cbind(C_1,C_2,C_3)
  #---------
  diag_result1 <- rowSums((Z %*% pi) * C)
  #diag_result2 <- rowSums((Z %*% pi2) * exp(C))
  # X <- Z %*% alpha + C %*% eta_1 + diag(Z %*% pi %*% t(C)) + diag(Z %*% pi2 %*% t(exp(C))) + U + epsi_x 
  X <- Z %*% alpha + C %*% eta_1 + U + epsi_x + diag_result1 #+ diag_result2
  # X <- Z %*% alpha + C %*% eta_1 + U + epsi_x 
  Z_return <- Z[,1:p_Z]
  #---------遗传度计算
  var_Z <- sum(apply(Z_return, 2, var))
  Var_X <- var(X)
  h_compute <- var_Z / Var_X 
  #---------summray数据
  res <- summary_compute(Z=Z_return,X=X)
  beta <- res$beta
  sd <- res$sd
  p <- res$p
  eaf <- res$eaf
  R2 <- res$R2
  F_stats <- res$F_stats
  #---------
  param_true <- list(alpha=alpha,pi=pi)
  #---------
  data_list <- list(Z=Z_return,C=C,X=X,beta=beta,sd=sd,p=p,eaf=eaf,param=param_true,h_c=h_compute,R2=R2,F_stats=F_stats)
  return(data_list)
}
################################################################################
data_generate_Y <- function(n,p_Z,p_Z_max,p_C,param_X,param_Y){
  #set.seed(seed)
  #---------param_X
  #h <- param_X$heritability
  alpha <- param_X$alpha
  pi <- param_X$pi
  pi2 <- param_X$pi2
  eta_1 <- param_X$eta_1
  #---------param_Y
  beta <- param_Y$beta
  eta_2 <- param_Y$eta_2
  gamma_random <- rbinom(p_Z_max,1,param_Y$ratio_pleiotropy)
  gamma <- param_Y$gamma * gamma_random
  heter_level <- heter_level
  #---------error
  U <- rnorm(n,0,1)
  epsi_x <- rnorm(n,0,1)
  epsi_y <- rnorm(n,0,1)
  #---------
  if(p_C == 1){
    p_C_list <- c(1,0,0)
  }else if(p_C == 2){
    p_C_list <- c(1,1,0)
  }else if(p_C >= 3){
    p_C_list <- c(1,1,(p_C-2))
  }
  #---------
  snp_pr_vector <- seq(0.15,0.35,(0.35-0.15)/(p_Z_max-1))
  Z <- matrix(rbinom(n = n*p_Z_max, size = 2, prob = rep(snp_pr_vector,each=n)), nrow = n)
  #---------
  
  if(heter_level != 0){
    # #---------分类+连续协变量产生
    # prob_1 <- logistic_function(0 + (heter_level)) 
    # C_1 <- matrix(rbinom(n = p_C_list[1] * n, size = 1, prob = prob_1), nrow = n)
    # C_2 <- abs((1+heter_level) * matrix(rnorm(n = p_C_list[2] * n, mean = 2, sd = 2), nrow = n))
    # C_3 <- abs((1+heter_level) * matrix(rnorm(n = p_C_list[3] * n, mean = 1, sd = 1), nrow = n))
    #---------连续协变量产生
    C_1 <- ( sqrt(exp(heter_level)) * abs(matrix(rnorm(n = p_C_list[1] * n, mean = 2, sd = 1), nrow = n)) )
    C_2 <- ( sqrt(exp(heter_level)) * abs(matrix(rnorm(n = p_C_list[2] * n, mean = 3, sd = 2), nrow = n)) )
    C_3 <- ( sqrt(exp(heter_level)) * abs(matrix(rnorm(n = p_C_list[3] * n, mean = 1.5, sd = 1), nrow = n)) )
    # #---------------分类协变量产生
    # prob_1 <- logistic_function(0 - (heter_level)) 
    # prob_2 <- logistic_function(0.5 - (heter_level)) 
    # prob_3 <- logistic_function(1 + (heter_level)) 
    # C_1 <- matrix(rbinom(n = p_C_list[1] * n, size = 1, prob = prob_1), nrow = n)
    # C_2 <- matrix(rbinom(n = p_C_list[2] * n, size = 2, prob = prob_2), nrow = n)
    # C_3 <- matrix(rbinom(n = p_C_list[3] * n, size = 3, prob = prob_3), nrow = n)
  }else{
    # #---------分类+连续协变量产生
    # prob_1 <- logistic_function(0) 
    # C_1 <- matrix(rbinom(n = p_C_list[1] * n, size = 1, prob = prob_1), nrow = n)
    # C_2 <- abs(matrix(rnorm(n = p_C_list[2] * n, mean = 2, sd = 2), nrow = n))
    # C_3 <- abs(matrix(rnorm(n = p_C_list[3] * n, mean = 1, sd = 1), nrow = n))
    #---------连续协变量产生
    C_1 <- abs(matrix(rnorm(n = p_C_list[1] * n, mean = 2, sd = 2), nrow = n))
    C_2 <- abs(matrix(rnorm(n = p_C_list[2] * n, mean = 3, sd = 1.5), nrow = n))
    C_3 <- abs(matrix(rnorm(n = p_C_list[3] * n, mean = 1.5, sd = 1), nrow = n))
    # #----------------分类协变量产生
    # prob_1 <- logistic_function(0) 
    # prob_2 <- logistic_function(0.5) 
    # prob_3 <- logistic_function(1) 
    # C_1 <- matrix(rbinom(n = p_C_list[1] * n, size = 1, prob = prob_1), nrow = n)
    # C_2 <- matrix(rbinom(n = p_C_list[2] * n, size = 2, prob = prob_2), nrow = n)
    # C_3 <- matrix(rbinom(n = p_C_list[3] * n, size = 3, prob = prob_3), nrow = n)
  }
  C <- cbind(C_1,C_2,C_3)
  #---------
  diag_result1 <- rowSums((Z %*% pi) * C)
  #diag_result2 <- rowSums((Z %*% pi2) * exp(C))
  #X <- Z %*% alpha + C %*% eta_1 + diag(Z %*% pi %*% t(C)) + diag(Z %*% pi2 %*% t(exp(C)))  + U + epsi_x 
  #X <- Z %*% alpha + C %*% eta_1 + diag(Z %*% pi %*% t(C)) + U + epsi_x 
  X <- Z %*% alpha + C %*% eta_1+ U + epsi_x + diag_result1 #+ diag_result2
  Y <- Z %*% gamma + X %*% beta + C %*% eta_2 + U  + epsi_y
  
  Z_return <- Z[,1:p_Z]
  #---------遗传度计算
  var_Z <- sum(apply(Z_return, 2, var))
  Var_X <- var(X)
  h_compute <- var_Z / Var_X 
  #---------summray数据
  res <- summary_compute(Z=Z_return,X=Y)
  beta <- res$beta
  sd <- res$sd
  p <- res$p
  eaf <- res$eaf
  R2 <- res$R2
  F_stats <- res$F_stats
  #---------
  param_true <- list(alpha=alpha,pi=pi)
  #---------
  data_list <- list(Z=Z_return,C=C,X=X, Y=Y,beta=beta,sd=sd,p=p,eaf=eaf,param=param_true,h_c=h_compute,R2=R2,F_stats=F_stats)
  return(data_list)
}

#-------------------------------------------------------------------------------
Total_F <- function(n, k, R2) {
  return(((n - k - 1) / k) * (sum(R2) / (1 - sum(R2))))
}

#---summary数据
summary_compute <- function(Z,X){
  N <- dim(Z)[1]
  p_Z <- dim(Z)[2]
  #------------------
  beta_summary <- rep(NA,p_Z)
  beta_summary_p <- rep(NA,p_Z)
  
  beta_summary_sd <- rep(NA,p_Z)
  R_2 <- rep(NA,p_Z)
  F_stats <- rep(NA,p_Z) 
  
  beta_summary_eaf <- colSums(Z)/ (2 * N)  #也是正确的
  #------------------
  for (z in 1:p_Z){
    data_ivw_x <- as.data.frame(cbind(X,Z[,z]))
    names(data_ivw_x) <- c('x','z')
    Beta_res <- feols(x ~ z, data=data_ivw_x)
    Beta_res_coef <- coef(Beta_res)  ##报错1
    if ('(Intercept)' %in% names(Beta_res_coef)){
      beta_hat <- Beta_res_coef[2]
    }else{
      beta_hat <- Beta_res_coef[1]
    }
    #------------------summary数据
    beta_summary[z] <- beta_hat 
    beta_summary_sd[z] <- Beta_res$coeftable[,2][2]  
    beta_summary_p[z] <- Beta_res$coeftable$`Pr(>|t|)`[2]
    # f_1 <- sum(Z[,z]==1)/N  #也可以
    # f_2 <- sum(Z[,z]==2)/N
    # beta_summary_eaf[z] <- f_1/2 + f_2
  }
  
  for(i in 1:p_Z){
    var_z <- 2 * beta_summary_eaf[i] *  (1-beta_summary_eaf[i])
    term1 <- beta_summary[i]^2 * var_z
    R_2[i] <- term1 / ( term1 + (beta_summary_sd[i]^2 * N * var_z ) )
  }
  
  for(j in 1:p_Z){
    F_stats[j] <- ( ( N-2 ) * R_2[j] ) / ( 1-R_2[j]  )
  }
  
  res_list <- list(beta = beta_summary, sd=beta_summary_sd, p=beta_summary_p,
                   eaf=beta_summary_eaf, R2=R_2, F_stats=F_stats )
  return(res_list)
}


logistic_function <- function(x) {
  return(1 / (1 + exp(-x)))
}


dimention_C <- function(p){
  list_1 <- rnorm(mean = p, sd = p/2 ,n = 2)
  list_2 <- round(list_1/sum(list_1) * p*(2/3))
  list_3 <- c(list_2,(p-sum(list_2)))
  return(list_3)
}
