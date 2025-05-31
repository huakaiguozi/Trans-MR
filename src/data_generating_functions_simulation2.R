# 
# alpha <- 1
# pi <- c(1,1)
# pi2 <- c(1,1)
# eta_1 <- c(1,1)
# eta_2 <- c(1,1)
# param_X <- list(
#   alpha= alpha,
#   pi = pi,
#   pi2 = pi2,
#   eta_1 = eta_1
#   
# )


data_generate_X <- function(n,p_Z,p_Z_max,p_C,param_X){
  #set.seed(seed)
  #U_md <- param_X$U_md
  #---------param_X
  #h <- param_X$heritability
  alpha <- param_X$alpha
  #---
  pi <- param_X$pi
  pi2 <- param_X$pi2
  #---
  eta_1 <- param_X$eta_1
  #U <- rnorm(n,0,(6/U_md))
  U <- rnorm(n,0,1)
  #epsi_x <- rnorm(n,0,1)
  #---------
  if(p_C == 1){
    p_C_list <- c(1,0,0)
  }else if(p_C == 2){
    p_C_list <- c(1,1,0)
  }else if(p_C >= 3){
    p_C_list <- c(1,1,(p_C-2))
  }
  #---------
  if(p_Z_max>=2){
    snp_pr_vector <- seq(0.2,0.4,(0.4-0.2)/(p_Z_max-1)) 
  }else{
    snp_pr_vector <- c(0.3)
  }
  
  Z <- matrix(rbinom(n = n*p_Z_max, size = 2, prob = rep(snp_pr_vector,each=n)), nrow = n)
  # #---------分类+连续协变量产生
  # prob_1 <- logistic_function(0) 
  # C_1 <- matrix(rbinom(n = p_C_list[1] * n, size = 1, prob = prob_1), nrow = n)
  # C_2 <- abs(matrix(rnorm(n = p_C_list[2] * n, mean = 2, sd = 2), nrow = n))
  # C_3 <- abs(matrix(rnorm(n = p_C_list[3] * n, mean = 1, sd = 1), nrow = n))
  # #---------分类+连续协变量产生
  #---------------分类协变量产生
  prob_1 <- logistic_function(0) 
  prob_2 <- logistic_function(0.5) 
  prob_3 <- logistic_function(1) 
  C_1 <- matrix(rbinom(n = p_C_list[1] * n, size = 1, prob = prob_1), nrow = n)
  C_2 <- matrix(rbinom(n = p_C_list[2] * n, size = 2, prob = prob_2), nrow = n)
  C_3 <- matrix(rbinom(n = p_C_list[3] * n, size = 3, prob = prob_3), nrow = n)
  #---------------分类协变量产生
  C <- cbind(C_1,C_2,C_3)
  #---------
  # Z <- Z ^ (floor(p_Z/3)/p_Z)
  # C <- C ^ ((p_Z - floor(p_Z/3))/p_Z)
  # Z <- Z * (1 / p_Z)
  # C <- C * (1 / p_Z)
  Z <- Z * (1 / sqrt(p_Z))
  C <- C * (1 / sqrt(p_Z))

  #---------
  diag_result1 <- rowSums((Z %*% pi) * C)
  #diag_result2 <- rowSums((Z %*% pi2) * exp(C))
  # X <- Z %*% alpha + C %*% eta_1 + diag(Z %*% pi %*% t(C)) + diag(Z %*% pi2 %*% t(exp(C))) + U + epsi_x 
  X <- Z %*% alpha + C %*% eta_1 + U + diag_result1 #+ diag_result2
  # X <- Z %*% alpha + C %*% eta_1 + U + epsi_x 
  Z_return <- Z[,1:p_Z]
  
  if(is.null(dim(Z_return))){
    Z_return <- matrix(Z_return,ncol = p_Z) 
  }
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
data_generate_Y_new2 <- function(n,p_Z,p_Z_max,p_C,param_X,param_Y){
  #set.seed(seed)
  #U_md <- param_X$U_md
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
  heter_level <- param_Y$heter_level
  #---------
  n_origin <- n
  #---------error
  U <- rnorm(n,0,0.1)
  epsi_x <- rnorm(n,0,0.1)
  #epsi_y <- rnorm(n,0,0.1)
  #---------
  if(p_C == 1){
    p_C_list <- c(1,0,0)
  }else if(p_C == 2){
    p_C_list <- c(1,1,0)
  }else if(p_C >= 3){
    p_C_list <- c(1,1,(p_C-2))
  }
  #---------
  if(p_Z_max>=2){
    snp_pr_vector <- seq(0.25,0.45,(0.45-0.25)/(p_Z_max-1)) 
  }else{
    snp_pr_vector <- c(0.35)
  }
  Z <- matrix(rbinom(n = n*p_Z_max, size = 2, prob = rep(snp_pr_vector,each=n)), nrow = n)
  #---------
  
  if(heter_level != 0){
    # #---------分类+连续协变量产生
    # prob_1 <- logistic_function(0 + (heter_level)) 
    # C_1 <- matrix(rbinom(n = p_C_list[1] * n, size = 1, prob = prob_1), nrow = n)
    # C_2 <- abs((1+heter_level) * matrix(rnorm(n = p_C_list[2] * n, mean = 2, sd = 2), nrow = n))
    # C_3 <- abs((1+heter_level) * matrix(rnorm(n = p_C_list[3] * n, mean = 1, sd = 1), nrow = n))
    # #---------分类+连续协变量产生
    #---------------分类协变量产生
    prob_1 <- logistic_function(0 - (heter_level)) 
    prob_2 <- logistic_function(0.5 - (heter_level)) 
    prob_3 <- logistic_function(1 + (heter_level)) 
    C_1 <- matrix(rbinom(n = p_C_list[1] * n, size = 1, prob = prob_1), nrow = n)
    C_2 <- matrix(rbinom(n = p_C_list[2] * n, size = 2, prob = prob_2), nrow = n)
    C_3 <- matrix(rbinom(n = p_C_list[3] * n, size = 3, prob = prob_3), nrow = n)
    #---------------分类协变量产生
  }else{
    # #---------分类+连续协变量产生
    # prob_1 <- logistic_function(0) 
    # C_1 <- matrix(rbinom(n = p_C_list[1] * n, size = 1, prob = prob_1), nrow = n)
    # C_2 <- abs(matrix(rnorm(n = p_C_list[2] * n, mean = 2, sd = 2), nrow = n))
    # C_3 <- abs(matrix(rnorm(n = p_C_list[3] * n, mean = 1, sd = 1), nrow = n))
    # #---------分类+连续协变量产生
    #----------------分类协变量产生
    prob_1 <- logistic_function(0) 
    prob_2 <- logistic_function(0.5) 
    prob_3 <- logistic_function(1) 
    C_1 <- matrix(rbinom(n = p_C_list[1] * n, size = 1, prob = prob_1), nrow = n)
    C_2 <- matrix(rbinom(n = p_C_list[2] * n, size = 2, prob = prob_2), nrow = n)
    C_3 <- matrix(rbinom(n = p_C_list[3] * n, size = 3, prob = prob_3), nrow = n)
    #---------------分类协变量产生
  }
  C <- cbind(C_1,C_2,C_3)
  #------------------
  Z <- Z * (1 / sqrt(p_Z))
  C <- C * (1 / sqrt(p_Z))
  # Z <- Z * (1 / p_Z)
  # C <- C * (1 / p_Z)
  #------------------
  diag_result1 <- rowSums((Z %*% pi) * C)
  
  X <- Z %*% alpha + C %*% eta_1+ U  + diag_result1 + epsi_x
  
  #Intercept <- -ceiling(max(Z %*% gamma + X %*% beta + C %*% eta_2 + U + epsi_y)*10)/10
  Intercept <- -ceiling(max(Z %*% gamma + X %*% beta + C %*% eta_2 + U)*10)/10
  #Intercept <- param_Y$Intercept
  
  Y_linear <- Z %*% gamma + X %*% beta + C %*% eta_2 + U +Intercept #+ epsi_y
  P_Y <- exp(Y_linear)
  
  Y <- rbinom(n,1,P_Y)
  
  # Y <- rep(NA,n)
  # 
  # Y[which(P_Y >= 1)] <- 1 
  # 
  # Y[which(P_Y < 1)] <- rbinom(sum(P_Y<1),1,P_Y[which(P_Y < 1)])
  
  Z_return <- Z[,1:p_Z]
  if(is.null(dim(Z_return))){
    Z_return <- matrix(Z_return,ncol = p_Z) 
  }
  #---------遗传度计算
  var_Z <- sum(apply(Z_return, 2, var))
  Var_X <- var(X)
  h_compute <- var_Z / Var_X 
  #---------summray数据
  res <- summary_compute(Z=Z_return,X=Y)
  beta_s <- res$beta
  sd <- res$sd
  p <- res$p
  eaf <- res$eaf
  R2 <- res$R2
  F_stats <- res$F_stats
  #---------
  param_true <- list(alpha=alpha,pi=pi)
  #---------
  data_list <- list(Z=Z_return,C=C,X=X, Y=Y,beta=beta_s,sd=sd,p=p,eaf=eaf,param=param_true,h_c=h_compute,R2=R2,F_stats=F_stats)
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
  indicator_bin <- if(n_distinct(X) < 5){TRUE}else{FALSE}  # TRUE:分类；FALSE:连续
  #------------------
  if(!indicator_bin){
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
  }else{
    beta_summary <- rep(NA,p_Z)
    beta_summary_p <- rep(NA,p_Z)
    
    beta_summary_sd <- rep(NA,p_Z)
    R_2 <- rep(NA,p_Z)
    F_stats <- rep(NA,p_Z) 
    
    beta_summary_eaf <- colSums(Z)/ (2 * N)  #也是正确的
    #------------------
    for (z in 1:p_Z){
      # data_ivw_x <- as.data.frame(cbind(X,Z[,z]))
      # names(data_ivw_x) <- c('x','z')
      # Beta_res <- feols(x ~ z, data=data_ivw_x)
      
      model_glm <- glm(X~Z[,z], family = binomial(link = "log"))
      summary_glm <- summary(model_glm)
      Beta_res_coef <- summary_glm$coefficients[,1]
      Beta_res_sd <- summary_glm$coefficients[,2]
      Beta_res_p <- summary_glm$coefficients[,4]
      #Beta_res_coef <- coef(Beta_res)  ##报错1
      if ('(Intercept)' %in% names(Beta_res_coef)){
        beta_hat <- as.numeric(Beta_res_coef[2])
        sd_hat <- as.numeric(Beta_res_sd[2])
        p_hat <- as.numeric(Beta_res_p[2])
      }else{
        beta_hat <- as.numeric(Beta_res_coef[1])
        sd_hat <- as.numeric(Beta_res_sd[1])
        p_hat <- as.numeric(Beta_res_p[1])
      }
      #------------------summary数据
      beta_summary[z] <- beta_hat 
      beta_summary_sd[z] <- sd_hat 
      beta_summary_p[z] <- p_hat
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
