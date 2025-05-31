#-------------------------------------------------------------------------------
Total_F <- function(n, k, R2) {
  return(((n - k - 1) / k) * (sum(R2) / (1 - sum(R2))))
}

#---summary数据
summary_compute <- function(Z,X,C){
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
    data_ivw_x <- as.data.frame(cbind(X,Z[,z],C[,1],C[,2]))
    names(data_ivw_x) <- c('x','z','c1','c2')
    Beta_res <- feols(x ~ z+c1+c2, data=data_ivw_x)
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

