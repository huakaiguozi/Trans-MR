data_result <- matrix(NA,nrow = 0,ncol = 7)
colnames(data_result) <- c('method','sample','outcome','beta','p','d1','d2')
data_result <- as.data.frame(data_result)




# 创建一个空向量用于存储结果
d1 <- rep(NA, nrow(data_result))
# 遍历每种method
for (m in unique(data_result$method)) {
  # 选出sample为0和1的两行
  subset_rows <- which(data_result$method == m)
  rows_sample_0 <- subset_rows[data_result$sample[subset_rows] == 0]
  rows_sample_1 <- subset_rows[data_result$sample[subset_rows] == 1]
  
  # 如果两个样本都存在
  if (length(rows_sample_0) > 0 & length(rows_sample_1) > 0) {
    # 比较beta的符号是否一致
    sign_0 <- sign(data_result$beta[rows_sample_0])
    sign_1 <- sign(data_result$beta[rows_sample_1])
    
    # 对每对sample 0 和 1 的 beta 比较符号是否一致（按 outcome 匹配）
    for (i in seq_along(rows_sample_0)) {
      outcome0 <- data_result$outcome[rows_sample_0[i]]
      j <- which(data_result$outcome[rows_sample_1] == outcome0)
      if (length(j) == 1) {
        d1[rows_sample_0[i]] <- as.integer(sign_0[i] == sign_1[j])
        d1[rows_sample_1[j]] <- as.integer(sign_0[i] == sign_1[j])
      }
    }
  }
}
data_result$d1 <- d1

# 初始化 d2 列为 NA
data_result$d2 <- NA

# 循环每种 method 和 sample
for (m in unique(data_result$method)) {
  for (s in unique(data_result$sample)) {
    # 提取对应行
    beta_fev1 <- data_result %>% filter(method == m, sample == s, outcome == "FEV1") %>% pull(beta)
    beta_fvc  <- data_result %>% filter(method == m, sample == s, outcome == "FVC")  %>% pull(beta)
    beta_ratio <- data_result %>% filter(method == m, sample == s, outcome == "FEV1/FVC") %>% pull(beta)
    
    # 若三者都存在（即长度为1），才计算
    if (length(beta_fev1) == 1 && length(beta_fvc) == 1 && length(beta_ratio) == 1) {
      k1 <- beta_fev1
      k2 <- beta_fvc
      k3 <- (1 + k1) / (1 + k2)
      
      # 符号一致性判断
      sign_k3 <- sign(k3 - 1)
      sign_ratio <- sign(beta_ratio)
      consistency <- ifelse(sign_k3 == sign_ratio, 1, 0)
      
      # 写入 d2 列中对应 FEVI/FVC 行
      data_result$d2[data_result$method == m & data_result$sample == s & data_result$outcome == "FEV1/FVC"] <- consistency
    }
  }
}

data_result$p_T <- as.integer(data_result$p < 0.05)

p_star_vec <- rep('-',nrow(data_result))
ind_star1 <- which(data_result$p < 0.05)
ind_star2 <- which(data_result$p < 0.005)
ind_star3 <- which(data_result$p < 0.001)
p_star_vec[ind_star1] <- '*'
p_star_vec[ind_star2] <- '**'
p_star_vec[ind_star3] <- '***'

data_result$p_star <- p_star_vec
#---------------------------------------------------------------------------------------------------
# data_result$T_all <- 



data_result_save_6 <- data_result
data_result_save_8 <- data_result

data_result_save_6[data_result_save_6 $method==method_names[1],]
data_result_save_6[data_result_save_6 $method==method_names[2],]
data_result_save_6[data_result_save_6 $method==method_names[3],]
data_result_save_6[data_result_save_6 $method==method_names[4],]
data_result_save_6[data_result_save_6 $method==method_names[5],]
data_result_save_6[data_result_save_6 $method==method_names[6],]
data_result_save_6[data_result_save_6 $method==method_names[7],]
data_result_save_6[data_result_save_6 $method==method_names[8],]
data_result_save_6[data_result_save_6 $method==method_names[9],]

data_result_save_8[data_result_save_8 $method==method_names[1],]
data_result_save_8[data_result_save_8 $method==method_names[2],]
data_result_save_8[data_result_save_8 $method==method_names[3],]
data_result_save_8[data_result_save_8 $method==method_names[4],]
data_result_save_8[data_result_save_8 $method==method_names[5],]
data_result_save_8[data_result_save_8 $method==method_names[6],]
data_result_save_8[data_result_save_8 $method==method_names[7],]
data_result_save_8[data_result_save_8 $method==method_names[8],]
data_result_save_8[data_result_save_8 $method==method_names[9],]


