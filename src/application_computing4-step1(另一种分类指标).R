#cluster::
library(foreach)
library(doParallel)


n_divides <- 5
indicator_print<-1

path_data_read2 <- 'data/50W全数据所有结局数据_250516-4.csv'
data_k00 <- fread(path_data_read2)

loc_1 <- data_k00$loc_1
remove(data_k00)

path_data_read <- 'data/50W全数据所有结局数据_250512.csv'
data_k0 <- fread(path_data_read)
#---------------------------------改人群划分
# G1_codes <- c(
#   10003,  # Stockport        -  3,794
#   11009,  # Newcastle        - 36,995
#   11010,  # Leeds            - 44,184
#   11016,  # Liverpool        - 32,799
#   #  11002,  # Oxford           - 14,054*
#   11003,  # Cardiff          - 17,875*
#   # 11007,  # Reading          - 29,399*
#   11020,  # Croydon          - 27,363
#   11011,  # Bristol          - 43,000*
#   11013  # Nottingham       - 33,870*
# )
# 
# G2_codes <- c(
#   11001,  # Manchester       - 13,937*
#   11004,  # Glasgow          - 18,644*
#   11005,  # Edinburgh        - 17,193*
#   11006,  # Stoke            - 19,426*
#   11008,  # Bury             - 28,311*
#   11017,   # Middlesborough   - 21,283*
#   11012,  # Barts            - 12,574
#   11014,  # Sheffield        - 30,381
#   11018,  # Hounslow         - 28,866
#   # 11020,  # Croydon          - 27,363
#   11007,  # Reading          - 29,399*
#   11021,  # Birmingham       - 25,493
#   11022,  # Swansea          - 2,280
#   11002,  # Oxford           - 14,054*
#   11023   # Wrexham          -    649      271804
# )
# 
# # 0：1：
# loc_2 <- ifelse(loc_1 %in% G1_codes, 0,
#                 ifelse(loc_1 %in% G2_codes, 1, NA))
# 
# data_k0$地理划分 <- loc_2

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%数据预处理
#--------
data_outcomes <- data_k0[,c(103:112)]

data_outcomes_continuous <- data_outcomes[,c(4:7)]
data_outcomes_binary <- data_outcomes[,c(1:3,8:10)]
remove(data_outcomes)

data_k1 <- data_k0[,c(1:102)]


#--------处理缺失
#"FVC"      "FEV1"     "FEV1/FVC" "PEF"  
data_outcome <- data_outcomes_continuous[,3]

data_k2 <- cbind(data_k1,data_outcome)
indice_row1 <- rep(FALSE,nrow(data_k2))
for (i in 1:ncol(data_k2)) {
  data_col <- data_k2[,..i]
  indice_row1 <- indice_row1 | is.na(data_col) | (data_col=='') 
}

data_k3 <- data_k2[as.vector(!indice_row1),]
data_k4 <- data_k3[data_k3$种族背景 %in% c(1,1001,1002,1003),]

#------------------
data_all <- data_k4
rm(data_k2, data_k3, data_k4)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# 使用 grep 查找所有 rs 开头的列
rs_columns <- grep("^rs", colnames(data_all), value = TRUE)

Z_all <- data_all[, ..rs_columns]

#-----------------
Z_all_save <- Z_all
#-----------------

# # 提取这些列的数据
# Z_1 <- data_all[, ..rs_columns]
# Z_1 <- as.matrix(Z_1)

data_all2 <- data_all %>% select(-all_of(rs_columns))

#C_1 <- data_all2[,c(2:ncol(data_all2)),with = FALSE]

C_all <- data_all2[,c(2:ncol(data_all2)),with = FALSE]

# C_all <- C_all %>% select(-all_of(c("种族背景",'BMI','腰围','臀围','FEV1','地理划分')))

# C_all <- C_all %>% select(-all_of(c("种族背景",'BMI','中度运动天数/周','PEF','地理划分')))
# C_all <- C_all %>% select(-all_of(c("种族背景",'BMI','中度运动天数/周','FVC','地理划分')))
# C_all <- C_all %>% select(-all_of(c("种族背景",'BMI','中度运动天数/周','FEV1','地理划分')))
C_all <- C_all %>% select(-all_of(c("种族背景",'BMI','中度运动天数/周','FEV1/FVC','地理划分')))

AS_all <- C_all[,c('性别', '年龄')]

X_all <- data_all2$BMI
# Y_all <- data_all2$`PEF`
# Y_all <- data_all2$`FVC`
# Y_all <- data_all2$`FEV1`
Y_all <- data_all2$`FEV1/FVC`

#C_1 <- C_1 %>% select(-all_of(c("种族背景",'BMI','第1秒用力呼气容积(FEV1)')))


#-----------------------------------
# Step 1: 计算 WHR
#df$WHR <- df$waist / df$hip
WHR <- data_all2$腰围 / data_all2$臀围

# Step 2: 定义 sigmoid 函数
sigmoid <- function(x) {
  1 / (1 + exp(-x))
}

# Step 3: 设置中点和陡峭程度
WHR_m <- median(WHR, na.rm = TRUE)  # 中心点
WHR_a <- 20  # 陡峭程度，10~30之间效果不错，可调

# Step 4: 计算每个个体被分到 group 2 的概率
set.seed(12345)
prob_group2 <- sigmoid(WHR_a * (WHR - WHR_m))

# Step 5: 按概率随机分组
group <- rbinom(n = nrow(data_all2), size = 1, prob = prob_group2)
table(group)

indice_population_0 <- group == 0  #sample0
indice_population_1 <- group == 1  #sample1

indice_population_0 <- group == 1
indice_population_1 <- group == 0

#-----------------------------------


# basis_cov_name <- '中度运动天数/周'
# #basis_cov_name <- '年龄'
# #basis_cov_name <- 'BMI'
# 
# indice_basis_name <- which(colnames(data_all2) == basis_cov_name)
# 
# # basis_cov2 <- data_all2[[indice_basis_name]]
# 
# basis_cov <- data_all2[[indice_basis_name]]
# 
# index_low <- quantile(basis_cov,0.45)
# index_up <- quantile(basis_cov,0.55)
# 
# #选一个运行
# indice_population_0 <- basis_cov < index_low  #sample0
# indice_population_1 <- basis_cov > index_up   #sample1
# 
# indice_population_0 <- basis_cov > index_up
# indice_population_1 <- basis_cov < index_low

# basis_cov <- C_all[[indice_basis_name]]
# # indice_population_0 <- basis_cov >= summary(basis_cov)[3]
# # indice_population_1 <- basis_cov < summary(basis_cov)[3]
# 
# indice_population_0 <- basis_cov < summary(basis_cov)[3]
# indice_population_1 <- basis_cov >= summary(basis_cov)[3]

Z_all <- Z_all_save
#----------------------------------------第二次从这开始
Z_0 <- Z_all[indice_population_0,]
Z_1 <- Z_all[indice_population_1,]

C_0 <- C_all[indice_population_0,]
C_1 <- C_all[indice_population_1,]

AS_0 <- AS_all[indice_population_0,]  #计算summary数据时调整年龄性别
AS_1 <- AS_all[indice_population_1,]


# AS1 <- C_1[,c('性别','年龄')]  #计算summary数据时调整年龄性别
# AS0 <- C_0[,c('性别','年龄')]

X_0 <- X_all[indice_population_0]
X_1 <- X_all[indice_population_1]
X_0 <- matrix(X_0,nrow = nrow(Z_0))
X_1 <- matrix(X_1,nrow = nrow(Z_1))


Y_0 <- Y_all[indice_population_0]
Y_1 <- Y_all[indice_population_1]
Y_0 <- matrix(Y_0,nrow = nrow(Z_0))
Y_1 <- matrix(Y_1,nrow = nrow(Z_1))

# C_0 <- C_0 %>% select(-all_of(c(basis_cov_name)))
# C_1 <- C_1 %>% select(-all_of(c(basis_cov_name)))
#------------------------- 选好的效应修饰因子

# C_0 <- C_0 %>% select(all_of(c('腰围', 'C反应蛋白', '睡眠时长' )))
# C_1 <- C_1 %>% select(all_of(c('腰围', 'C反应蛋白', '睡眠时长' )))
#-------------------------
quantile_vec_all <- seq(0,1,1/n_divides)
quantile_vec <- quantile_vec_all[-c(1,length(quantile_vec_all))]
#------------------------------
C_1 <- as.data.frame(C_1)
C_1_cl <- C_1
for (i in c(1:dim(C_1)[2])) {
  #print(i)
  if( length(unique(C_1[,i])) > 10 ){
    #C_qt_vector <- quantile(C_1[,i],quantile_vec)
    #C_qt_label <- cut( C_0[,i],c(-Inf,C_qt_vector,Inf),labels=c(1:n_divides) )
    #C_0_cl[,i] <- C_qt_label
    # 计算分位数
    C_qt_vector <- unique(quantile(C_1[, i], quantile_vec))
    
    # 检查去重后是否有足够的分位数
    if (length(C_qt_vector) >= n_divides - 1) {
      # 使用去重后的分位数作为 cut 的 breaks 参数
      C_qt_label <- cut(C_1[, i], breaks = c(-Inf, C_qt_vector, Inf), 
                        labels = 1:n_divides, include.lowest = TRUE)
      C_1_cl[, i] <- C_qt_label
    } else {
      C_qt_label <- cut(C_1[, i], breaks = c(-Inf, C_qt_vector, Inf), 
                        labels = 1:(length(C_qt_vector)+1), include.lowest = TRUE)
      C_1_cl[, i] <- C_qt_label
      
      warning(paste("第", i, "列的分位数不足以进行划分,","现分为",(length(C_qt_vector)+1),"类"))
    }
    
    
    # C_qt_label <- cut( C_1[,i],c(-Inf,C_qt_vector,Inf),labels=c(1:n_divides),include.lowest=TRUE )
    # C_1_cl[,i] <- C_qt_label
  }
}
# 替换成分类变量；原始数据存储至origin
#C_0_origin <- C_0
C_1_origin <- C_1
#C_0 <- C_0_cl
C_1 <- C_1_cl
C_1 <- as.matrix(C_1)
C_1_cl <- as.matrix(C_1_cl)
#------------------------------
p_Z <- ncol(Z_1)
p_C <- ncol(C_1)

C_1 <- apply(C_1, 2, function(x) as.numeric(as.character(x)))
Z_1 <- apply(Z_1, 2, function(x) as.numeric(as.character(x)))
C_1_origin <- apply(C_1_origin, 2, function(x) as.numeric(as.character(x)))
#C_0_cl <- apply(C_0_cl, 2, function(x) as.numeric(as.character(x)))
################################################################################
#------------------------------
# X_1 <- data_all[["BMI"]]
# X_1 <- matrix(X_1,nrow = nrow(data_all))
#
res_1 <- summary_compute(Z=Z_1,X=X_1,C=AS_1)

beta_1 <- res_1$beta
sd_1 <- res_1$sd
p_1 <- res_1$p
eaf_1 <- res_1$eaf
R2_1 <- res_1$R2
F_stats_1 <- res_1$F_stats

#!!!!!!!!!!!!!!!!!!!!!!!从头开始,第二次忽略
sum(p_1 <= 5e-3)
sum(p_1 <= 5e-4)
sum(p_1 <= 5e-6)
sum(p_1 <= 5e-8)

indice_e4 <- which(p_1 <= 5e-4)
indice_e6 <- which(p_1 <= 5e-6)
indice_e8 <- which(p_1 <= 5e-8)


Z_all <- Z_all_save[,..indice_e6]

#!!!!!!!!!!!!!!!!!!!!!!!
Z_0 <- Z_all[indice_population_0,]
Z_1 <- Z_all[indice_population_1,]

C_0 <- C_all[indice_population_0,]
C_1 <- C_all[indice_population_1,]

AS_0 <- AS_all[indice_population_0,]  #计算summary数据时调整年龄性别
AS_1 <- AS_all[indice_population_1,]


# AS1 <- C_1[,c('性别','年龄')]  #计算summary数据时调整年龄性别
# AS0 <- C_0[,c('性别','年龄')]

X_0 <- X_all[indice_population_0]
X_1 <- X_all[indice_population_1]
X_0 <- matrix(X_0,nrow = nrow(Z_0))
X_1 <- matrix(X_1,nrow = nrow(Z_1))


Y_0 <- Y_all[indice_population_0]
Y_1 <- Y_all[indice_population_1]
Y_0 <- matrix(Y_0,nrow = nrow(Z_0))
Y_1 <- matrix(Y_1,nrow = nrow(Z_1))

# C_0 <- C_0 %>% select(-all_of(c(basis_cov_name)))
# C_1 <- C_1 %>% select(-all_of(c(basis_cov_name)))
#------------------------- 选好的效应修饰因子

# C_0 <- C_0 %>% select(all_of(c('腰围', 'C反应蛋白', '睡眠时长' )))
# C_1 <- C_1 %>% select(all_of(c('腰围', 'C反应蛋白', '睡眠时长' )))
#-------------------------
quantile_vec_all <- seq(0,1,1/n_divides)
quantile_vec <- quantile_vec_all[-c(1,length(quantile_vec_all))]
#------------------------------
C_1 <- as.data.frame(C_1)
C_1_cl <- C_1
for (i in c(1:dim(C_1)[2])) {
  #print(i)
  if( length(unique(C_1[,i])) > 10 ){
    #C_qt_vector <- quantile(C_1[,i],quantile_vec)
    #C_qt_label <- cut( C_0[,i],c(-Inf,C_qt_vector,Inf),labels=c(1:n_divides) )
    #C_0_cl[,i] <- C_qt_label
    # 计算分位数
    C_qt_vector <- unique(quantile(C_1[, i], quantile_vec))
    
    # 检查去重后是否有足够的分位数
    if (length(C_qt_vector) >= n_divides - 1) {
      # 使用去重后的分位数作为 cut 的 breaks 参数
      C_qt_label <- cut(C_1[, i], breaks = c(-Inf, C_qt_vector, Inf), 
                        labels = 1:n_divides, include.lowest = TRUE)
      C_1_cl[, i] <- C_qt_label
    } else {
      C_qt_label <- cut(C_1[, i], breaks = c(-Inf, C_qt_vector, Inf), 
                        labels = 1:(length(C_qt_vector)+1), include.lowest = TRUE)
      C_1_cl[, i] <- C_qt_label
      
      warning(paste("第", i, "列的分位数不足以进行划分,","现分为",(length(C_qt_vector)+1),"类"))
    }
    
    
    # C_qt_label <- cut( C_1[,i],c(-Inf,C_qt_vector,Inf),labels=c(1:n_divides),include.lowest=TRUE )
    # C_1_cl[,i] <- C_qt_label
  }
}
# 替换成分类变量；原始数据存储至origin
#C_0_origin <- C_0
C_1_origin <- C_1
#C_0 <- C_0_cl
C_1 <- C_1_cl
C_1 <- as.matrix(C_1)
C_1_cl <- as.matrix(C_1_cl)
#------------------------------
p_Z <- ncol(Z_1)
p_C <- ncol(C_1)

C_1 <- apply(C_1, 2, function(x) as.numeric(as.character(x)))
Z_1 <- apply(Z_1, 2, function(x) as.numeric(as.character(x)))
C_1_origin <- apply(C_1_origin, 2, function(x) as.numeric(as.character(x)))
#C_0_cl <- apply(C_0_cl, 2, function(x) as.numeric(as.character(x)))
################################################################################
#------------------------------
# X_1 <- data_all[["BMI"]]
# X_1 <- matrix(X_1,nrow = nrow(data_all))
#
res_1 <- summary_compute(Z=Z_1,X=X_1,C=AS_1)

beta_1 <- res_1$beta
sd_1 <- res_1$sd
p_1 <- res_1$p
eaf_1 <- res_1$eaf
R2_1 <- res_1$R2
F_stats_1 <- res_1$F_stats


data_1 <- list(Z=Z_1,C=C_1,X=X_1,beta=beta_1,sd=sd_1,p=p_1,eaf=eaf_1,R2=R2_1,F_stats=F_stats_1)
################################################################################
#------------------------------
# Y_0 <- data_all[["第1秒用力呼气容积(FEV1)"]]
# Y_0 <- matrix(Y_0,nrow = nrow(data_all))

#
res_0 <- summary_compute(Z=Z_0,X=Y_0,C=AS_0)
beta_0 <- res_0$beta
sd_0 <- res_0$sd
p_0 <- res_0$p
eaf_0 <- res_0$eaf
R2_0 <- res_0$R2
F_stats_0 <- res_0$F_stats
#------------------
C_0 <- as.data.frame(C_0)
C_0_cl <- C_0
for (i in c(1:dim(C_0)[2])) {
  if( length(unique(C_0[,i])) > 10 ){
    # 计算分位数
    C_qt_vector <- unique(quantile(C_0[, i], quantile_vec))
    
    # 检查去重后是否有足够的分位数
    if (length(C_qt_vector) >= n_divides - 1) {
      # 使用去重后的分位数作为 cut 的 breaks 参数
      C_qt_label <- cut(C_0[, i], breaks = c(-Inf, C_qt_vector, Inf), 
                        labels = 1:n_divides, include.lowest = TRUE)
      C_0_cl[, i] <- C_qt_label
    } else {
      C_qt_label <- cut(C_0[, i], breaks = c(-Inf, C_qt_vector, Inf), 
                        labels = 1:(length(C_qt_vector)+1), include.lowest = TRUE)
      C_0_cl[, i] <- C_qt_label
      
      warning(paste("第", i, "列的分位数不足以进行划分,","现分为",(length(C_qt_vector)+1),"类"))
    }
    
  }
}


C_0_origin <- C_0
#C_0 <- C_0_cl
C_0 <- C_0_cl
C_0 <- as.matrix(C_0)
C_0_cl <- as.matrix(C_0_cl)
#------------------------------
# p_Z <- ncol(Z_1)
# p_C <- ncol(C_1)

C_0 <- apply(C_0, 2, function(x) as.numeric(as.character(x)))
Z_0 <- apply(Z_0, 2, function(x) as.numeric(as.character(x)))
C_0_origin <- apply(C_0_origin, 2, function(x) as.numeric(as.character(x)))


#-------
# Z_0<-Z_1
# C_0<-C_1
# C_0_origin <-C_1_origin
#------------------
data_0 <- list(Z=Z_0,C=C_0,X=X_0, Y=Y_0,beta=beta_0,sd=sd_0,p=p_0,eaf=eaf_0,R2=R2_0,F_stats=F_stats_0)
################################################################################
beta_hat_1_list <- Twosamle_package(p_Z=p_Z,p_C=p_C,data_0=data_0,data_1=data_1,beta=1)

beta_hat_1 <- beta_hat_1_list$beta

beta_hat_1_list
# #-----------------------------------若不计算交互矩阵，就到这为止
#=================================================================================================
# 并行设置：检测并注册可用核心
n_cores <- parallel::detectCores() - 2  # 保留一个核心
cl <- makeCluster(n_cores)
registerDoParallel(cl)

# 创建空矩阵
modifier_matrix <- matrix(0, nrow = p_Z, ncol = p_C)
modifier_p <- matrix(NA, nrow = p_Z, ncol = p_C)

# 外层对 k 并行
results_list <- foreach(k = 1:p_Z, .packages = c("metafor")) %dopar% {
  modifier_row <- numeric(p_C)
  modifier_p_row <- numeric(p_C)
  
  for (j in 1:p_C) {
    data_chow_test <- data.frame(Z = Z_1[,k], X = X_1, C = C_1[,j])
    data_chow_test <- data_chow_test[order(data_chow_test$C), ]
    C_group_vec <- table(data_chow_test$C)
    
    end_indice <- cumsum(C_group_vec)
    start_indice <- c(1, head(end_indice, -1) + 1)
    
    beta_group_vec <- c()
    se_group_vec <- c()
    
    for (i in 1:length(C_group_vec)) {
      if (as.numeric(C_group_vec[i]) < 2) next
      data_chow_test_sub <- data_chow_test[start_indice[i]:end_indice[i], ]
      model_group_C <- lm(X ~ Z, data = data_chow_test_sub)
      res_model_group <- summary(model_group_C)$coefficients
      if (nrow(res_model_group) < 2) next
      beta_group_vec <- c(beta_group_vec, res_model_group[2,1])
      se_group_vec <- c(se_group_vec, res_model_group[2,2])
    }
    
    if (length(beta_group_vec) > 1) {
      res_meta_analysis <- metafor::rma.uni(yi = beta_group_vec, sei = se_group_vec, method = "REML")
      p_Q <- res_meta_analysis$QEp
      modifier_p_row[j] <- p_Q
      modifier_row[j] <- as.numeric(p_Q < 0.001)
    } else {
      modifier_p_row[j] <- NA
      modifier_row[j] <- 0
    }
  }
  
  list(modifier_row = modifier_row, modifier_p_row = modifier_p_row)
}

# 汇总结果
for (k in 1:p_Z) {
  modifier_matrix[k, ] <- results_list[[k]]$modifier_row
  modifier_p[k, ] <- results_list[[k]]$modifier_p_row
}

# 释放资源
stopCluster(cl)


C_0_save <- C_0
C_1_save <- C_1
C_0_origin_save <- C_0_origin
C_1_origin_save <- C_1_origin
#-----------------
indice_Cnouse <- which(colSums(modifier_matrix) <1)
#indice_Cnouse <- c(11,12,which(colSums(modifier_matrix) <1))
C_0 <- C_0_save
C_1 <- C_1_save
C_0_origin <- C_0_origin_save
C_1_origin <- C_1_origin_save

#indice_Cuse <- 
setdiff(c(1:p_C),indice_Cnouse)
C_1 <- C_1[,-indice_Cnouse]
C_1_origin <- C_1_origin[,-indice_Cnouse]
C_0 <- C_0[,-indice_Cnouse]
C_0_origin <- C_0_origin[,-indice_Cnouse]
#=================================================================================================
# #--------------------------------------画图1
# modifier_p <- as.data.table(modifier_p)
# 
# rownames(modifier_p) <- colnames(Z_1)
# colnames(modifier_p) <- colnames(C_1)
# 
# colnames(modifier_p) <- c("Gender", "Smoking", "Alcohol Frequency", "Glucose", 
#                           "Total Cholesterol", "HDL", "LDL", "Triglycerides", 
#                           "DBP", "SBP", "Waist Circumference", "Hip Circumference", 
#                           "Height", "Moderate Exercise Days/Week", "C-Reactive Protein", 
#                           "Sleep Duration", "Age", "Population Density", "Income Score", 
#                           "Cooked Vegetable Intake", "Raw Vegetable Intake", 
#                           "Natural Environment Percentage", "SNP")
# 
# 
# modifier_p$SNP <- colnames(Z_1)
# 
# modifier_p_long <- melt(modifier_p, id.vars = "SNP", variable.name = "Trait", value.name = "Value")
# 
# modifier_p_tf1 <- log10(modifier_p_long$Value)
# 
# modifier_p_tf2 <- (modifier_p_tf1- (-16)  )/(0-(-16))
# 
# modifier_p_long$Value2 <- modifier_p_tf2
# 
# modifier_p_long$SNP <- factor(modifier_p_long$SNP, levels = unique(modifier_p_long$SNP))
# modifier_p_long$Trait <- factor(modifier_p_long$Trait, levels = sort(unique(modifier_p_long$Trait)))
# 
# ###手动调色
# my_palette1 <- colorRampPalette(c('#DC1623','white'))(100000)
# my_palette2 <- colorRampPalette(c( 'white','white'))(100000)
# my_palette_combine <- c(my_palette1,my_palette2)
# #my_palette_combine <- c(my_palette1)
# value_breaks <- c(seq(0, 0.875, length.out = 100000), seq(0.875, 1, length.out = 100000))
# #value_breaks <- c(seq(0, 1, length.out = 100000))
# ###
# p1 <- ggplot(modifier_p_long, aes(Trait,SNP))+
#   #热图
#   geom_tile(aes(fill=Value2)
#             , color = "grey90"
#             
#   )+
#   #基于Importance(%)结果绘制圆圈
#   #geom_point(data=df_Importance2[df_Importance2$value>0,], aes(variable, OTUID, size = value), shape = 1)+
#   #主题设置
#   theme_bw()+
#   theme(panel.grid = element_blank(),
#         axis.text.y = element_text(size=10, color = "black"),
#         axis.text.x = element_text(size=10, color = "black", angle = 45, hjust = 1, vjust = 1))+
#   labs(x = NULL, y = NULL, fill = "Significance(P-value)"
#        #, size = "Importance(%)"
#        
#   )+
#   scale_fill_gradientn(limit = c(0,1)
#                        , colors = my_palette_combine
#                        , values = value_breaks
#                        , breaks = c(0,0.875,1)
#                        , labels = c('1e-16','1e-2','1')
#   )  #values在设置时应该在0-1之间
# #scale_fill_gradientn(limit = c(0,1), colors = my_palette_combine)
# 
# #+scale_size_continuous(range = c(2,7))
# p1
# ggsave(paste('交互_英文',date_mark,'_',ver,'.png',sep = ''),plot = p1, bg="white", path = path_output_plots,width = 20,height = 9,dpi = 600)
# ggsave(paste('交互_英文',date_mark,'_',ver,'.pdf',sep = ''),plot = p1, bg="white", path = path_output_plots,width = 20,height = 9,dpi = 600)
# 
# 
# #--------------
# 
# 
# #-----------
# 
# 
# #------
# indice_Cnouse <- which(colSums(modifier_matrix) <5)
# C_1 <- C_1[,-indice_Cnouse]
# C_1_origin <- C_1_origin[,-indice_Cnouse]
# C_0 <- C_0[,-indice_Cnouse]
# C_0_origin <- C_0_origin[,-indice_Cnouse]
#----------------------------如果不做热图，就从这里开始
b_Y_list2 <- c()
sd_Y_list2 <- c()
# b_Y_list3 <- c()
# sd_Y_list3 <- c()
# b_Y_list4 <- c()
# sd_Y_list4 <- c()

# modifier_matrix_save <- modifier_matrix
# modifier_matrix <- modifier_matrix[,setdiff(c(1:22),indice_Cnouse)]

indice_modify <- which(rowSums(modifier_matrix) > 0)
indice_notmodify <- setdiff(c(1:p_Z),which(rowSums(modifier_matrix) > 0))
p_Zm <- length(indice_modify)
p_Znm <- length(indice_notmodify)

pval_i <- 0.05
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
  
  beta_hat_res <- Twosamle_package4(p_Z=length(indice_notmodify),p_C=p_C,data_0=data3_0,data_1=data3_1)
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
    X_hat_matrix_OR <- Trans_OLS2(Z_0=Z_0m,Z_1=Z_1m,C_0=C_0_origin,C_1=C_1_origin,X_1=X_1)$X_hat_matrix
    #Trans_IPW
    # X_hat_matrix_IPW <- IPW2(Z_0=Z_0m,Z_1=Z_1m,C_0=C_0,C_1=C_1,X_1=X_1)$X_hat_matrix
    # #Trans_AIPW
    # X_hat_matrix_AIPW <- AIPW2(Z_0=Z_0m,Z_1=Z_1m,C_0=C_0,C_1=C_1,X_1=X_1
    # )$X_hat_matrix
    #-----------------------------------------
    
    #---------------------------------------------------------------Trans_OR0         21
    if(p_Zm >= 1){
      for (i in 1:p_Zm) {  
        model2 <- lm(Y_0~X_hat_matrix_OR[,i] + C_0_origin + Z_0m[,i])  #有Z
        summary_model2 <- summary(model2)
        pval_model2 <- summary_model2$coefficients[,4]
        Z_pval <- pval_model2[length(pval_model2)]
        if(Z_pval < pval_i){
          b_Y_list2[i] <- model2$coefficients[2]
          sd_Y_list2[i] <- summary(model2)$coefficients[,'Std. Error'][2]
        }else{
          model2_1 <- lm(Y_0~X_hat_matrix_OR[,i] + C_0_origin)  #无Z
          b_Y_list2[i] <- model2_1$coefficients[2]
          sd_Y_list2[i] <- summary(model2_1)$coefficients[,'Std. Error'][2]
        }
      }
    }
    b_Y_list2 <- c(b_Y_list2, beta_hat_degrade)
    sd_Y_list2 <- c(sd_Y_list2, se_degrade)
    
    sd_Y_list2_inv <- 1/sd_Y_list2 
    beta_hat_2 <- sum(sd_Y_list2_inv^2 * b_Y_list2)/sum(sd_Y_list2_inv^2)
    # #---------------------------------------------------------------Trans_IPW0          31
    # if(p_Zm>=1){
    #   
    #   for (i in 1:p_Zm) {
    #     
    #     model3 <- lm(Y_0~X_hat_matrix_IPW[,i] + C_0 + Z_0m[,i])  #有Z
    #     summary_model3 <- summary(model3)
    #     pval_model3 <- summary_model3$coefficients[,4]
    #     Z_pval <- pval_model3[length(pval_model3)]
    #     print(as.numeric(Z_pval))
    #     if(Z_pval < pval_i){
    #       b_Y_list3[i] <- model3$coefficients[2]
    #       sd_Y_list3[i] <- summary(model3)$coefficients[,'Std. Error'][2]
    #     }else{
    #       model3_1 <- lm(Y_0~X_hat_matrix_IPW[,i] + C_0)  #无Z
    #       b_Y_list3[i] <- model3_1$coefficients[2]
    #       sd_Y_list3[i] <- summary(model3_1)$coefficients[,'Std. Error'][2]
    #     }
    #   }
    # }
    # 
    # b_Y_list3 <- c(b_Y_list3, beta_hat_degrade)
    # sd_Y_list3 <- c(sd_Y_list3, se_degrade)
    # 
    # sd_Y_list3_inv <- 1/sd_Y_list3 
    # beta_hat_3 <- sum(sd_Y_list3_inv^2 * b_Y_list3)/sum(sd_Y_list3_inv^2)
    # #---------------------------------------------------------------Trans_AIPW0          41
    # if(p_Zm>=1){
    #   
    #   for (i in 1:p_Zm) {
    #     model4 <- lm(Y_0~X_hat_matrix_AIPW[,i] + C_0 + Z_0m[,i])  #有Z
    #     summary_model4 <- summary(model4)
    #     pval_model4 <- summary_model4$coefficients[,4]
    #     Z_pval <- pval_model4[length(pval_model4)]
    #     if(Z_pval<pval_i){
    #       b_Y_list4[i] <- model4$coefficients[2]
    #       sd_Y_list4[i] <- summary(model4)$coefficients[,'Std. Error'][2]
    #     }else{
    #       model4_1 <- lm(Y_0~X_hat_matrix_AIPW[,i] + C_0)  #无Z
    #       b_Y_list4[i] <- model4_1$coefficients[2]
    #       sd_Y_list4[i] <- summary(model4_1)$coefficients[,'Std. Error'][2]
    #     }
    #   }
    # }
    # b_Y_list4 <- c(b_Y_list4, beta_hat_degrade)
    # sd_Y_list4 <- c(sd_Y_list4, se_degrade)
    # 
    # sd_Y_list4_inv <- 1/sd_Y_list4
    # beta_hat_4 <- sum(sd_Y_list4_inv^2 * b_Y_list4)/sum(sd_Y_list4_inv^2)
  }
}else{
  #===============================================================迁移估计量
  #Trans_OR
  X_hat_matrix_OR <- Trans_OLS2(Z_0=Z_0,Z_1=Z_1,C_0=C_0,C_1=C_1,X_1=X_1)$X_hat_matrix
  # #Trans_IPW
  # X_hat_matrix_IPW <- IPW2(Z_0=Z_0,Z_1=Z_1,C_0=C_0,C_1=C_1,X_1=X_1)$X_hat_matrix
  # #Trans_AIPW
  # X_hat_matrix_AIPW <- AIPW2(Z_0=Z_0,Z_1=Z_1,C_0=C_0,C_1=C_1,X_1=X_1
  # )$X_hat_matrix
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
  # #---------------------------------------------------------------Trans_IPW0          31
  # for (i in 1:p_Z) {
  #   
  #   model3 <- lm(Y_0~X_hat_matrix_IPW[,i] + C_0 + Z_0[,i])  #有Z
  #   summary_model3 <- summary(model3)
  #   pval_model3 <- summary_model3$coefficients[,4]
  #   Z_pval <- pval_model3[length(pval_model3)]
  #   print(as.numeric(Z_pval))
  #   if(Z_pval < pval_i){
  #     b_Y_list3[i] <- model3$coefficients[2]
  #     sd_Y_list3[i] <- summary(model3)$coefficients[,'Std. Error'][2]
  #   }else{
  #     model3_1 <- lm(Y_0~X_hat_matrix_IPW[,i] + C_0)  #无Z
  #     b_Y_list3[i] <- model3_1$coefficients[2]
  #     sd_Y_list3[i] <- summary(model3_1)$coefficients[,'Std. Error'][2]
  #   }
  #   
  #   
  # }
  # sd_Y_list3_inv <- 1/sd_Y_list3 
  # beta_hat_3 <- sum(sd_Y_list3_inv^2 * b_Y_list3)/sum(sd_Y_list3_inv^2)
  # #---------------------------------------------------------------Trans_AIPW0          41
  # for (i in 1:p_Z) {
  #   model4 <- lm(Y_0~X_hat_matrix_AIPW[,i] + C_0 + Z_0[,i])  #有Z
  #   summary_model4 <- summary(model4)
  #   pval_model4 <- summary_model4$coefficients[,4]
  #   Z_pval <- pval_model4[length(pval_model4)]
  #   if(Z_pval < pval_i){
  #     b_Y_list4[i] <- model4$coefficients[2]
  #     sd_Y_list4[i] <- summary(model4)$coefficients[,'Std. Error'][2]
  #   }else{
  #     model4_1 <- lm(Y_0~X_hat_matrix_AIPW[,i] + C_0)  #无Z
  #     b_Y_list4[i] <- model4_1$coefficients[2]
  #     sd_Y_list4[i] <- summary(model4_1)$coefficients[,'Std. Error'][2]
  #   }
  # }
  # sd_Y_list4_inv <- 1/sd_Y_list4
  # beta_hat_4 <- sum(sd_Y_list4_inv^2 * b_Y_list4)/sum(sd_Y_list4_inv^2)
}

#----------------------------结果保存
#beta_hat_vec <- c(beta_hat_1, beta_hat_2, beta_hat_3,beta_hat_4)
beta_hat_vec <- c(beta_hat_1, beta_hat_2)
#--------------------------------------------------IVW方差计算
# library(Rmpfr)
# library(Ryacas0)
sd_OR <- sqrt(1/sum(sd_Y_list2_inv^2))
# sd_IPW <- sqrt(1/sum(sd_Y_list3_inv^2))
# sd_AIPW <- sqrt(1/sum(sd_Y_list4_inv^2))

p_power_OR <- 2 * (1 - pnorm(abs(beta_hat_2 / sd_OR)))
# p_power_IPW <- 2 * (1 - pnorm(abs(beta_hat_3 / sd_IPW)))
# p_power_AIPW <- 2 * (1 - pnorm(abs(beta_hat_4 / sd_AIPW)))
# p_power_mm <- c(p_power_OR,p_power_IPW,p_power_AIPW)
p_power_mm <- c(p_power_OR)
p_power_OR
p_vec <- c(beta_hat_1_list$p,p_power_mm)
p_vec


#------------------手动存数据
data_result1 <- matrix(NA,nrow = 9,ncol = 7)
colnames(data_result1) <- c('method','sample','outcome','beta','p','d1','d2')
data_result1 <- as.data.frame(data_result1)

# data_result1 <- data_result[1:9,]
data_result1$method <- method_names
data_result1$sample <- 1
data_result1$outcome <- 'FEV1/FVC'  #'FEV1/FVC'
data_result1$beta <- beta_hat_vec
data_result1$p <- p_vec

data_result <- rbind(data_result,data_result1)

##
##
##
data_result[data_result$method==method_names[1],]
#------------------
#cbind(b_Y_list2, sd_Y_list2,p_1)[order(sd_Y_list2), ]

#sd_hat_vec <- c(beta_hat_1_list$se,sd_OR,sd_IPW,sd_AIPW)
sd_hat_vec <- c(beta_hat_1_list$se,sd_OR)
p_vec <- c(beta_hat_1_list$p,p_power_mm)
uci_vec <- beta_hat_vec + 1.95 * sd_hat_vec
lci_vec <- beta_hat_vec - 1.95 * sd_hat_vec

beta_hat_vec_t <- round(beta_hat_vec,3)
sd_hat_vec_t <-round(sd_hat_vec,4)
uci_vec_t <- round(uci_vec,4)
lci_vec_t <- round(lci_vec,4)


betaCI <- paste0( beta_hat_vec_t,'(',lci_vec_t,',',uci_vec_t,')'  )
betaCI2 <- paste0('[',lci_vec_t,',',uci_vec_t,']'  )
#-----OR值
beta_OR <- exp(beta_hat_vec)
sd_OR <- beta_OR * sd_hat_vec
uci_OR <- exp(uci_vec)
lci_OR <- exp(lci_vec)

beta_OR <- round(beta_OR,3)
sd_OR <-round(sd_OR,4)
uci_OR <- round(uci_OR,4)
lci_OR <- round(lci_OR,4)
ORCI <- paste0( beta_OR,'(',lci_OR,',',uci_OR,')'  )
ORCI2 <- paste0( '[',lci_OR,',',uci_OR,']'  )
#method_names <- c("IVW","Egger","Weighted median","Weighted mode","MR-RAPS","Contamination Mixture", "MR-Robust", "MR-Lasso",'Trans-OR','Trans-IPW','Trans-AIPW')
method_names <- c("IVW","Egger","Weighted median","Weighted mode","MR-RAPS","MR-Conmix", "MR-Robust", "MR-Lasso",'TLMR')

target_population <- rep('1',length(beta_hat_vec))
target_population <- rep('0',length(beta_hat_vec))

dt_forest <- data.frame(p=round(p_vec,3)
                        ,OR=beta_OR
                        ,loci_OR=lci_OR
                        ,upci_OR=uci_OR
                        ,sd_OR=sd_OR
                        ,ORCI=ORCI
                        ,ORCI2=ORCI2
                        ,beta=beta_hat_vec_t
                        ,loci=lci_vec_t
                        ,upci=uci_vec_t
                        ,sd=sd_hat_vec_t
                        ,betaCI=betaCI
                        ,betaCI2=betaCI2
                        ,Method = method_names
                        ,Exposure = 'BMI'
                        ,Outcome = 'Pulmonary function (FEV1)'
                        ,target = target_population
)

dt_forest_save1 <- dt_forest   #1
dt_forest_save2 <- dt_forest   #0

dt_forest_all <- rbind(dt_forest_save1,dt_forest_save2)

dt_forest_all_sorted <- dt_forest_all %>%
  mutate(Method = factor(Method, levels = method_names)) %>%
  arrange(Method, factor(target, levels = c("north", "south")))  # 按 target 排序

dt_forest <- dt_forest_all_sorted

# levels(dt_forest$Method) <- c(levels(dt_forest$Method), "MR-Conmix")
# #dt_forest$Method[is.na(dt_forest$Method)] <- 'MR-Conmix'
# dt_forest$Method[dt_forest$Method=='Contamination Mixture'] <- 'MR-Conmix'

dt_forest<-dt_forest %>% select(
  Exposure,Outcome,target,Method,beta,betaCI2,sd,p,loci,upci
)

#
#
#
#到此为止
#------------------------------------------------------------------------------------------


dt_forest <- dt_forest %>% 
  mutate(empty_column = '                                 ') %>%  # 创建空列
  select(Exposure,Outcome,target,Method,empty_column,beta,betaCI2,sd,p,loci,upci) %>%
  rename(` ` = empty_column)  # 将列名改为 "空格"


library('forestploter')   
tm <- forest_theme(base_size = 10,  #文本的大小
                   # Confidence interval point shape, line type/color/width
                   ci_pch = 20,   #可信区间点的形状
                   ci_col = "#762a83",    #CI的颜色
                   ci_fill = "blue",     #ci颜色填充
                   ci_alpha = 0.8,        #ci透明度
                   ci_lty = 1,            #CI的线型
                   ci_lwd = 1.5,          #CI的线宽
                   ci_Theight = 0.2, # Set an T end at the end of CI  ci的高度，默认是NULL
                   # Reference line width/type/color   参考线默认的参数，中间的竖的虚线
                   refline_lwd = 1,       #中间的竖的虚线
                   refline_lty = "dashed",
                   refline_col = "grey20",
                   # Vertical line width/type/color  垂直线宽/类型/颜色   可以添加一条额外的垂直线，如果没有就不显示
                   vertline_lwd = 1,              #可以添加一条额外的垂直线，如果没有就不显示
                   vertline_lty = "dashed",
                   vertline_col = "grey20",
                   # Change summary color for filling and borders   更改填充和边框的摘要颜色
                   summary_fill = "yellow",       #汇总部分大菱形的颜色
                   summary_col = "#4575b4",
                   # Footnote font size/face/color  脚注字体大小/字体/颜色
                   footnote_cex = 0.6,
                   footnote_fontface = "italic",
                   footnote_col = "red")+ 
  theme(plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm"))




x_limits_1 <- as.vector(
  c(floor(quantile(as.numeric(dt_forest$loci),0.05)*10)/10-0.1
    ,ceiling(quantile(as.numeric(dt_forest$upci),0.95)*10)/10+0.1)
)
x_step_1 <- floor((x_limits_1[2]-x_limits_1[1])/5 * 10)/10
x_ticks_1 <- seq(x_limits_1[1],x_limits_1[2],x_step_1)

colnames(dt_forest) <- c("Exposure", "Outcome", "Target population", "Method",   " "  ,      "Beta",     "95%CI",  "SE",      
                         "Pvalue"   ,     "loci"    , "upci"   )
dt_forest$Pvalue[dt_forest$Pvalue == 0] <- '<0.001'

plot_forest <- forestploter::forest(dt_forest[,c(1:9)]
                                    ,est = as.numeric(dt_forest$Beta)
                                    ,lower = as.numeric(dt_forest$loci)
                                    ,upper = as.numeric(dt_forest$upci)
                                    ,sizes = as.numeric(dt_forest$SE)
                                    ,ci_column =5
                                    ,ref_line = 0
                                    ,arrow_lab = c("Negative", "Positive")
                                    ,xlim = x_limits_1
                                    ,ticks_at=c(-0.5,0,0.4)
                                    #,footnote = "This is the demo data, Please feel free to change\nanything you want."
                                    ,theme = tm
                                    #,column_widths =c (1,1,1,3,1,1)
)
plot_forest
# ggsave(paste('森林图_全人群',date_mark,'_',ver,'.png',sep = ''),plot = plot_forest, bg="white", path = path_output_plots,width = 20,height = 9,dpi = 600)
# ggsave(paste('森林图_全人群',date_mark,'_',ver,'.pdf',sep = ''),plot = plot_forest, bg="white", path = path_output_plots,width = 20,height = 9,dpi = 600)

ggsave(paste('森林图_分组250304_',date_mark,'_',ver,'.png',sep = ''),plot = plot_forest, bg="white", path = path_output_plots, height = 13.5,width = 18,dpi = 500)
ggsave(paste('森林图_分组250304_',date_mark,'_',ver,'.pdf',sep = ''),plot = plot_forest, bg="white", path = path_output_plots, height = 13.5,width = 18,dpi = 500)
#--------------------------------------------------

#p_cover_vec <- c(beta_hat_1_list$p_cover,p_cover_mm) 
dt_forest <- data.frame()

