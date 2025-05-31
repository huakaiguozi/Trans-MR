
n_divides <- 3
indicator_print<-1

path_data_read <- 'data/较高收入人群_10W_241112.csv'
data_1 <- fread(path_data_read)
path_data_read2 <- 'data/较低收入人群_10W_241112.csv'
data_0 <- fread(path_data_read2)

# 使用 grep 查找所有 rs 开头的列
rs_columns <- grep("^rs", colnames(data_1), value = TRUE)

# 提取这些列的数据
Z_1 <- data_1[, ..rs_columns]
Z_1 <- as.matrix(Z_1)

Z_0 <- data_0[, ..rs_columns]
Z_0 <- as.matrix(Z_0)

data_12 <- data_1 %>% select(-all_of(rs_columns))
data_02 <- data_0 %>% select(-all_of(rs_columns))

C_1 <- data_12[,c(2:ncol(data_12)),with = FALSE]
C_1 <- C_1 %>% select(-all_of(c("种族背景",'BMI','第1秒用力呼气容积(FEV1)','收入得分')))
AS1 <- C_1[,c('性别','年龄')]

C_0 <- data_02[,c(2:ncol(data_02)),with = FALSE]
C_0 <- C_0 %>% select(-all_of(c("种族背景",'BMI','第1秒用力呼气容积(FEV1)','收入得分')))
AS0 <- C_0[,c('性别','年龄')]

quantile_vec_all <- seq(0,1,1/n_divides)
quantile_vec <- quantile_vec_all[-c(1,length(quantile_vec_all))]

#------------------------------
C_1 <- as.data.frame(C_1)
C_0 <- as.data.frame(C_0)
C_all <- rbind(C_0,C_1)

C_0_cl <- C_0
C_1_cl <- C_1
for (i in c(1:dim(C_0)[2])) {
  if( length(unique(C_0[,i])) > 5 ){
    
    #C_qt_vector <- quantile(C_all[,i],quantile_vec)
    C_qt_vector <- unique(quantile(C_1[, i], quantile_vec))
    # 检查去重后是否有足够的分位数
    if (length(C_qt_vector) >= n_divides - 1) {
      # 使用去重后的分位数作为 cut 的 breaks 参数
      C_qt_label <- cut(C_0[, i], breaks = c(-Inf, C_qt_vector, Inf), 
                        labels = 1:n_divides, include.lowest = TRUE)
      C_0_cl[, i] <- C_qt_label
      C_qt_label <- cut(C_1[, i], breaks = c(-Inf, C_qt_vector, Inf), 
                        labels = 1:n_divides, include.lowest = TRUE)
      C_1_cl[, i] <- C_qt_label
    } else {
      C_qt_label <- cut(C_0[, i], breaks = c(-Inf, C_qt_vector, Inf), 
                        labels = 1:(length(C_qt_vector)+1), include.lowest = TRUE)
      C_0_cl[, i] <- C_qt_label
      C_qt_label <- cut(C_1[, i], breaks = c(-Inf, C_qt_vector, Inf), 
                        labels = 1:(length(C_qt_vector)+1), include.lowest = TRUE)
      C_1_cl[, i] <- C_qt_label
      warning(paste("第", i, "列的分位数不足以进行划分,","现分为",(length(C_qt_vector)+1),"类"))
    }
  }
}
# 替换成分类变量；原始数据存储至origin
C_0_origin <- C_0
C_1_origin <- C_1
C_0 <- C_0_cl
C_1 <- C_1_cl

#------------------------------
# C_1 <- as.data.frame(C_1)
# C_1_cl <- C_1
# for (i in c(1:dim(C_1)[2])) {
#   #print(i)
#   if( length(unique(C_1[,i])) > 30 ){
#     #C_qt_vector <- quantile(C_1[,i],quantile_vec)
#     #C_qt_label <- cut( C_0[,i],c(-Inf,C_qt_vector,Inf),labels=c(1:n_divides) )
#     #C_0_cl[,i] <- C_qt_label
#     # 计算分位数
#     C_qt_vector <- unique(quantile(C_1[, i], quantile_vec))
#     
#     # 检查去重后是否有足够的分位数
#     if (length(C_qt_vector) >= n_divides - 1) {
#       # 使用去重后的分位数作为 cut 的 breaks 参数
#       C_qt_label <- cut(C_1[, i], breaks = c(-Inf, C_qt_vector, Inf), 
#                         labels = 1:n_divides, include.lowest = TRUE)
#       C_1_cl[, i] <- C_qt_label
#     } else {
#       C_qt_label <- cut(C_1[, i], breaks = c(-Inf, C_qt_vector, Inf), 
#                         labels = 1:(length(C_qt_vector)+1), include.lowest = TRUE)
#       C_1_cl[, i] <- C_qt_label
#       
#       warning(paste("第", i, "列的分位数不足以进行划分,","现分为",(length(C_qt_vector)+1),"类"))
#     }
#     
#     
#     # C_qt_label <- cut( C_1[,i],c(-Inf,C_qt_vector,Inf),labels=c(1:n_divides),include.lowest=TRUE )
#     # C_1_cl[,i] <- C_qt_label
#   }
# }
# # 替换成分类变量；原始数据存储至origin
# #C_0_origin <- C_0
# C_1_origin <- C_1
# #C_0 <- C_0_cl
# C_1 <- C_1_cl
C_1 <- as.matrix(C_1)
C_1_cl <- as.matrix(C_1_cl)
C_0 <- as.matrix(C_0)
C_0_cl <- as.matrix(C_0_cl)
#------------------------------
p_Z <- ncol(Z_1)
p_C <- ncol(C_1)

C_1 <- apply(C_1, 2, function(x) as.numeric(as.character(x)))
Z_1 <- apply(Z_1, 2, function(x) as.numeric(as.character(x)))
C_1_origin <- apply(C_1_origin, 2, function(x) as.numeric(as.character(x)))

C_0 <- apply(C_0, 2, function(x) as.numeric(as.character(x)))
Z_0 <- apply(Z_0, 2, function(x) as.numeric(as.character(x)))
C_0_origin <- apply(C_0_origin, 2, function(x) as.numeric(as.character(x)))
#C_0_cl <- apply(C_0_cl, 2, function(x) as.numeric(as.character(x)))
################################################################################
#------------------------------
X_1 <- data_1[["BMI"]]
X_1 <- matrix(X_1,nrow = nrow(data_1))
#
res_1 <- summary_compute(Z=Z_1,X=X_1,C=AS1)
beta_1 <- res_1$beta
sd_1 <- res_1$sd
p_1 <- res_1$p
eaf_1 <- res_1$eaf
R2_1 <- res_1$R2
F_stats_1 <- res_1$F_stats

data_1 <- list(Z=Z_1,C=C_1,X=X_1,beta=beta_1,sd=sd_1,p=p_1,eaf=eaf_1,R2=R2_1,F_stats=F_stats_1)
################################################################################
#------------------------------
X_0 <- data_0[["BMI"]]
X_0 <- matrix(X_0,nrow = nrow(data_0))
Y_0 <- data_0[["第1秒用力呼气容积(FEV1)"]]
Y_0 <- matrix(Y_0,nrow = nrow(data_0))

#
res_0 <- summary_compute(Z=Z_0,X=Y_0,C=AS0)
beta_0 <- res_0$beta
sd_0 <- res_0$sd
p_0 <- res_0$p
eaf_0 <- res_0$eaf
R2_0 <- res_0$R2
F_stats_0 <- res_0$F_stats
#------------------
# Z_0<-Z_1
# C_0<-C_1
# C_0_origin <-C_1_origin
#------------------
data_0 <- list(Z=Z_0,C=C_0,X=X_0, Y=Y_0,beta=beta_0,sd=sd_0,p=p_0,eaf=eaf_0,R2=R2_0,F_stats=F_stats_0)
################################################################################
beta_hat_1_list <- Twosamle_package(p_Z=p_Z,p_C=p_C,data_0=data_0,data_1=data_1,beta=1)
beta_hat_1 <- beta_hat_1_list$beta

# #----------------------------看是否存在修饰作用的代码
C_all <- rbind(C_0,C_1)
X_all <- rbind(X_0,X_1)
Z_all <- rbind(Z_0,Z_1)
#-----------
modifier_matrix <- matrix(nrow = p_Z,ncol = p_C)
modifier_p <- matrix(nrow = p_Z,ncol = p_C)
for(k in 1:p_Z){
  for (j in 1:p_C) {
    data_chow_test <- data.frame(Z=Z_all[,k],X=X_all,C=C_all[,j])
    data_chow_test <- data_chow_test[order(data_chow_test$C), ]
    C_group_vec <- table(data_chow_test$C)

    end_indice <- cumsum(C_group_vec)
    start_indice <- c(1,head(end_indice,-1)+1)

    beta_group_vec <- c()
    se_group_vec <- c()
    for(i in 1:length(C_group_vec)){
      print(c(k,j,i))
      data_chow_test_sub <- data_chow_test[(start_indice[i]:end_indice[i]),]
      if(nrow(data_chow_test_sub) < 2){
        next
      }
      model_group_C <- lm(X~Z,data = data_chow_test_sub)
      res_model_group <- summary(model_group_C)$coefficients
      if(nrow(res_model_group)<2){
        next
      }
      beta_group_vec <- c(beta_group_vec,res_model_group[2,1])
      se_group_vec <- c(se_group_vec,res_model_group[2,2])
    }
    res_meta_analysis <- metafor::rma.uni(yi=beta_group_vec,sei = se_group_vec,method = 'REML')
    Q_statistics_group <- res_meta_analysis$QE
    p_Q_statistics_group <- res_meta_analysis$QEp
    modifier_p[k,j] <- p_Q_statistics_group
    if(p_Q_statistics_group <0.05/p_C){
      modifier_matrix[k,j] <- 1
    }else{
      modifier_matrix[k,j] <- 0
    }
  }
}
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

#--------------


#-----------


#------
indice_Cnouse <- which(colSums(modifier_matrix) <5)
#16 睡眠市场
indice_Cnouse <-c(1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 13, 14, 16, 17, 18, 19, 20, 21)
C_1 <- C_1[,-indice_Cnouse]
C_1_origin <- C_1_origin[,-indice_Cnouse]
C_0 <- C_0[,-indice_Cnouse]
C_0_origin <- C_0_origin[,-indice_Cnouse]
#----------------------------
b_Y_list2 <- c()
sd_Y_list2 <- c()
b_Y_list3 <- c()
sd_Y_list3 <- c()
b_Y_list4 <- c()
sd_Y_list4 <- c()

modifier_matrix_save <- modifier_matrix
modifier_matrix <- modifier_matrix[,setdiff(c(1:21),indice_Cnouse)]
p_C <- ncol(C_1)

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
    X_hat_matrix_OR <- Trans_OLS2(Z_0=Z_0m,Z_1=Z_1m,C_0=C_0_origin,C_1=C_1_origin,X_1=X_1)$X_hat_matrix
    #Trans_IPW
    X_hat_matrix_IPW <- IPW2(Z_0=Z_0m,Z_1=Z_1m,C_0=C_0,C_1=C_1,X_1=X_1)$X_hat_matrix
    #Trans_AIPW
    X_hat_matrix_AIPW <- AIPW2(Z_0=Z_0m,Z_1=Z_1m,C_0=C_0,C_1=C_1,X_1=X_1
                               ,C_0_origin=C_0_origin,C_1_origin=C_1_origin)$X_hat_matrix
    
    write.csv(X_hat_matrix_OR,file=paste(path_output_tables, '/X_OR2_',as.character(n_divides),'_',date_mark,'_',ver,'.csv',sep = ''))
    write.csv(X_hat_matrix_IPW,file=paste(path_output_tables, '/X_IPW2__',as.character(n_divides),'_',date_mark,'_',ver,'.csv',sep = ''))
    write.csv(X_hat_matrix_AIPW,file=paste(path_output_tables, '/X_AIPW2_',as.character(n_divides),'_',date_mark,'_',ver,'.csv',sep = ''))
    
    #-----------------------------------------
    
    #---------------------------------------------------------------Trans_OR0         21
    if(p_Zm >= 1){
      for (i in 1:p_Zm) {  
        model2 <- lm(Y_0~X_hat_matrix_OR[,i] + C_0 + Z_0m[,i])  #有Z
        summary_model2 <- summary(model2)
        pval_model2 <- summary_model2$coefficients[,4]
        Z_pval <- pval_model2[length(pval_model2)]
        if(Z_pval < pval_i){
          b_Y_list2[i] <- model2$coefficients[2]
          sd_Y_list2[i] <- summary(model2)$coefficients[,'Std. Error'][2]
        }else{
          model2_1 <- lm(Y_0~X_hat_matrix_OR[,i] + C_0)  #无Z
          b_Y_list2[i] <- model2_1$coefficients[2]
          sd_Y_list2[i] <- summary(model2_1)$coefficients[,'Std. Error'][2]
        }
      }
    }
    b_Y_list2 <- c(b_Y_list2, beta_hat_degrade)
    sd_Y_list2 <- c(sd_Y_list2, se_degrade)
    
    sd_Y_list2_inv <- 1/sd_Y_list2 
    beta_hat_2 <- sum(sd_Y_list2_inv^2 * b_Y_list2)/sum(sd_Y_list2_inv^2)
    #---------------------------------------------------------------Trans_IPW0          31
    if(p_Zm>=1){
      
      for (i in 1:p_Zm) {
        
        model3 <- lm(Y_0~X_hat_matrix_IPW[,i] + C_0 + Z_0m[,i])  #有Z
        summary_model3 <- summary(model3)
        pval_model3 <- summary_model3$coefficients[,4]
        Z_pval <- pval_model3[length(pval_model3)]
        if(Z_pval < pval_i){
          b_Y_list3[i] <- model3$coefficients[2]
          sd_Y_list3[i] <- summary(model3)$coefficients[,'Std. Error'][2]
        }else{
          model3_1 <- lm(Y_0~X_hat_matrix_IPW[,i] + C_0 + Z_0[,i])  #无Z
          b_Y_list3[i] <- model3_1$coefficients[2]
          sd_Y_list3[i] <- summary(model3_1)$coefficients[,'Std. Error'][2]
        }
      }
    }
    
    b_Y_list3 <- c(b_Y_list3, beta_hat_degrade)
    sd_Y_list3 <- c(sd_Y_list3, se_degrade)
    
    sd_Y_list3_inv <- 1/sd_Y_list3 
    beta_hat_3 <- sum(sd_Y_list3_inv^2 * b_Y_list3)/sum(sd_Y_list3_inv^2)
    #---------------------------------------------------------------Trans_AIPW0          41
    if(p_Zm>=1){
      
      for (i in 1:p_Zm) {
        model4 <- lm(Y_0~X_hat_matrix_AIPW[,i] + C_0 + Z_0m[,i])  #有Z
        summary_model4 <- summary(model4)
        pval_model4 <- summary_model4$coefficients[,4]
        Z_pval <- pval_model4[length(pval_model4)]
        if(Z_pval<pval_i){
          b_Y_list4[i] <- model4$coefficients[2]
          sd_Y_list4[i] <- summary(model4)$coefficients[,'Std. Error'][2]
        }else{
          model4_1 <- lm(Y_0~X_hat_matrix_AIPW[,i] + C_0 + Z_0[,i])  #无Z
          b_Y_list4[i] <- model4_1$coefficients[2]
          sd_Y_list4[i] <- summary(model4_1)$coefficients[,'Std. Error'][2]
        }
      }
    }
    b_Y_list4 <- c(b_Y_list4, beta_hat_degrade)
    sd_Y_list4 <- c(sd_Y_list4, se_degrade)
    
    sd_Y_list4_inv <- 1/sd_Y_list4
    beta_hat_4 <- sum(sd_Y_list4_inv^2 * b_Y_list4)/sum(sd_Y_list4_inv^2)
    
  }
}else{
  #===============================================================迁移估计量
  #Trans_OR
  X_hat_matrix_OR <- Trans_OLS2(Z_0=Z_0,Z_1=Z_1,C_0=C_0_origin,C_1=C_1_origin,X_1=X_1)$X_hat_matrix
  #Trans_IPW
  X_hat_matrix_IPW <- IPW2(Z_0=Z_0,Z_1=Z_1,C_0=C_0,C_1=C_1,X_1=X_1)$X_hat_matrix
  #Trans_AIPW
  X_hat_matrix_AIPW <- AIPW2(Z_0=Z_0,Z_1=Z_1,C_0=C_0,C_1=C_1,X_1=X_1
                             ,C_0_origin=C_0_origin,C_1_origin=C_1_origin)$X_hat_matrix
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
}

#----------------------------结果保存
beta_hat_vec <- c(beta_hat_1, beta_hat_2, beta_hat_3,beta_hat_4)

#--------------------------------------------------IVW方差计算
# library(Rmpfr)
# library(Ryacas0)
sd_OR <- sqrt(1/sum(sd_Y_list2_inv^2))
sd_IPW <- sqrt(1/sum(sd_Y_list3_inv^2))
sd_AIPW <- sqrt(1/sum(sd_Y_list4_inv^2))
#p_power_OR <- 2 * (1 - pnorm(abs(beta_hat_2 / sd_OR)))

p_power_OR <- 2 * (1 - pnorm(abs(beta_hat_2 / sd_OR)))
p_power_IPW <- 2 * (1 - pnorm(abs(beta_hat_3 / sd_IPW)))
p_power_AIPW <- 2 * (1 - pnorm(abs(beta_hat_4 / sd_AIPW)))
p_power_mm <- c(p_power_OR,p_power_IPW,p_power_AIPW)

sd_hat_vec <- c(beta_hat_1_list$se,sd_OR,sd_IPW,sd_AIPW)
p_vec <- c(beta_hat_1_list$p,p_power_mm)
uci_vec <- beta_hat_vec + 1.95 * sd_hat_vec
lci_vec <- beta_hat_vec - 1.95 * sd_hat_vec
#-----OR
beta_OR <- exp(beta_hat_vec)
sd_OR <- beta_OR * sd_hat_vec
uci_OR <- exp(uci_vec)
lci_OR <- exp(lci_vec)

beta_OR <- round(beta_OR,3)
sd_OR <-round(sd_OR,4)
uci_OR <- round(uci_OR,4)
lci_OR <- round(lci_OR,4)
ORCI <- paste0( beta_OR,'(',lci_OR,',',uci_OR,')'  )

method_names <- c("IVW","Egger","Weighted median","Weighted mode","MR-RAPS","Contamination Mixture", "MR-Robust", "MR-Lasso",'Trans-OR','Trans-IPW','Trans-AIPW')
dt_forest <- data.frame(p=round(p_vec,3)
                        ,OR=beta_OR
                        ,loci_OR=lci_OR
                        ,upci_OR=uci_OR
                        ,se=sd_OR
                        ,ORCI=ORCI
                        ,Method = method_names
                        ,Exposure = 'BMI'
                        ,Outcome = 'Pulmonary function (FEV1)'
)
dt_forest<-dt_forest %>% select(
  Exposure,Outcome,Method,ORCI,p,OR,loci_OR,upci_OR,se
)
dt_forest <- dt_forest %>% 
  mutate(empty_column = '                                 ') %>%  # 创建空列
  select(Exposure, Outcome, Method, empty_column, ORCI, p, OR, loci_OR, upci_OR, se) %>%
  rename(` ` = empty_column)  # 将列名改为 "空格"


library('forestploter')   
tm <- forest_theme(base_size = 10,  #文本的大小
                   # Confidence interval point shape, line type/color/width
                   ci_pch = 15,   #可信区间点的形状
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
                   footnote_col = "red")

x_limits_1 <- as.vector(
  c(floor(quantile(as.numeric(dt_forest$loci_OR),0.05)*10)/10-0.2
    ,ceiling(quantile(as.numeric(dt_forest$upci_OR),0.95)*10)/10+0.2)
)
x_step_1 <- floor((x_limits_1[2]-x_limits_1[1])/5 * 10)/10
x_ticks_1 <- seq(x_limits_1[1],x_limits_1[2],x_step_1)

plot_forest <- forestploter::forest(dt_forest[,c(1:6)]
                                    ,est = as.numeric(dt_forest$OR)
                                    ,lower = as.numeric(dt_forest$loci_OR)
                                    ,upper = as.numeric(dt_forest$upci_OR)
                                    ,sizes = as.numeric(dt_forest$se)
                                    ,ci_column =4
                                    ,ref_line = 1
                                    ,arrow_lab = c("Negative", "Positive")
                                    ,xlim = x_limits_1
                                    ,ticks_at=c(0.7,1.0,1.4)
                                    #,footnote = "This is the demo data, Please feel free to change\nanything you want."
                                    ,theme = tm
                                    #,column_widths =c (1,1,1,3,1,1)
)
plot_forest
ggsave(paste('森林图_低收入人群',date_mark,'_',ver,'.png',sep = ''),plot = plot_forest, bg="white", path = path_output_plots,width = 20,height = 9,dpi = 600)
ggsave(paste('森林图_低收入人群',date_mark,'_',ver,'.pdf',sep = ''),plot = plot_forest, bg="white", path = path_output_plots,width = 20,height = 9,dpi = 600)

#--------------------------------------------------

#p_cover_vec <- c(beta_hat_1_list$p_cover,p_cover_mm) 
dt_forest <- data.frame()

