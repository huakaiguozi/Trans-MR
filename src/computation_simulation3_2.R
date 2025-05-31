#time_start <- Sys.time()
# 设置并行计算的核心数量
cores <- detectCores()
registerDoParallel(cores = cores-2)

#参数空间
pi_data1 <- matrix(round(with_seed(1, runif(90000,0,2)),2),nrow = 300)
pi_data2 <- matrix(round(with_seed(1, runif(90000,0,2)),2),nrow = 300)
omega_data <- round(with_seed(1, runif(90000,0,2)),2)
alpha_data <- round(with_seed(1,runif(1000,0.2,0.8)),2)
gamma_data <- round(with_seed(1,runif(1000,0.2,0.4)),2)
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

inter_strength <- 0
alpha_modifier <- 1
eta_modifier <- 1

ratio_pleiotropy <- 1

#--------------------------
p_Z <- 1
p_Z_max =1
#-----------
p_C <- 2

eta_1 <- eta1_data[1:p_C] * eta_modifier
eta_2 <- eta2_data[1:p_C] * eta_modifier

#-----
#alpha <- seq(2,4,2/(p_Z-1))*alpha_modifier
alpha <- alpha_data[1:p_Z_max] * alpha_modifier
pi <- rep(1,p_C) * inter_strength
pi2 <- pi_data2[1:p_Z_max,1:p_C] * inter_strength
omega <- omega_data[1:p_C]

gamma <- gamma_data[1:p_Z_max] * gamma_multi
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
param_Y <- list(beta=beta,eta_2=eta_2,gamma=gamma,ratio_pleiotropy=ratio_pleiotropy
                ,heter_level=heter_level,omega=omega
                #,U_md = U_md
)

pval_i <- 0.05

#===========================================================================
beta_hat_ite <- foreach(ite_k=1:ite_times, .export =c('heter_level')
                        , .packages = c('fixest','TwoSampleMR','glmnet','mr.raps','penalized','ridge','FNN','data.table'))%dopar%
  {
    #for(ite_k in 1:ite_times){
    tryCatch({
      #indicator_print <- 1
      
      seed <- ite_k * 100
      set.seed(seed)
      #-----------------------------------------------生成随机数
      data_0 <- data_generate_Y(n=n_0,p_Z=p_Z,p_Z_max = p_Z_max
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
          data_0 <- data_generate_Y(n=n_0,p_Z=p_Z,p_Z_max = p_Z_max
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
      
      #Trans_OR
      X_hat_matrix_OR <- Trans_OLS2(Z_0=Z_0,Z_1=Z_1,C_0=C_0,C_1=C_1,X_1=X_1)$X_hat_matrix
      #Trans_IPW
      X_hat_matrix_IPW <- IPW2(Z_0=Z_0,Z_1=Z_1,C_0=C_0,C_1=C_1,X_1=X_1)$X_hat_matrix
      #Trans_AIPW
      X_hat_matrix_AIPW <- AIPW2(Z_0=Z_0,Z_1=Z_1,C_0=C_0,C_1=C_1,X_1=X_1
      )$X_hat_matrix
      
      # Xinter_OR <- C_0 * as.numeric(X_hat_matrix_OR)
      # Xinter_IPW <- C_0 * as.numeric(X_hat_matrix_IPW) 
      # Xinter_AIPW <- C_0 * as.numeric(X_hat_matrix_AIPW) 
      # 
      # model2_1 <- lm(Y_0~Xinter_OR)  #无Z
      # model2_2 <- lm(Y_0~Xinter_IPW)
      # model2_3 <- lm(Y_0~Xinter_AIPW)
      
      # model2_1 <- lm(Y_0~X_hat_matrix_OR + C_0+ X_hat_matrix_OR * C_0)  #无Z
      # model2_2 <- lm(Y_0~X_hat_matrix_IPW + C_0+ X_hat_matrix_IPW * C_0)
      # model2_3 <- lm(Y_0~X_hat_matrix_AIPW + C_0+ X_hat_matrix_AIPW * C_0)
      
      model2_1 <- lm(Y_0~X_hat_matrix_OR + C_0+ Z_0+ X_hat_matrix_OR * C_0 )  #有Z
      model2_2 <- lm(Y_0~X_hat_matrix_IPW + C_0+ Z_0+ X_hat_matrix_IPW * C_0)
      model2_3 <- lm(Y_0~X_hat_matrix_AIPW + C_0+ Z_0+ X_hat_matrix_AIPW * C_0)
      
      omega_1 <- model2_1$coefficients[c((length(model2_1$coefficients)-p_C+1):length(model2_1$coefficients))]
      omega_2 <- model2_2$coefficients[c((length(model2_1$coefficients)-p_C+1):length(model2_1$coefficients))]
      omega_3 <- model2_3$coefficients[c((length(model2_1$coefficients)-p_C+1):length(model2_1$coefficients))]
      
      #omega_hat <- c(omega_1,omega_2,omega_3)
      
      return_list <- list(omega1 = omega_1,omega2=omega_2,omega3=omega_3)
      
    }, error = function(e) {
      #warnings(paste("Iteration", ite_k, "failed:", e$message))
      return(NULL)
      # message_list <- list(message = paste("Iteration", ite_k, "failed:", e$message)) 
      # return(message_list)
    })
    
  }
gc()

# 构建空列表用于存储各个 omega 值
omega1_list <- list()
omega2_list <- list()
omega3_list <- list()

# 遍历 beta_hat_ite，提取 omega1, omega2, omega3
for (i in 1:length(beta_hat_ite)) {
  omega1_list[[i]] <- as.numeric(beta_hat_ite[[i]]$omega1)
  omega2_list[[i]] <- as.numeric(beta_hat_ite[[i]]$omega2)
  omega3_list[[i]] <- as.numeric(beta_hat_ite[[i]]$omega3)
}

# 将列表拼接成矩阵，每行为一次取值
omega1_matrix <- do.call(rbind, omega1_list)
omega2_matrix <- do.call(rbind, omega2_list)
omega3_matrix <- do.call(rbind, omega3_list)

apply(sweep(omega1_matrix, 2, omega, FUN = "-"),2,mean)
#gamma 0.5   -0.04468899 -0.14231651
#      0.25  -0.02296020 -0.07177885
summary(omega1_matrix)
summary(omega2_matrix)

bias_1 <- apply(sweep(omega1_matrix, 2, omega, FUN = "-"),2,mean)
bias_2 <- apply(sweep(omega2_matrix, 2, omega, FUN = "-"),2,mean)
bias_3 <- apply(sweep(omega3_matrix, 2, omega, FUN = "-"),2,mean)

sd_1 <- apply(omega1_matrix, 2, sd)
sd_2 <- apply(omega2_matrix, 2, sd)
sd_3 <- apply(omega3_matrix, 2, sd)

mse_1 <- bias_1^2 + sd_1^2
mse_2 <- bias_2^2 + sd_2^2
mse_3 <- bias_3^2 + sd_3^2

dt1 <- c()

for (i in 1:p_C) {
  dt2 <- c()
  name_vec1 <- rep(paste0('Omega_',i),3)
  name_vec2 <- c('bias','sd','mse')
  num_vec1 <- c(bias_1[i],sd_1[i],mse_1[i])
  num_vec2 <- c(bias_2[i],sd_2[i],mse_2[i])
  num_vec3 <- c(bias_3[i],sd_3[i],mse_3[i])
  dt2 <- rbind(name_vec1,name_vec2,num_vec1,num_vec2,num_vec3)
  dt1 <- cbind(dt1,dt2)
}
rnames <- c('param','index','TLMR-OR','TLMR-IPW','TLMR-AIPW')
dt1 <- cbind(rnames,dt1)
fwrite(dt1,file=paste(path_output_tables, '/mse_pi=0',date_mark,'_',ver,'.csv',sep = ''))
#-----------------------------------------------------------------------

methods_ordered <- factor(c('TLMR-OR','TLMR-IPW','TLMR-AIPW'),levels = c('TLMR-OR','TLMR-IPW','TLMR-AIPW'))

omega.1 <- c(omega1_matrix[,1],omega2_matrix[,1],omega3_matrix[,1])
methods_vec <- c(rep(methods_ordered[1],dim(omega1_matrix)[1]) 
                  ,rep(methods_ordered[2],dim(omega2_matrix)[1])
                  ,rep(methods_ordered[3],dim(omega3_matrix)[1]))

omega.2 <- c(omega1_matrix[,2],omega2_matrix[,2],omega3_matrix[,2])

# methods_vec <- c(rep(methods_ordered[1],dim(omega1_matrix)[1]) 
#                  ,rep(methods_ordered[2],dim(omega2_matrix)[1])
#                  ,rep(methods_ordered[3],dim(omega3_matrix)[1]))

#omega_names <- factor(c('',''),levels = c(''))


data_omega.1 <- data.frame(omega=omega.1,method=methods_vec) 
data_omega.2 <- data.frame(omega=omega.2,method=methods_vec) 


color_setting <- c("#66C2A5", "#FC8D62", "#8DA0CB")
zihao1 <- 15 #x轴标题字号
zihao2 <- 12 #x轴刻度字号
zihao3 <- 15 #图例文本字号
kedu_up <- 0.8
kedu_down <- 0.2
kedu_fen <- 0.1
breaks_set <- round(seq(kedu_down,kedu_up,kedu_fen),1)

hline1 <- omega[1]
hline2 <- omega[2]

P1<-ggplot(data=data_omega.1,aes( y = omega)) +
  geom_boxplot(aes(fill = method),outlier.alpha=0.5) +
  scale_fill_manual(values=color_setting)+
  geom_hline(aes(yintercept = hline1), linetype = "dashed", color = "darkred") +
  theme_bw() +
  scale_y_continuous(limits = c(kedu_down, kedu_up), 
                     breaks = breaks_set)+
  xlab(expression(Omega[1])) +
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


kedu_up2 <- 1
kedu_down2 <- 0.5
kedu_fen2 <- 0.1
breaks_set2 <- round(seq(kedu_down2,kedu_up2,kedu_fen2),1)

P2<-ggplot(data=data_omega.2,aes( y = omega)) +
  geom_boxplot(aes(fill = method),outlier.alpha=0.5) +
  scale_fill_manual(values=color_setting)+
  geom_hline(aes(yintercept = hline2), linetype = "dashed", color = "darkred") +
  theme_bw() +
  scale_y_continuous(limits = c(kedu_down2, kedu_up2), 
                     breaks = breaks_set2)+
  xlab(expression(Omega[2])) +
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

# P1
# P2

P3 <- ggarrange(P1
                ,P2,
                ncol=1,nrow=2,
                labels = c('A','B'),
                # labels.x = 1,
                hjust= c(-6,-5.5),
                vjust = 2,
                #widths = c(1,1),
                common.legend = TRUE,
                legend="right") 

P3

ggsave(paste('pi=0_',date_mark,'_',ver,'.png',sep = ''),plot = P3, bg="white", path = path_output_plots,width = 12,height = 9,dpi = 600)
ggsave(paste('pi=0_',date_mark,'_',ver,'.pdf',sep = ''),plot = P3, bg="white", path = path_output_plots,width = 12,height = 9,dpi = 600)


#===============================================================================绘图1
# time_end <- Sys.time()
# time_consume <- difftime(time_end,time_start)
# time_consume 
stopImplicitCluster()