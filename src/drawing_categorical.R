#--------------------------
beta1 <- 0.05
beta2 <- 0.05


zihao1 <- 13 #x轴标题字号
zihao2 <- 12 #x轴刻度字号
zihao3 <- 13 #图例文本字号

zihao4 <- 3.4 #柱状图上面数字字号
#--------------------------
shapes_1 <- c(1:len_methods)-1

# colour_table1 <- c(  '#e6194B', '#3cb44b', '#ffe119', '#4363d8', '#f58231'
#                      , '#42d4f4', '#f032e6', '#fabed4', '#469990', '#dcbeff'
#                      , '#9A6324', '#fffac8', '#800000', '#aaffc3', '#000075'
#                      , '#a9a9a9', '#ffffff', '#000000')
# 
# colour_table2 <- c(  '#e6194B', '#3cb44b', '#ffe119', '#4363d8', '#f58231'
#                      , '#911eb4', '#42d4f4', '#f032e6', '#bfef45', '#fabed4'
#                      , '#469990', '#dcbeff', '#9A6324', '#fffac8', '#800000'
#                      , '#aaffc3', '#808000', '#ffd8b1', '#000075', '#a9a9a9'
#                      , '#ffffff', '#000000')
# 
# 
# colour_table_select <- c( '#e6194B', '#ffe119', '#4363d8', '#f58231'
#                           , '#911eb4', '#42d4f4', '#f032e6'
#                           , '#469990', '#9A6324', '#fffac8', '#800000'
#                           , '#808000', '#000075', '#a9a9a9'
#                           , '#000000','#fabed4','#dcbeff','#ffd8b1')
# 
# #coloo_1 <- c('#aaffc3','#bfef45','#3cb44b', colour_table_select[1:(len_methods-3)])
# coloo_1 <- c(colour_table2[1:len_methods])
#----------
#------------------
# color_df <- Genshinpalette::colorlist
# "#D5E3ED" "#9BC0E1" "#9A9294" "#463542"  #KAMISATOAYAKA 绫华
# "#F8FAFB" "#BCACBF" "#4A6AAF" "#4BDCFB"  #SANGONOMIIYAKOKOMI 心海
# "#EADDF2" "#A785C5" "#6A488E" "#1F1F56"  #RAIDENSHOGUN 雷神
# "#F5DEDB" "#CEA8A0" "#B84C33" "#442A29"  #YANFEI 烟菲
# "#E1FAFF" "#A6E0F8" "#71A0C9" "#27293C"  #EULA 优拉
# "#AFD5F1" "#276BCB" "#1A439D" "#141F45"  #YELAN 夜兰
# "#D0E6F1" "#5C86C6" "#32355A" "#BA85B3"  #QIQI 七七
# "#FAE1BE" "#EF9A56" "#B46349" "#452821"  #HUTAO 胡桃
# "#FFFFF8" "#DEA573" "#B7421E" "#471C12"  #可莉

# color_vec_1 <- c(Genshinpalette::Genshinpalette('KAMISATOAYAKA',4)
#                  ,Genshinpalette::Genshinpalette('SANGONOMIIYAKOKOMI',4)
#                  ,Genshinpalette::Genshinpalette('HUTAO',(len_methods- len_other_methods))
# )
# color_vec_2 <- c(RColorBrewer::brewer.pal(8,"Pastel2"),RColorBrewer::brewer.pal((len_methods- len_other_methods),"Accent"))
# color_vec_3 <- c(RColorBrewer::brewer.pal(8,"Pastel2"),Genshinpalette::Genshinpalette('HUTAO',(len_methods- len_other_methods+1))[-1])
# color_vec_4 <- c(RColorBrewer::brewer.pal(8,"Pastel2"),RColorBrewer::brewer.pal((len_methods- len_other_methods),"Pastel1"))
# color_vec_5 <- c(RColorBrewer::brewer.pal(8,"Pastel2"),RColorBrewer::brewer.pal((len_methods- len_other_methods),"Set3"))
# 
# color_vec_6 <- c(RColorBrewer::brewer.pal(8,"Set3"),Genshinpalette::Genshinpalette('HUTAO',(len_methods- len_other_methods+1))[-1])
# color_vec_7 <- c(RColorBrewer::brewer.pal(8,"Set3"),Genshinpalette::Genshinpalette('HUTAO',(len_methods- len_other_methods)))
# 
# color_vec_8 <- c(RColorBrewer::brewer.pal(8,"Set2"),Genshinpalette::Genshinpalette('HUTAO',(len_methods- len_other_methods+1))[-1])
# color_vec_9 <- c(RColorBrewer::brewer.pal(8,"Set2"),Genshinpalette::Genshinpalette('HUTAO',(len_methods- len_other_methods)))
# 
# color_vec_10 <- c(RColorBrewer::brewer.pal(8,"Set1"),Genshinpalette::Genshinpalette('HUTAO',(len_methods- len_other_methods+1))[-1])
# color_vec_11 <- c(RColorBrewer::brewer.pal(8,"Set1"),Genshinpalette::Genshinpalette('HUTAO',(len_methods- len_other_methods)))
# 
# color_vec_12 <- c(RColorBrewer::brewer.pal(8,"Accent"),Genshinpalette::Genshinpalette('HUTAO',(len_methods- len_other_methods+1))[-1])
# color_vec_13 <- c(RColorBrewer::brewer.pal(8,"Accent"),Genshinpalette::Genshinpalette('HUTAO',(len_methods- len_other_methods)))
# 
# color_vec_14 <- c(RColorBrewer::brewer.pal(8,"Dark2"),Genshinpalette::Genshinpalette('HUTAO',(len_methods- len_other_methods+1))[-1])
# color_vec_15 <- c(RColorBrewer::brewer.pal(8,"Dark2"),Genshinpalette::Genshinpalette('HUTAO',(len_methods- len_other_methods)))

#RColorBrewer::brewer.pal(10,'Paired')
#c("#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072" ,"#80B1D3" ,"#FDB462" ,"#B3DE69", "#FCCDE5" ,"#D9D9D9", "#BC80BD" ,"#CCEBC5" ,"#FFED6F")
color_setting <- c("#1B9E77", "#D95F02" ,"#7570B3", "#E7298A" ,"#66A61E" ,"#E6AB02" ,"#A6761D" ,"#666666","#66C2A5", "#FC8D62", "#8DA0CB")
#===============================================================================
#& save_df2$n1 == n_1_list[k_min] 
# data_gamma <- save_df2[save_df2$heterogeneity == heter_level_list[k_heter] & save_df2$n_0 == n_0_list[k_n0] & save_df2$inter_ratio == ratio_inter_list[k_inter] & save_df2$p_Z == p_Z_list[k_pz] & save_df2$inter_strength==inter_strength_list[k_inter_strength]& save_df2$alpha_modifier == alpha_list[k_alpha] & save_df2$eta_modifier == eta_list[k_eta], ]
# data_heterogeneity <- save_df2[save_df2$gamma == gamma_list[k_gamma] & save_df2$n_0 == n_0_list[k_n0] & save_df2$inter_ratio == ratio_inter_list[k_inter] & save_df2$p_Z == p_Z_list[k_pz] & save_df2$inter_strength==inter_strength_list[k_inter_strength]& save_df2$alpha_modifier == alpha_list[k_alpha] & save_df2$eta_modifier == eta_list[k_eta], ]
# data_n0 <- save_df2[save_df2$heterogeneity == heter_level_list[k_heter] & save_df2$gamma == gamma_list[k_gamma] & save_df2$inter_ratio == ratio_inter_list[k_inter] & save_df2$p_Z == p_Z_list[k_pz] & save_df2$inter_strength==inter_strength_list[k_inter_strength]& save_df2$alpha_modifier == alpha_list[k_alpha] & save_df2$eta_modifier == eta_list[k_eta], ]
# data_inter_ratio <- save_df2[save_df2$heterogeneity == heter_level_list[k_heter] & save_df2$gamma == gamma_list[k_gamma]& save_df2$n_0 == n_0_list[k_n0] & save_df2$p_Z == p_Z_list[k_pz] & save_df2$inter_strength==inter_strength_list[k_inter_strength]& save_df2$alpha_modifier == alpha_list[k_alpha] & save_df2$eta_modifier == eta_list[k_eta], ]
# data_pz <- save_df2[save_df2$heterogeneity == heter_level_list[k_heter] & save_df2$gamma == gamma_list[k_gamma]& save_df2$n_0 == n_0_list[k_n0] & save_df2$inter_ratio == ratio_inter_list[k_inter] & save_df2$inter_strength==inter_strength_list[k_inter_strength]& save_df2$alpha_modifier == alpha_list[k_alpha] & save_df2$eta_modifier == eta_list[k_eta], ]
# data_inter_strength <- save_df2[save_df2$heterogeneity == heter_level_list[k_heter] & save_df2$gamma == gamma_list[k_gamma2]& save_df2$n_0 == n_0_list[k_n0] & save_df2$p_Z == p_Z_list[k_pz] & save_df2$inter_ratio == ratio_inter_list[k_inter]& save_df2$alpha_modifier == alpha_list[k_alpha2] & save_df2$eta_modifier == eta_list[k_eta], ]
# data_alpha <- save_df2[save_df2$heterogeneity == heter_level_list[k_heter] & save_df2$gamma == gamma_list[k_gamma]& save_df2$n_0 == n_0_list[k_n0] & save_df2$p_Z == p_Z_list[k_pz] & save_df2$inter_ratio == ratio_inter_list[k_inter] & save_df2$inter_strength==inter_strength_list[k_inter_strength] & save_df2$eta_modifier == eta_list[k_eta], ]
# data_eta <- save_df2[save_df2$heterogeneity == heter_level_list[k_heter] & save_df2$gamma == gamma_list[k_gamma]& save_df2$n_0 == n_0_list[k_n0] & save_df2$p_Z == p_Z_list[k_pz] & save_df2$inter_ratio == ratio_inter_list[k_inter]& save_df2$inter_strength==inter_strength_list[k_inter_strength] & save_df2$alpha_modifier == alpha_list[k_alpha], ]
# pyt3 <- ggplot(data=data_n0, aes(x=n_0,y=num, colour=method, group=method)) +
#   geom_point(size=2,aes(shape=method,colour=method)) +
#   labs(x='sample size of two datasets')+
#   #geom_point(size=2,aes(x=heterogeneity,shape=method,colour=method)) +
#   scale_shape_manual(values=shapes_1) +
#   geom_line(linewidth=1.0,aes( colour=method)) +
#   #geom_line(linewidth=0.7,aes(colour=method)) +
#   scale_colour_manual(values=coloo_1)+
#   facet_grid(vars(),vars(index),scales="free_x")+
#   theme_bw() +
#   theme(plot.title = element_text(size = 15,hjust = 0.5,vjust = 0.5),
#         axis.title.y = element_blank()) +
#   ylim(-ylim,ylim)

data_b1 <- save_df2[save_df2$heterogeneity == heter_level_list[k_heter] & save_df2$gamma == gamma_list[k_gamma] & save_df2$inter_ratio == ratio_inter_list[k_inter] & save_df2$p_Z == p_Z_list[k_pz] & save_df2$inter_strength==inter_strength_list[k_inter_strength]& save_df2$alpha_modifier == alpha_list[k_alpha] & save_df2$eta_modifier == eta_list[k_eta] & save_df2$beta == beta_list[1], ]
data_b2 <- save_df2[save_df2$heterogeneity == heter_level_list[k_heter] & save_df2$gamma == gamma_list[k_gamma] & save_df2$inter_ratio == ratio_inter_list[k_inter] & save_df2$p_Z == p_Z_list[k_pz] & save_df2$inter_strength==inter_strength_list[k_inter_strength]& save_df2$alpha_modifier == alpha_list[k_alpha] & save_df2$eta_modifier == eta_list[k_eta] & save_df2$beta == beta_list[2], ]
data_b3 <- save_df2[save_df2$heterogeneity == heter_level_list[k_heter] & save_df2$gamma == gamma_list[k_gamma] & save_df2$inter_ratio == ratio_inter_list[k_inter] & save_df2$p_Z == p_Z_list[k_pz] & save_df2$inter_strength==inter_strength_list[k_inter_strength]& save_df2$alpha_modifier == alpha_list[k_alpha] & save_df2$eta_modifier == eta_list[k_eta] & save_df2$beta == beta_list[3], ]
data_b4 <- save_df2[save_df2$heterogeneity == heter_level_list[k_heter] & save_df2$gamma == gamma_list[k_gamma] & save_df2$inter_ratio == ratio_inter_list[k_inter] & save_df2$p_Z == p_Z_list[k_pz] & save_df2$inter_strength==inter_strength_list[k_inter_strength]& save_df2$alpha_modifier == alpha_list[k_alpha] & save_df2$eta_modifier == eta_list[k_eta] & save_df2$beta == beta_list[4], ]
data_b5 <- save_df2[save_df2$heterogeneity == heter_level_list[k_heter] & save_df2$gamma == gamma_list[k_gamma] & save_df2$inter_ratio == ratio_inter_list[k_inter] & save_df2$p_Z == p_Z_list[k_pz] & save_df2$inter_strength==inter_strength_list[k_inter_strength]& save_df2$alpha_modifier == alpha_list[k_alpha] & save_df2$eta_modifier == eta_list[k_eta] & save_df2$beta == beta_list[5], ]
#---------------------------------------------------------------
data_b1_bias_raw <- data_b1[data_b1$index == 'Bias',][,c('beta','n_0','method','index','num')]
data_b1_bias <- data_b1_bias_raw %>%
  arrange(n_0) %>%                      # 按 n_0 从小到大排序
  select(method, n_0, num) %>%         # 选择相关列
  pivot_wider(names_from = n_0, values_from = num)

data_b1_sd_raw <- data_b1[data_b1$index == 'SD',][,c('beta','n_0','method','index','num')]
data_b1_sd <- data_b1_bias_raw %>%
  arrange(n_0) %>%                      # 按 n_0 从小到大排序
  select(method, n_0, num) %>%         # 选择相关列
  pivot_wider(names_from = n_0, values_from = num)

data_b1_mse_raw <- data_b1[data_b1$index == 'MSE',][,c('beta','n_0','method','index','num')]
data_b1_mse <- data_b1_bias_raw %>%
  arrange(n_0) %>%                      # 按 n_0 从小到大排序
  select(method, n_0, num) %>%         # 选择相关列
  pivot_wider(names_from = n_0, values_from = num)
#------
data_b2_bias_raw <- data_b2[data_b2$index == 'Bias',][,c('beta','n_0','method','index','num')]
data_b2_bias <- data_b2_bias_raw %>%
  arrange(n_0) %>%                      # 按 n_0 从小到大排序
  select(method, n_0, num) %>%         # 选择相关列
  pivot_wider(names_from = n_0, values_from = num)

data_b2_sd_raw <- data_b2[data_b2$index == 'SD',][,c('beta','n_0','method','index','num')]
data_b2_sd <- data_b2_bias_raw %>%
  arrange(n_0) %>%                      # 按 n_0 从小到大排序
  select(method, n_0, num) %>%         # 选择相关列
  pivot_wider(names_from = n_0, values_from = num)

data_b2_mse_raw <- data_b2[data_b2$index == 'MSE',][,c('beta','n_0','method','index','num')]
data_b2_mse <- data_b2_bias_raw %>%
  arrange(n_0) %>%                      # 按 n_0 从小到大排序
  select(method, n_0, num) %>%         # 选择相关列
  pivot_wider(names_from = n_0, values_from = num)
#------
data_b3_bias_raw <- data_b3[data_b3$index == 'Bias',][,c('beta','n_0','method','index','num')]
data_b3_bias <- data_b3_bias_raw %>%
  arrange(n_0) %>%                      # 按 n_0 从小到大排序
  select(method, n_0, num) %>%         # 选择相关列
  pivot_wider(names_from = n_0, values_from = num)

data_b3_sd_raw <- data_b3[data_b3$index == 'SD',][,c('beta','n_0','method','index','num')]
data_b3_sd <- data_b3_bias_raw %>%
  arrange(n_0) %>%                      # 按 n_0 从小到大排序
  select(method, n_0, num) %>%         # 选择相关列
  pivot_wider(names_from = n_0, values_from = num)

data_b3_mse_raw <- data_b3[data_b3$index == 'MSE',][,c('beta','n_0','method','index','num')]
data_b3_mse <- data_b3_bias_raw %>%
  arrange(n_0) %>%                      # 按 n_0 从小到大排序
  select(method, n_0, num) %>%         # 选择相关列
  pivot_wider(names_from = n_0, values_from = num)
#------
data_b4_bias_raw <- data_b4[data_b4$index == 'Bias',][,c('beta','n_0','method','index','num')]
data_b4_bias <- data_b4_bias_raw %>%
  arrange(n_0) %>%                      # 按 n_0 从小到大排序
  select(method, n_0, num) %>%         # 选择相关列
  pivot_wider(names_from = n_0, values_from = num)

data_b4_sd_raw <- data_b4[data_b4$index == 'SD',][,c('beta','n_0','method','index','num')]
data_b4_sd <- data_b4_bias_raw %>%
  arrange(n_0) %>%                      # 按 n_0 从小到大排序
  select(method, n_0, num) %>%         # 选择相关列
  pivot_wider(names_from = n_0, values_from = num)

data_b4_mse_raw <- data_b4[data_b4$index == 'MSE',][,c('beta','n_0','method','index','num')]
data_b4_mse <- data_b4_bias_raw %>%
  arrange(n_0) %>%                      # 按 n_0 从小到大排序
  select(method, n_0, num) %>%         # 选择相关列
  pivot_wider(names_from = n_0, values_from = num)
#------
data_b5_bias_raw <- data_b5[data_b5$index == 'Bias',][,c('beta','n_0','method','index','num')]
data_b5_bias <- data_b5_bias_raw %>%
  arrange(n_0) %>%                      # 按 n_0 从小到大排序
  select(method, n_0, num) %>%         # 选择相关列
  pivot_wider(names_from = n_0, values_from = num)

data_b5_sd_raw <- data_b5[data_b5$index == 'SD',][,c('beta','n_0','method','index','num')]
data_b5_sd <- data_b5_bias_raw %>%
  arrange(n_0) %>%                      # 按 n_0 从小到大排序
  select(method, n_0, num) %>%         # 选择相关列
  pivot_wider(names_from = n_0, values_from = num)

data_b5_mse_raw <- data_b5[data_b5$index == 'MSE',][,c('beta','n_0','method','index','num')]
data_b5_mse <- data_b5_bias_raw %>%
  arrange(n_0) %>%                      # 按 n_0 从小到大排序
  select(method, n_0, num) %>%         # 选择相关列
  pivot_wider(names_from = n_0, values_from = num)
#---------
fwrite(data_b1_bias,file=paste(path_output_tables,'/bias_b1_',date_mark,'_',ver,'.csv',sep = ''))
fwrite(data_b1_sd,file=paste(path_output_tables,'/sd_b1_',date_mark,'_',ver,'.csv',sep = ''))
fwrite(data_b1_mse,file=paste(path_output_tables,'/mse_b1_',date_mark,'_',ver,'.csv',sep = ''))

fwrite(data_b2_bias,file=paste(path_output_tables,'/bias_b2_',date_mark,'_',ver,'.csv',sep = ''))
fwrite(data_b2_sd,file=paste(path_output_tables,'/sd_b2_',date_mark,'_',ver,'.csv',sep = ''))
fwrite(data_b2_mse,file=paste(path_output_tables,'/mse_b2_',date_mark,'_',ver,'.csv',sep = ''))

fwrite(data_b3_bias,file=paste(path_output_tables,'/bias_b3_',date_mark,'_',ver,'.csv',sep = ''))
fwrite(data_b3_sd,file=paste(path_output_tables,'/sd_b3_',date_mark,'_',ver,'.csv',sep = ''))
fwrite(data_b3_mse,file=paste(path_output_tables,'/mse_b3_',date_mark,'_',ver,'.csv',sep = ''))

fwrite(data_b4_bias,file=paste(path_output_tables,'/bias_b4_',date_mark,'_',ver,'.csv',sep = ''))
fwrite(data_b4_sd,file=paste(path_output_tables,'/sd_b4_',date_mark,'_',ver,'.csv',sep = ''))
fwrite(data_b4_mse,file=paste(path_output_tables,'/mse_b4_',date_mark,'_',ver,'.csv',sep = ''))

fwrite(data_b5_bias,file=paste(path_output_tables,'/bias_b5_',date_mark,'_',ver,'.csv',sep = ''))
fwrite(data_b5_sd,file=paste(path_output_tables,'/sd_b5_',date_mark,'_',ver,'.csv',sep = ''))
fwrite(data_b5_mse,file=paste(path_output_tables,'/mse_b5_',date_mark,'_',ver,'.csv',sep = ''))
#----
#data_pz <- data_pz[data_pz$p_Z != p_Z_fix,]
#------------------
ylim_up <- 0.1
ylim_down <- -0.2
breaks_set <- n_0_list

pyt1 <- ggplot(data = data_b1, aes(x = n_0, y = num, colour = method, group = method)) +
  geom_point(size = 2, aes(shape = method, colour = method)) +
  labs(y = paste0('β*=', as.character(beta_list[1]))) +
  scale_x_continuous(breaks = breaks_set) +
  scale_shape_manual(values = shapes_1) +
  geom_line(linewidth = 1.0, aes(colour = method)) +
  scale_colour_manual(values = color_setting) +
  facet_wrap(~ index, ncol = 3, nrow = 1, scales = "free_x") +  # 设置 ncol = 3 和 nrow = 1 实现横向排列
  theme_bw() +
  theme(
    plot.title = element_text(size = 15, hjust = 0.5, vjust = 0.5),
    #axis.title.y = element_blank()
    #,strip.text.x = element_blank()  # 隐藏分面标题
    axis.title.x = element_blank()  # 隐藏 x 轴标签
    ,axis.text.x = element_blank()    # 隐藏 x 轴坐标
    ,axis.title.y = element_text(size = zihao1)
    ,axis.text.y = element_text(size = zihao2)
  ) +
  ylim(ylim_down, ylim_up)

pyt2 <- ggplot(data = data_b2, aes(x = n_0, y = num, colour = method, group = method)) +
  geom_point(size = 2, aes(shape = method, colour = method)) +
  labs(y = paste0('β*=',as.character(beta_list[2]))) +
  scale_x_continuous(breaks = breaks_set) +
  scale_shape_manual(values = shapes_1) +
  geom_line(linewidth = 1.0, aes(colour = method)) +
  scale_colour_manual(values = color_setting) +
  facet_wrap(~ index, ncol = 3, nrow = 1, scales = "free_x") +  # 设置 ncol = 3 和 nrow = 1 实现横向排列
  theme_bw() +
  theme(plot.title = element_text(size = 15, hjust = 0.5, vjust = 0.5),
        #axis.title.y = element_blank()
        strip.text.x = element_blank()  # 隐藏分面标题
        ,axis.title.x = element_blank()  # 隐藏 x 轴标签
        ,axis.text.x = element_blank()    # 隐藏 x 轴坐标
        ,axis.title.y = element_text(size = zihao1)
        ,axis.text.y = element_text(size = zihao2)
  ) +
  ylim(ylim_down, ylim_up)
pyt3 <- ggplot(data = data_b3, aes(x = n_0, y = num, colour = method, group = method)) +
  geom_point(size = 2, aes(shape = method, colour = method)) +
  labs(y = paste0('β*=',as.character(beta_list[3]))) +
  scale_x_continuous(breaks = breaks_set) +
  scale_shape_manual(values = shapes_1) +
  geom_line(linewidth = 1.0, aes(colour = method)) +
  scale_colour_manual(values = color_setting) +
  facet_wrap(~ index, ncol = 3, nrow = 1, scales = "free_x") +  # 设置 ncol = 3 和 nrow = 1 实现横向排列
  theme_bw() +
  theme(plot.title = element_text(size = 15, hjust = 0.5, vjust = 0.5),
        #axis.title.y = element_blank()
        strip.text.x = element_blank()  # 隐藏分面标题
        ,axis.title.x = element_blank()  # 隐藏 x 轴标签
        ,axis.text.x = element_blank()    # 隐藏 x 轴坐标
        ,axis.title.y = element_text(size = zihao1)
        ,axis.text.y = element_text(size = zihao2)
  ) +
  ylim(ylim_down, ylim_up)

pyt4 <- ggplot(data = data_b4, aes(x = n_0, y = num, colour = method, group = method)) +
  geom_point(size = 2, aes(shape = method, colour = method)) +
  labs(y = paste0('β*=',as.character(beta_list[4]))
       ,x = 'sample size of two samples'
  ) +
  scale_x_continuous(breaks = breaks_set) +
  scale_shape_manual(values = shapes_1) +
  geom_line(linewidth = 1.0, aes(colour = method)) +
  scale_colour_manual(values = color_setting) +
  facet_wrap(~ index, ncol = 3, nrow = 1, scales = "free_x") +  # 设置 ncol = 3 和 nrow = 1 实现横向排列
  theme_bw() +
  theme(plot.title = element_text(size = 15, hjust = 0.5, vjust = 0.5),
        #axis.title.y = element_blank()
        strip.text.x = element_blank()  # 隐藏分面标题
        ,axis.title.x = element_blank()  # 隐藏 x 轴标签
        ,axis.text.x = element_blank()    # 隐藏 x 轴坐标
        ,axis.title.y = element_text(size = zihao1)
        ,axis.text.y = element_text(size = zihao2)
  ) +
  ylim(ylim_down, ylim_up)

pyt5 <- ggplot(data = data_b5, aes(x = n_0, y = num, colour = method, group = method)) +
  geom_point(size = 2, aes(shape = method, colour = method)) +
  labs(y = paste0('β*=',as.character(beta_list[5]))
       ,x = 'Sample size of two samples'
  ) +
  scale_x_continuous(breaks = breaks_set) +
  scale_shape_manual(values = shapes_1) +
  geom_line(linewidth = 1.0, aes(colour = method)) +
  scale_colour_manual(values = color_setting) +
  facet_wrap(~ index, ncol = 3, nrow = 1, scales = "free_x") +  # 设置 ncol = 3 和 nrow = 1 实现横向排列
  theme_bw() +
  theme(plot.title = element_text(size = 15, hjust = 0.5, vjust = 0.5),
        #axis.title.y = element_blank()
        strip.text.x = element_blank()  # 隐藏分面标题
        ,axis.title.y = element_text(size = zihao1)
        ,axis.text.y = element_text(size = zihao2)
        ,axis.title.x = element_text(size = zihao1)
  ) +
  ylim(ylim_down, ylim_up)


pyt_combine <- ggarrange(pyt1,pyt2,pyt3,pyt4,pyt5,
                         ncol=1,nrow=5,
                         #labels = c('(A)','(B)','(C)','(D)','(E)','(F)'),
                         common.legend = TRUE,
                         legend="right")

pyt_combine <- pyt_combine + theme(plot.margin = unit(c(1, 2, 1, 1), "cm")) 

ggsave(paste('MSE_',date_mark,'_',ver,'.pdf',sep = ''),plot = pyt_combine, bg="white", path = path_output_plots, height = 13.5,width = 18,dpi = 500)
ggsave(paste('MSE_',date_mark,'_',ver,'.png',sep = ''),plot = pyt_combine, bg="white", path = path_output_plots, height = 13.5,width = 18,dpi = 500)


#ggsave(paste('MSE',ver,'.png',sep = ''),plot = pyt_combine, bg="white", path = path_output_plots)
#===============================================================================绘图2
# parameter_plot <- function(data,x,y,color_setting,x_lab='',y_lab=''
#                            ,kedu_up,kedu_down,kedu_fen){
#   data_plot <- 
#   
#   ggplot(data=data_heterogeneity_beta1,aes(x = heterogeneity, y = beta_hat)) +
#     geom_boxplot(aes(fill = method),outlier.alpha=0.5) +
#     scale_fill_manual(values=color_setting)+
#     geom_hline(aes(yintercept = 1), linetype = "dashed", color = "darkred") +
#     theme_bw() +
#     scale_y_continuous(limits = c(kedu_down-gap_1, kedu_up+gap_1), 
#                        breaks = seq(kedu_down,kedu_up,kedu_fen))+
#     xlab("heterogeneity level between the two samples") +
#     ylab("Estimation") +
#     theme(axis.text.x = element_text(angle = 0,vjust = 0.5,hjust = 0.5))+ 
#     theme(strip.text.x = element_text(size = 35),
#           strip.text.y = element_text(size = 32)) +
#     theme(plot.title = element_text(size = 25,hjust = 0.5,vjust = 0.5),
#           axis.text = element_text(size = 12),
#           axis.title = element_text(size = 15),
#           legend.text = element_text(size = 18),
#           legend.title = element_text(size = 18))+  # 调整图例文本大小
#     guides(fill=guide_legend(title="Method"))
# }



#------------------------------------------
#-----------
data_gamma_beta <- save_df[save_df$heterogeneity == heter_level_list[k_heter] & save_df$n_0 == n_0_list[k_n0] & save_df$inter_ratio == ratio_inter_list[k_inter] & save_df$p_Z == p_Z_list[k_pz] & save_df$inter_strength == inter_strength_list[k_inter_strength] & save_df$alpha_modifier == alpha_list[k_alpha] & save_df$eta_modifier == eta_list[k_eta] & save_df$beta == beta_list[k_beta], ]
data_heterogeneity_beta <- save_df[save_df$gamma == gamma_list[k_gamma] & save_df$n_0 == n_0_list[k_n0] & save_df$inter_ratio == ratio_inter_list[k_inter] & save_df$p_Z == p_Z_list[k_pz] & save_df$inter_strength == inter_strength_list[k_inter_strength] & save_df$alpha_modifier == alpha_list[k_alpha] & save_df$eta_modifier == eta_list[k_eta]& save_df$beta == beta_list[k_beta], ]
data_n0_beta <- save_df[save_df$heterogeneity == heter_level_list[k_heter] & save_df$gamma == gamma_list[k_gamma] & save_df$inter_ratio == ratio_inter_list[k_inter] & save_df$p_Z == p_Z_list[k_pz] & save_df$inter_strength == inter_strength_list[k_inter_strength] & save_df$alpha_modifier == alpha_list[k_alpha] & save_df$eta_modifier == eta_list[k_eta]& save_df$beta == beta_list[k_beta], ]

condition3 <- (heter_level == heter_level_list[k_heter] & gamma_multi == gamma_list[k_gamma]  & ratio_inter == ratio_inter_list[k_inter] & p_Z==p_Z_list[k_pz] & inter_strength==inter_strength_list[k_inter_strength] & alpha_modifier == alpha_list[k_alpha] & eta_modifier == eta_list[k_eta]& beta == beta_list[k_beta]& n_02 == n_0_list2[k_n02])

data_inter_ratio_beta <- save_df[save_df$heterogeneity == heter_level_list[k_heter] & save_df$gamma == gamma_list[k_gamma]& save_df$n_0 == n_0_list[k_n0] & save_df$p_Z == p_Z_list[k_pz] & save_df$inter_strength == inter_strength_list[k_inter_strength] & save_df$alpha_modifier == alpha_list[k_alpha] & save_df$eta_modifier == eta_list[k_eta]& save_df$beta == beta_list[k_beta], ]
data_pz_beta <- save_df[save_df$heterogeneity == heter_level_list[k_heter] & save_df$gamma == gamma_list[k_gamma]& save_df$n_0 == n_0_list[k_n0] & save_df$inter_ratio == ratio_inter_list[k_inter] & save_df$inter_strength == inter_strength_list[k_inter_strength] & save_df$alpha_modifier == alpha_list[k_alpha] & save_df$eta_modifier == eta_list[k_eta]& save_df$beta == beta_list[k_beta], ]
data_inter_strength_beta <- save_df[save_df$heterogeneity == heter_level_list[k_heter] & save_df$gamma == gamma_list[k_gamma2]& save_df$n_0 == n_0_list[k_n0] & save_df$inter_ratio == ratio_inter_list[k_inter] & save_df$p_Z == p_Z_list[k_pz] & save_df$alpha_modifier == alpha_list[k_alpha2] & save_df$eta_modifier == eta_list[k_eta]& save_df$beta == beta_list[k_beta], ]
data_alpha_beta <- save_df[save_df$heterogeneity == heter_level_list[k_heter] & save_df$gamma == gamma_list[k_gamma]& save_df$n_0 == n_0_list[k_n0] & save_df$inter_ratio == ratio_inter_list[k_inter] & save_df$p_Z == p_Z_list[k_pz]& save_df$inter_strength == inter_strength_list[k_inter_strength]& save_df$eta_modifier == eta_list[k_eta]& save_df$beta == beta_list[k_beta], ]
data_eta_beta <- save_df[save_df$heterogeneity == heter_level_list[k_heter] & save_df$gamma == gamma_list[k_gamma]& save_df$n_0 == n_0_list[k_n0] & save_df$inter_ratio == ratio_inter_list[k_inter] & save_df$p_Z == p_Z_list[k_pz]& save_df$inter_strength == inter_strength_list[k_inter_strength]& save_df$alpha_modifier == alpha_list[k_alpha]& save_df$beta == beta_list[k_beta], ]
#---
data_pz_beta <- data_pz_beta[data_pz_beta$p_Z!=p_Z_fix,]

data_gamma_beta1 <- subset(data_gamma_beta, abs(beta_hat) < 8)
data_heterogeneity_beta1<- subset(data_heterogeneity_beta, abs(beta_hat) < 8)
data_n0_beta1 <- subset(data_n0_beta, abs(beta_hat) < 8)
data_inter_ratio_beta1 <- subset(data_inter_ratio_beta, abs(beta_hat) < 8)
data_pz_beta1 <- subset(data_pz_beta, abs(beta_hat) < 8)
data_inter_strength_beta1 <- subset(data_inter_strength_beta, abs(beta_hat) < 8)
data_alpha_beta1 <- subset(data_alpha_beta, abs(beta_hat) < 8)
data_eta_beta1 <- subset(data_eta_beta, abs(beta_hat) < 8)

data_heterogeneity_beta1$heterogeneity <- factor(data_heterogeneity_beta1$heterogeneity,
                                                 levels=unique(data_heterogeneity_beta1$heterogeneity),
                                                 ordered = T)
data_pz_beta1$p_Z<- factor(data_pz_beta1$p_Z,
                           levels=unique(data_pz_beta1$p_Z),
                           ordered = T)
data_gamma_beta1$gamma <- factor(data_gamma_beta1$gamma,
                                 levels = unique(data_gamma_beta1$gamma),
                                 ordered = T)
data_n0_beta1$n_0 <- factor(data_n0_beta1$n_0,
                            levels=unique(data_n0_beta1$n_0),
                            ordered = T)
data_inter_strength_beta1$inter_strength <- factor(data_inter_strength_beta1$inter_strength,
                                                   levels=unique(data_inter_strength_beta1$inter_strength),
                                                   ordered = T)
data_alpha_beta1$alpha_modifier <- factor(data_alpha_beta1$alpha_modifier,
                                          levels=unique(data_alpha_beta1$alpha_modifier),
                                          ordered = T)
data_eta_beta1$eta_modifier <- factor(data_eta_beta1$eta_modifier,
                                      levels=unique(data_eta_beta1$eta_modifier),
                                      ordered = T)

gap_1 <- 0.003

kedu_up <- 0.06
kedu_down <- -0.01
kedu_fen <- 0.01
breaks_set <- round(seq(kedu_down,kedu_up,kedu_fen),2)
#------------------------------
P1<-ggplot(data=data_heterogeneity_beta1,aes(x = heterogeneity, y = beta_hat)) +
  geom_boxplot(aes(fill = method),outlier.alpha=0.5) +
  scale_fill_manual(values=color_setting)+
  geom_hline(aes(yintercept = beta1), linetype = "dashed", color = "darkred") +
  theme_bw() +
  scale_y_continuous(limits = c(kedu_down-gap_1, kedu_up+gap_1), 
                     breaks = breaks_set)+
  xlab("Heterogeneity level between the two samples") +
  ylab("Estimation") +
  theme(axis.text.x = element_text(angle = 0,vjust = 0.5,hjust = 0.5))+ 
  theme(strip.text.x = element_text(size = 35),
        strip.text.y = element_text(size = 32)) +
  theme(plot.title = element_text(size = 25,hjust = 0.5,vjust = 0.5),
        axis.text = element_text(size = zihao2),
        axis.title = element_text(size = zihao1),
        legend.text = element_text(size = zihao3),
        legend.title = element_text(size = zihao3))+  # 调整图例文本大小
  guides(fill=guide_legend(title="Method"))

P2<-ggplot(data=data_pz_beta1,aes(x = p_Z, y = beta_hat)) +
  geom_boxplot(aes(fill = method),outlier.alpha=0.5) +
  scale_fill_manual(values=color_setting)+
  geom_hline(aes(yintercept = beta1), linetype = "dashed", color = "darkred") +
  theme_bw() +
  scale_y_continuous(limits = c(kedu_down-gap_1, kedu_up+gap_1), 
                     breaks = breaks_set)+
  xlab("Number of instrumental variables") +
  ylab("") +
  theme(axis.title.y = element_blank(),  # 去掉 y 轴标题
        axis.text.y = element_blank(),   # 去掉 y 轴文本
        axis.ticks.y = element_blank()) + # 去掉 y 轴刻度
  theme(axis.text.x = element_text(angle = 0,vjust = 0.5,hjust = 0.5))+ 
  theme(strip.text.x = element_text(size = 35),
        strip.text.y = element_text(size = 32)) +
  theme(plot.title = element_text(size = 25,hjust = 0.5,vjust = 0.5),
        axis.text = element_text(size = zihao2),
        axis.title = element_text(size = zihao1),
        legend.text = element_text(size = zihao3),
        legend.title = element_text(size = zihao3))+  # 调整图例文本大小
  guides(fill=guide_legend(title="Method"))

P2a<-ggplot(data=data_pz_beta1,aes(x = p_Z, y = beta_hat)) +
  geom_boxplot(aes(fill = method),outlier.alpha=0.5) +
  scale_fill_manual(values=color_setting)+
  geom_hline(aes(yintercept = beta1), linetype = "dashed", color = "darkred") +
  theme_bw() +
  scale_y_continuous()+
  xlab("Number of instrumental variables") +
  ylab("Estimation") +
  # theme(axis.title.y = element_blank(),  # 去掉 y 轴标题
  #       axis.text.y = element_blank(),   # 去掉 y 轴文本
  #       axis.ticks.y = element_blank()) + # 去掉 y 轴刻度
  theme(axis.text.x = element_text(angle = 0,vjust = 0.5,hjust = 0.5))+ 
  theme(strip.text.x = element_text(size = 35),
        strip.text.y = element_text(size = 32)) +
  theme(plot.title = element_text(size = 25,hjust = 0.5,vjust = 0.5),
        axis.text = element_text(size = zihao2),
        axis.title = element_text(size = zihao1),
        legend.text = element_text(size = zihao3),
        legend.title = element_text(size = zihao3))+  # 调整图例文本大小
  guides(fill=guide_legend(title="Method"))

P3<-ggplot(data=data_gamma_beta1,aes(x = gamma, y = beta_hat)) +
  geom_boxplot(aes(fill = method),outlier.alpha=0.5) +
  scale_fill_manual(values=color_setting)+
  geom_hline(aes(yintercept = beta1), linetype = "dashed", color = "darkred") +
  theme_bw() +
  scale_y_continuous(limits = c(kedu_down-gap_1, kedu_up+gap_1), 
                     breaks = breaks_set)+
  xlab("Relative magnitude of horizontal pleiotropy effect") +
  ylab("Estimation") +
  #ylab("") +
  # theme(axis.title.y = element_blank(),  # 去掉 y 轴标题
  #       axis.text.y = element_blank(),   # 去掉 y 轴文本
  #       axis.ticks.y = element_blank()) + # 去掉 y 轴刻度
  theme(axis.text.x = element_text(angle = 0,vjust = 0.5,hjust = 0.5))+ 
  theme(strip.text.x = element_text(size = 35),
        strip.text.y = element_text(size = 32)) +
  theme(plot.title = element_text(size = 25,hjust = 0.5,vjust = 0.5),
        axis.text = element_text(size = zihao2),
        axis.title = element_text(size = zihao1),
        legend.text = element_text(size = zihao3),
        legend.title = element_text(size = zihao3))+  # 调整图例文本大小
  guides(fill=guide_legend(title="Method"))

P4<-ggplot(data=data_n0_beta1,aes(x = n_0, y = beta_hat)) +
  geom_boxplot(aes(fill = method),outlier.alpha=0.5) +
  scale_fill_manual(values=color_setting)+
  geom_hline(aes(yintercept = beta1), linetype = "dashed", color = "darkred") +
  theme_bw() +
  scale_y_continuous(limits = c(kedu_down-gap_1, kedu_up+gap_1), 
                     breaks = breaks_set)+
  xlab("Sample size of two samples") +
  #ylab("") +
  ylab("Estimation") +
  # theme(axis.title.y = element_blank(),  # 去掉 y 轴标题
  #       axis.text.y = element_blank(),   # 去掉 y 轴文本
  #       axis.ticks.y = element_blank()) + # 去掉 y 轴刻度
  theme(axis.text.x = element_text(angle = 0,vjust = 0.5,hjust = 0.5))+ # 设置 x 轴文本的角度和位置，angle=0 表示水平显示，vjust 和 hjust 分别表示垂直和水平对齐方式
  theme(strip.text.x = element_text(size = 35),
        strip.text.y = element_text(size = 32)) +  # 设置 facet 网格中 x\y 方向标签的字号
  theme(plot.title = element_text(size = 25,hjust = 0.5,vjust = 0.5), # 设置图标题的字号、水平对齐（hjust）和垂直对齐（vjust）
        axis.text = element_text(size = zihao2), # 设置 x 轴和 y 轴刻度文本的字号
        axis.title = element_text(size = zihao1), # 设置 x 轴和 y 轴标题的字号
        legend.text = element_text(size = zihao3), # 设置图例项文本的字号
        legend.title = element_text(size = zihao3))+  # 调整图例文本大小
  guides(fill=guide_legend(title="Method"))

P5<-ggplot(data=data_inter_strength_beta1,aes(x = inter_strength, y = beta_hat)) +
  geom_boxplot(aes(fill = method),outlier.alpha=0.5) +
  scale_fill_manual(values=color_setting)+
  geom_hline(aes(yintercept = beta1), linetype = "dashed", color = "darkred") +
  theme_bw() +
  scale_y_continuous(limits = c(kedu_down-gap_1, kedu_up+gap_1), 
                     breaks = breaks_set)+
  xlab("Relative magnitude of interaction effects") +
  #ylab("") +
  ylab("Estimation") +
  theme(axis.title.y = element_blank(),  # 去掉 y 轴标题
        axis.text.y = element_blank(),   # 去掉 y 轴文本
        axis.ticks.y = element_blank()) + # 去掉 y 轴刻度
  theme(axis.text.x = element_text(angle = 0,vjust = 0.5,hjust = 0.5))+ 
  theme(strip.text.x = element_text(size = 35),
        strip.text.y = element_text(size = 32)) +
  theme(plot.title = element_text(size = 25,hjust = 0.5,vjust = 0.5),
        axis.text = element_text(size = zihao2),
        axis.title = element_text(size = zihao1),
        legend.text = element_text(size = zihao3),
        legend.title = element_text(size = zihao3))+  # 调整图例文本大小
  guides(fill=guide_legend(title="Method"))

kedu_down1_5 <- floor(summary(data_inter_strength_beta1$beta_hat)[1]*10)/10  
kedu_up1_5 <- ceiling(summary(data_inter_strength_beta1$beta_hat)[6]*10)/10  
kedu_down1_5 <- 0
kedu_up1_5 <- 0.1
breaks_set1_5 <- round(seq(kedu_down1_5,kedu_up1_5,kedu_fen),2)
P5a<-ggplot(data=data_inter_strength_beta1,aes(x = inter_strength, y = beta_hat)) +
  geom_boxplot(aes(fill = method),outlier.alpha=0.5) +
  scale_fill_manual(values=color_setting)+
  geom_hline(aes(yintercept = beta1), linetype = "dashed", color = "darkred") +
  theme_bw() +
  scale_y_continuous(limits = c(kedu_down1_5, kedu_up1_5), 
                     breaks = breaks_set1_5)+
  xlab("Relative magnitude of interaction effects") +
  #ylab("") +
  ylab("Estimation") +
  # theme(axis.title.y = element_blank(),  # 去掉 y 轴标题
  #       axis.text.y = element_blank(),   # 去掉 y 轴文本
  #       axis.ticks.y = element_blank()) + # 去掉 y 轴刻度
  theme(axis.text.x = element_text(angle = 0,vjust = 0.5,hjust = 0.5))+ 
  theme(strip.text.x = element_text(size = 35),
        strip.text.y = element_text(size = 32)) +
  theme(plot.title = element_text(size = 25,hjust = 0.5,vjust = 0.5),
        axis.text = element_text(size = zihao2),
        axis.title = element_text(size = zihao1),
        legend.text = element_text(size = zihao3),
        legend.title = element_text(size = zihao3))+  # 调整图例文本大小
  guides(fill=guide_legend(title="Method"))


P6<-ggplot(data=data_alpha_beta1,aes(x = alpha_modifier, y = beta_hat)) +
  geom_boxplot(aes(fill = method),outlier.alpha=0.5) +
  scale_fill_manual(values=color_setting)+
  geom_hline(aes(yintercept = beta1), linetype = "dashed", color = "darkred") +
  theme_bw() +
  scale_y_continuous(limits = c(kedu_down-gap_1, kedu_up+gap_1), 
                     breaks = breaks_set)+
  xlab("Relative magnitude of the effect of Z on X") +
  #ylab("Estimation") +
  ylab("") +
  theme(axis.title.y = element_blank(),  # 去掉 y 轴标题
        axis.text.y = element_blank(),   # 去掉 y 轴文本
        axis.ticks.y = element_blank()) + # 去掉 y 轴刻度
  theme(axis.text.x = element_text(angle = 0,vjust = 0.5,hjust = 0.5))+ 
  theme(strip.text.x = element_text(size = 35),
        strip.text.y = element_text(size = 32)) +
  theme(plot.title = element_text(size = 25,hjust = 0.5,vjust = 0.5),
        axis.text = element_text(size = zihao2),
        axis.title = element_text(size = zihao1),
        legend.text = element_text(size = zihao3),
        legend.title = element_text(size = zihao3))+  # 调整图例文本大小
  guides(fill=guide_legend(title="Method"))

P6a<-ggplot(data=data_alpha_beta1,aes(x = alpha_modifier, y = beta_hat)) +
  geom_boxplot(aes(fill = method),outlier.alpha=0.5) +
  scale_fill_manual(values=color_setting)+
  geom_hline(aes(yintercept = beta1), linetype = "dashed", color = "darkred") +
  theme_bw() +
  scale_y_continuous(limits = c(kedu_down-gap_1, kedu_up+gap_1), 
                     breaks = breaks_set)+
  xlab("Relative magnitude of the effect of Z on X") +
  #ylab("Estimation") +
  ylab("") +
  # theme(axis.title.y = element_blank(),  # 去掉 y 轴标题
  #       axis.text.y = element_blank(),   # 去掉 y 轴文本
  #       axis.ticks.y = element_blank()) + # 去掉 y 轴刻度
  theme(axis.text.x = element_text(angle = 0,vjust = 0.5,hjust = 0.5))+ 
  theme(strip.text.x = element_text(size = 35),
        strip.text.y = element_text(size = 32)) +
  theme(plot.title = element_text(size = 25,hjust = 0.5,vjust = 0.5),
        axis.text = element_text(size = zihao2),
        axis.title = element_text(size = zihao1),
        legend.text = element_text(size = zihao3),
        legend.title = element_text(size = zihao3))+  # 调整图例文本大小
  guides(fill=guide_legend(title="Method"))

P7<-ggplot(data=data_eta_beta1,aes(x = eta_modifier, y = beta_hat)) +
  geom_boxplot(aes(fill = method),outlier.alpha=0.5) +
  scale_fill_manual(values=color_setting)+
  geom_hline(aes(yintercept = beta1), linetype = "dashed", color = "darkred") +
  theme_bw() +
  scale_y_continuous(limits = c(kedu_down-gap_1, kedu_up+gap_1), 
                     breaks = breaks_set)+
  xlab("Relative magnitude of the effect of C on X and Y") +
  #ylab("Estimation") +
  ylab("") +
  theme(axis.text.x = element_text(angle = 0,vjust = 0.5,hjust = 0.5))+ 
  theme(strip.text.x = element_text(size = 35),
        strip.text.y = element_text(size = 32)) +
  theme(plot.title = element_text(size = 25,hjust = 0.5,vjust = 0.5),
        axis.text = element_text(size = zihao2),
        axis.title = element_text(size = zihao1),
        legend.text = element_text(size = zihao3),
        legend.title = element_text(size = zihao3))+  # 调整图例文本大小
  guides(fill=guide_legend(title="Method"))


pyt_combine2_1 <- ggarrange(P1
                            ,P5,
                            ncol=2,nrow=1,
                            labels = c('(A)','(B)'),
                            # labels.x = 1,
                            hjust= -c(2.4,1),
                            vjust = 2,
                            widths = c(1,1),
                            common.legend = TRUE,
                            legend="right")

pyt_combine2_2 <- ggarrange(P4
                            ,P2,
                            ncol=2,nrow=1,
                            labels = c('(A)','(B)'),
                            # labels.x = 1,
                            hjust= -c(2.4,1),
                            vjust = 2,
                            widths = c(1,1),
                            common.legend = TRUE,
                            legend="right")

pyt_combine2_3 <- ggarrange(P3
                            ,P6,
                            ncol=2,nrow=1,
                            labels = c('(A)','(B)'),
                            # labels.x = 1,
                            hjust= -c(2.4,1),
                            vjust = 2,
                            widths = c(1,1),
                            common.legend = TRUE,
                            legend="right")

pyt_combine2a <- ggarrange(P1,
                           P5a,
                           P4,
                           ncol=1,nrow=3,
                           labels = c('(A)','(B)','(C)'),
                           # labels.x = 1,
                           hjust= -2.4,
                           vjust = 2,
                           #widths = c(1,1),
                           common.legend = TRUE,
                           legend="right")
pyt_combine2b <- ggarrange(P2a,
                           P3,
                           P6a,
                           ncol=1,nrow=3,
                           labels = c('(A)','(B)','(C)'),
                           # labels.x = 1,
                           hjust= -2.4,
                           vjust = 2,
                           #widths = c(1,1),
                           common.legend = TRUE,
                           legend="right")

pyt_combine2i <- ggarrange(P1,
                           P5a,
                           ncol=1,nrow=2,
                           labels = c('(A)','(B)','(C)'),
                           # labels.x = 1,
                           hjust= -2.4,
                           vjust = 2,
                           #widths = c(1,1),
                           common.legend = TRUE,
                           legend="right")
pyt_combine2j <- ggarrange(P4,
                           P2a,
                           ncol=1,nrow=2,
                           labels = c('(A)','(B)','(C)'),
                           # labels.x = 1,
                           hjust= -2.4,
                           vjust = 2,
                           #widths = c(1,1),
                           common.legend = TRUE,
                           legend="right")
pyt_combine2k <- ggarrange(P3,
                           P6a,
                           ncol=1,nrow=2,
                           labels = c('(A)','(B)','(C)'),
                           # labels.x = 1,
                           hjust= -2.4,
                           vjust = 2,
                           #widths = c(1,1),
                           common.legend = TRUE,
                           legend="right")

pyt_combine2 <- ggarrange(P1,P5
                          ,P4,P2
                          ,P3,P6,
                          ncol=2,nrow=3,
                          labels = c('(A)','(B)','(C)','(D)','(E)','(F)'),
                          # labels.x = 1,
                          hjust= -c(2.4,1,2.4,1,2.4,1),
                          vjust = 2,
                          widths = c(1,1),
                          common.legend = TRUE,
                          legend="right")

pyt_combine2

# pyt_combine2 <- pyt_combine2 + theme(plot.margin = unit(c(1, 2, 1, 1), "cm")) 
# pyt_combine2_1 <- pyt_combine2_1 + theme(plot.margin = unit(c(1, 2, 1, 1), "cm")) 
# pyt_combine2_2 <- pyt_combine2_2 + theme(plot.margin = unit(c(1, 2, 1, 1), "cm")) 
# pyt_combine2_3 <- pyt_combine2_3 + theme(plot.margin = unit(c(1, 2, 1, 1), "cm")) 
# 
# pyt_combine2_a <- pyt_combine2_a + theme(plot.margin = unit(c(1, 2, 1, 1), "cm")) 
# pyt_combine2_b <- pyt_combine2_b + theme(plot.margin = unit(c(1, 2, 1, 1), "cm")) 
# pyt_combine2_i <- pyt_combine2_i + theme(plot.margin = unit(c(1, 2, 1, 1), "cm")) 
# pyt_combine2_j <- pyt_combine2_j + theme(plot.margin = unit(c(1, 2, 1, 1), "cm")) 
# pyt_combine2_k <- pyt_combine2_k + theme(plot.margin = unit(c(1, 2, 1, 1), "cm")) 

pyt_combine2 <- pyt_combine2 + theme(plot.margin = unit(c(0, 1, 0, 0), "cm")) 
pyt_combine2_1 <- pyt_combine2_1 + theme(plot.margin = unit(c(0, 1, 0, 0), "cm")) 
pyt_combine2_2 <- pyt_combine2_2 + theme(plot.margin = unit(c(0, 1, 0, 0), "cm")) 
pyt_combine2_3 <- pyt_combine2_3 + theme(plot.margin = unit(c(0, 1, 0, 0), "cm")) 

pyt_combine2a <- pyt_combine2a + theme(plot.margin = unit(c(0, 1, 0, 0), "cm")) 
pyt_combine2b <- pyt_combine2b + theme(plot.margin = unit(c(0, 1, 0, 0), "cm")) 
pyt_combine2i <- pyt_combine2i + theme(plot.margin = unit(c(0, 1, 0, 0), "cm")) 
pyt_combine2j <- pyt_combine2j + theme(plot.margin = unit(c(0, 1, 0, 0), "cm")) 
pyt_combine2k <- pyt_combine2k + theme(plot.margin = unit(c(0, 1, 0, 0), "cm")) 


ggsave(paste('parameters_all_',date_mark,'_',ver,'.png',sep = ''),plot = pyt_combine2, bg="white", path = path_output_plots,width = 18,height = 13.5,dpi = 500)
ggsave(paste('parameters_all_',date_mark,'_',ver,'.pdf',sep = ''),plot = pyt_combine2, bg="white", path = path_output_plots,width = 18,height = 13.5,dpi = 500)

ggsave(paste('parameters1_',date_mark,'_',ver,'.png',sep = ''),plot = pyt_combine2_1, bg="white", path = path_output_plots,width = 18,height = 8.1,dpi = 500)
ggsave(paste('parameters1_',date_mark,'_',ver,'.pdf',sep = ''),plot = pyt_combine2_1, bg="white", path = path_output_plots,width = 18,height = 8.1,dpi = 500)

ggsave(paste('parameters2_',date_mark,'_',ver,'.png',sep = ''),plot = pyt_combine2_2, bg="white", path = path_output_plots,width = 18,height = 8.1,dpi = 500)
ggsave(paste('parameters2_',date_mark,'_',ver,'.pdf',sep = ''),plot = pyt_combine2_2, bg="white", path = path_output_plots,width = 18,height = 8.1,dpi = 500)

ggsave(paste('parameters3_',date_mark,'_',ver,'.png',sep = ''),plot = pyt_combine2_3, bg="white", path = path_output_plots,width = 18,height = 8.1,dpi = 500)
ggsave(paste('parameters3_',date_mark,'_',ver,'.pdf',sep = ''),plot = pyt_combine2_3, bg="white", path = path_output_plots,width = 18,height = 8.1,dpi = 500)

ggsave(paste('parameters_a_',date_mark,'_',ver,'.png',sep = ''),plot = pyt_combine2a, bg="white", path = path_output_plots,width = 18,height = 18,dpi = 500)
ggsave(paste('parameters_a_',date_mark,'_',ver,'.pdf',sep = ''),plot = pyt_combine2a, bg="white", path = path_output_plots,width = 18,height = 18,dpi = 500)

ggsave(paste('parameters_b_',date_mark,'_',ver,'.png',sep = ''),plot = pyt_combine2b, bg="white", path = path_output_plots,width = 18,height = 18,dpi = 500)
ggsave(paste('parameters_b_',date_mark,'_',ver,'.pdf',sep = ''),plot = pyt_combine2b, bg="white", path = path_output_plots,width = 18,height = 18,dpi = 500)

ggsave(paste('parameters_i_',date_mark,'_',ver,'.png',sep = ''),plot = pyt_combine2i, bg="white", path = path_output_plots,width = 18,height = 13.5,dpi = 500)
ggsave(paste('parameters_i_',date_mark,'_',ver,'.pdf',sep = ''),plot = pyt_combine2i, bg="white", path = path_output_plots,width = 18,height = 13.5,dpi = 500)

ggsave(paste('parameters_j_',date_mark,'_',ver,'.png',sep = ''),plot = pyt_combine2j, bg="white", path = path_output_plots,width = 18,height = 13.5,dpi = 500)
ggsave(paste('parameters_j_',date_mark,'_',ver,'.pdf',sep = ''),plot = pyt_combine2j, bg="white", path = path_output_plots,width = 18,height = 13.5,dpi = 500)

ggsave(paste('parameters_k_',date_mark,'_',ver,'.png',sep = ''),plot = pyt_combine2k, bg="white", path = path_output_plots,width = 18,height = 13.5,dpi = 500)
ggsave(paste('parameters_k_',date_mark,'_',ver,'.pdf',sep = ''),plot = pyt_combine2k, bg="white", path = path_output_plots,width = 18,height = 13.5,dpi = 500)

###
ggsave(paste('heter_',date_mark,'_',ver,'.png',sep = ''),plot = P1, bg="white", path = path_output_plots,width = 18,height = 13.5,dpi = 500)
ggsave(paste('heter_',date_mark,'_',ver,'.pdf',sep = ''),plot = P1, bg="white", path = path_output_plots,width = 18,height = 13.5,dpi = 500)

ggsave(paste('pz_',date_mark,'_',ver,'.png',sep = ''),plot = P2a, bg="white", path = path_output_plots,width = 18,height = 13.5,dpi = 500)
ggsave(paste('pz_',date_mark,'_',ver,'.pdf',sep = ''),plot = P2a, bg="white", path = path_output_plots,width = 18,height = 13.5,dpi = 500)

ggsave(paste('gamma_',date_mark,'_',ver,'.png',sep = ''),plot = P3, bg="white", path = path_output_plots,width = 18,height = 13.5,dpi = 500)
ggsave(paste('gamma_',date_mark,'_',ver,'.pdf',sep = ''),plot = P3, bg="white", path = path_output_plots,width = 18,height = 13.5,dpi = 500)

ggsave(paste('n0_',date_mark,'_',ver,'.png',sep = ''),plot = P4, bg="white", path = path_output_plots,width = 18,height = 13.5,dpi = 500)
ggsave(paste('n0_',date_mark,'_',ver,'.pdf',sep = ''),plot = P4, bg="white", path = path_output_plots,width = 18,height = 13.5,dpi = 500)

ggsave(paste('pi_',date_mark,'_',ver,'.png',sep = ''),plot = P5a, bg="white", path = path_output_plots,width = 18,height = 13.5,dpi = 500)
ggsave(paste('pi_',date_mark,'_',ver,'.pdf',sep = ''),plot = P5a, bg="white", path = path_output_plots,width = 18,height = 13.5,dpi = 500)

ggsave(paste('alpha_',date_mark,'_',ver,'.png',sep = ''),plot = P6a, bg="white", path = path_output_plots,width = 18,height = 13.5,dpi = 500)
ggsave(paste('alpha_',date_mark,'_',ver,'.pdf',sep = ''),plot = P6a, bg="white", path = path_output_plots,width = 18,height = 13.5,dpi = 500)

#--------------
# data_inter_strength_beta1<-data_inter_strength_beta1[data_inter_strength_beta1$inter_strength != 0,]
# 
# P51<-ggplot(data=data_inter_strength_beta1,aes(x = inter_strength, y = beta_hat)) +
#   geom_boxplot(aes(fill = method),outlier.alpha=0.5) +
#   scale_fill_manual(values=color_setting)+
#   geom_hline(aes(yintercept = 1), linetype = "dashed", color = "darkred") +
#   theme_bw() +
#   scale_y_continuous(limits = c(kedu_down-gap_1, kedu_up+gap_1), 
#                      breaks = seq(kedu_down,kedu_up,kedu_fen))+
#   xlab("Relative magnitude of interaction effects") +
#   #ylab("Estimation") +
#   ylab("") +
#   theme(axis.text.x = element_text(angle = 0,vjust = 0.5,hjust = 0.5))+ 
#   theme(strip.text.x = element_text(size = 35),
#         strip.text.y = element_text(size = 32)) +
#   theme(plot.title = element_text(size = 25,hjust = 0.5,vjust = 0.5),
#         axis.text = element_text(size = zihao2),
#         axis.title = element_text(size = zihao1),
#         legend.text = element_text(size = zihao3),
#         legend.title = element_text(size = zihao3))+  # 调整图例文本大小
#   guides(fill=guide_legend(title="Method"))
# pyt_combine2 <- ggarrange(P1,P3,P51,P6,P7,
#                           ncol=3,nrow=2,
#                           #labels = c('(A)','(B)','(C)','(D)','(E)','(F)'),
#                           common.legend = TRUE,
#                           legend="right")
# print(pyt_combine2)
#--------------
# ggsave(paste('gamma_2_',date_mark,'_',ver,'.png',sep = ''),plot = P3, bg="white", path = path_output_plots)
# ggsave(paste('heterogeneity_2_',date_mark,'_',ver,'.png',sep = ''),plot = P1, bg="white", path = path_output_plots)
# ggsave(paste('n_0_2_',date_mark,'_',ver,'.png',sep = ''),plot = P4, bg="white", path = path_output_plots)
# #ggsave(paste('inter_ratio_2_',date_mark,'_',ver,'.png',sep = ''),plot = pyt5, bg="white", path = path_output_plots)
# ggsave(paste('inter_strength_2_',date_mark,'_',ver,'.png',sep = ''),plot = P5, bg="white", path = path_output_plots)
# ggsave(paste('pZ_2_',date_mark,'_',ver,'.png',sep = ''),plot = P2, bg="white", path = path_output_plots)
# ggsave(paste('alpha_2_',date_mark,'_',ver,'.png',sep = ''),plot = P6, bg="white", path = path_output_plots)
# ggsave(paste('eta_2_',date_mark,'_',ver,'.png',sep = ''),plot = P7, bg="white", path = path_output_plots)
# ggsave(paste('parameters_2_',date_mark,'_',ver,'.png',sep = ''),plot = pyt_combine2, bg="white", path = path_output_plots,width = 20,height = 15,dpi = 600)
# ggsave(paste('parameters_2_',date_mark,'_',ver,'.pdf',sep = ''),plot = pyt_combine2, bg="white", path = path_output_plots,width = 20,height = 15,dpi = 600)
#===============================================================================绘图3
hist.plot <- function(power,legend_label,x_label,y_label,coloo){
  po <- data.frame(Method=factor(legend_label,
                                 levels =legend_label),
                   pow=power)
  gt <- ggplot(po,aes(x=Method,y=pow,fill=Method))+
    geom_bar(stat="identity",width=0.5)+
    theme_bw()+
    theme(panel.background = element_blank()) +
    scale_fill_manual(values=coloo)+
    scale_y_continuous(limits = c(0, 1), 
                       breaks = c(0,0.2,0.4,0.6,0.8,1))+
    #theme(panel.grid =element_blank()) + 
    theme(axis.text= element_text(size = 18),
          axis.title.x=element_text(size=19),
          axis.title.y=element_text(size=19),
          plot.title=element_text(size=23)) +
    labs(x = x_label, y = y_label) +
    theme(legend.title = element_text(size = 15),
          legend.text=element_text(size=15))+  
    geom_text(mapping = aes(label = pow),vjust = -0.5)
  return(gt)
}

#注意：method一栏需要是排好序的factor
power_cover_plot <- function(data_test,color_setting,x_exist=TRUE,y_exist=TRUE,x_label,zihao1,zihao2,zihao3,zihao4,lable_list=c("A","B")){  
  #x轴y轴都有
  if(x_exist & y_exist){                     
    power_plot_v2 <- ggplot(data_test,aes(x=method,y=power,fill=method))+
      geom_bar(stat="identity",width=0.5)+
      theme_bw()+
      theme(panel.background = element_blank()) +
      scale_fill_manual(values=color_setting)+
      scale_y_continuous(limits = c(0, 1), 
                         breaks = c(0.95))+
      #theme(panel.grid =element_blank()) + 
      theme(axis.text= element_text(size = zihao2),
            axis.title.x = element_blank(),  # 移除 x 轴标题
            axis.text.x = element_blank(),   # 移除 x 轴刻度
            axis.ticks.x = element_blank(),  # 移除 x 轴刻度线
            axis.title.y=element_text(size=zihao1),
            plot.title=element_text(size=23)
            #axis.text.y = element_text(color = "darkred")# 修改所有 y 轴刻度的颜色
      ) +
      labs(y = 'Power') +
      theme(legend.title = element_text(size = zihao3),
            legend.text=element_text(size=zihao3))+  
      geom_text(mapping = aes(label = power),vjust = -0.5,size=zihao4)
    cover_plot_v2 <- ggplot(data_test,aes(x=method,y=cover_ratio,fill=method))+
      geom_bar(stat="identity",width=0.5)+
      theme_bw()+
      theme(panel.background = element_blank()) +
      scale_fill_manual(values=color_setting)+
      scale_y_continuous(limits = c(0, 1), 
                         breaks = c(0.95))+
      #theme(panel.grid =element_blank()) + 
      theme(axis.text= element_text(size = zihao2),
            #axis.title.x = element_blank(),  # 移除 x 轴标题
            axis.title.x=element_text(size=zihao1),
            axis.text.x = element_blank(),   # 移除 x 轴刻度
            axis.ticks.x = element_blank(),  # 移除 x 轴刻度线
            axis.title.y=element_text(size=zihao1),
            plot.title=element_text(size=23)) +
      labs(y = 'Cover Ratio',
           x = x_label) +
      theme(legend.title = element_text(size = zihao3),
            legend.text=element_text(size=zihao3))+  
      geom_text(mapping = aes(label = cover_ratio),vjust = -0.5,size=zihao4)
    
    # pyt_combine3 <- ggarrange(power_plot_v2,cover_plot_v2,
    #                           ncol=1,nrow=2,
    #                           #labels = c('(A)','(B)','(C)','(D)','(E)','(F)'),
    #                           common.legend = TRUE,
    #                           legend="none")
    #只有x轴没有y轴
  }else if(x_exist & (!y_exist) ){  
    power_plot_v2 <- ggplot(data_test,aes(x=method,y=power,fill=method))+
      geom_bar(stat="identity",width=0.5)+
      theme_bw()+
      theme(panel.background = element_blank()) +
      scale_fill_manual(values=color_setting)+
      scale_y_continuous(limits = c(0, 1), 
                         breaks = c(0.95))+
      #theme(panel.grid =element_blank()) + 
      theme(axis.text= element_text(size = zihao2),
            axis.title.x = element_blank(),  # 移除 x 轴标题
            axis.text.x = element_blank(),   # 移除 x 轴刻度
            axis.ticks.x = element_blank(),  # 移除 x 轴刻度线
            
            axis.title.y = element_blank(),  # 移除 y 轴标题
            axis.text.y = element_blank(),   # 移除 y 轴刻度
            axis.ticks.y = element_blank(),  # 移除 y 轴刻度线
            
            #axis.title.y=element_text(size=zihao1),
            plot.title=element_text(size=23)) +
      labs(y = '') +
      theme(legend.title = element_text(size = zihao3),
            legend.text=element_text(size=zihao3))+  
      geom_text(mapping = aes(label = power),vjust = -0.5,size=zihao4)
    
    cover_plot_v2 <- ggplot(data_test,aes(x=method,y=cover_ratio,fill=method))+
      geom_bar(stat="identity",width=0.5)+
      theme_bw()+
      theme(panel.background = element_blank()) +
      scale_fill_manual(values=color_setting)+
      scale_y_continuous(limits = c(0, 1), 
                         breaks = c(0.95))+
      #theme(panel.grid =element_blank()) + 
      theme(axis.text= element_text(size = zihao2),
            #axis.title.x = element_blank(),  # 移除 x 轴标题
            axis.title.x=element_text(size=zihao1),
            axis.text.x = element_blank(),   # 移除 x 轴刻度
            axis.ticks.x = element_blank(),  # 移除 x 轴刻度线
            
            axis.title.y = element_blank(),  # 移除 y 轴标题
            axis.text.y = element_blank(),   # 移除 y 轴刻度
            axis.ticks.y = element_blank(),  # 移除 y 轴刻度线
            
            #axis.title.y=element_text(size=zihao1),
            plot.title=element_text(size=23)) +
      labs(y = '',
           x = x_label) +
      theme(legend.title = element_text(size = zihao3),
            legend.text=element_text(size=zihao3))+  
      geom_text(mapping = aes(label = cover_ratio),vjust = -0.5,size=zihao4)
    
    # pyt_combine3 <- ggarrange(power_plot_v2,cover_plot_v2,
    #                           ncol=1,nrow=2,
    #                           #labels = c('(A)','(B)','(C)','(D)','(E)','(F)'),
    #                           common.legend = TRUE,
    #                           legend="none")
    #没有x轴只有y轴
  }else if( (!x_exist) & y_exist){  
    power_plot_v2 <- ggplot(data_test,aes(x=method,y=power,fill=method))+
      geom_bar(stat="identity",width=0.5)+
      theme_bw()+
      theme(panel.background = element_blank()) +
      scale_fill_manual(values=color_setting)+
      scale_y_continuous(limits = c(0, 1), 
                         breaks = c(0.95))+
      #theme(panel.grid =element_blank()) + 
      theme(axis.text= element_text(size = zihao2),
            axis.title.x = element_blank(),  # 移除 x 轴标题
            axis.text.x = element_blank(),   # 移除 x 轴刻度
            axis.ticks.x = element_blank(),  # 移除 x 轴刻度线
            
            axis.title.y=element_text(size=zihao1),
            plot.title=element_text(size=23)) +
      labs(y = 'Power') +
      theme(legend.title = element_text(size = zihao3),
            legend.text=element_text(size=zihao3))+  
      geom_text(mapping = aes(label = power),vjust = -0.5,size=zihao4)
    
    cover_plot_v2 <- ggplot(data_test,aes(x=method,y=cover_ratio,fill=method))+
      geom_bar(stat="identity",width=0.5)+
      theme_bw()+
      theme(panel.background = element_blank()) +
      scale_fill_manual(values=color_setting)+
      scale_y_continuous(limits = c(0, 1), 
                         breaks = c(0.95))+
      #theme(panel.grid =element_blank()) + 
      theme(axis.text= element_text(size = zihao2),
            #axis.title.x = element_blank(),  # 移除 x 轴标题
            axis.title.x=element_text(size=zihao1),
            axis.text.x = element_blank(),   # 移除 x 轴刻度
            axis.ticks.x = element_blank(),  # 移除 x 轴刻度线
            
            axis.title.y=element_text(size=zihao1),
            plot.title=element_text(size=23)) +
      labs(y = 'Cover Ratio',
           x = x_label) +
      theme(legend.title = element_text(size = zihao3),
            legend.text=element_text(size=zihao3))+  
      geom_text(mapping = aes(label = cover_ratio),vjust = -0.5,size=zihao4)
    
    # pyt_combine3 <- ggarrange(power_plot_v2,cover_plot_v2,
    #                           ncol=1,nrow=2,
    #                           #labels = c('(A)','(B)','(C)','(D)','(E)','(F)'),
    #                           common.legend = TRUE,
    #                           legend="none")
    # x轴没有y轴也没有
  }else if( (!x_exist) & (!y_exist)){  
    power_plot_v2 <- ggplot(data_test,aes(x=method,y=power,fill=method))+
      geom_bar(stat="identity",width=0.5)+
      theme_bw()+
      theme(panel.background = element_blank()) +
      scale_fill_manual(values=color_setting)+
      scale_y_continuous(limits = c(0, 1), 
                         breaks = c(0.95))+
      #theme(panel.grid =element_blank()) + 
      theme(axis.text= element_text(size = zihao2),
            axis.title.x = element_blank(),  # 移除 x 轴标题
            axis.text.x = element_blank(),   # 移除 x 轴刻度
            axis.ticks.x = element_blank(),  # 移除 x 轴刻度线
            
            axis.title.y = element_blank(),  # 移除 y 轴标题
            axis.text.y = element_blank(),   # 移除 y 轴刻度
            axis.ticks.y = element_blank(),  # 移除 y 轴刻度线
            
            #axis.title.y=element_text(size=zihao1),
            plot.title=element_text(size=23)) +
      labs(y = '') +
      theme(legend.title = element_text(size = zihao3),
            legend.text=element_text(size=zihao3))+  
      geom_text(mapping = aes(label = power),vjust = -0.5,size=zihao4)
    
    cover_plot_v2 <- ggplot(data_test,aes(x=method,y=cover_ratio,fill=method))+
      geom_bar(stat="identity",width=0.5)+
      theme_bw()+
      theme(panel.background = element_blank()) +
      scale_fill_manual(values=color_setting)+
      scale_y_continuous(limits = c(0, 1), 
                         breaks = c(0.95))+
      #theme(panel.grid =element_blank()) + 
      theme(axis.text= element_text(size = zihao2),
            #axis.title.x = element_blank(),  # 移除 x 轴标题
            axis.title.x=element_text(size=zihao1),
            axis.text.x = element_blank(),   # 移除 x 轴刻度
            axis.ticks.x = element_blank(),  # 移除 x 轴刻度线
            
            axis.title.y = element_blank(),  # 移除 y 轴标题
            axis.text.y = element_blank(),   # 移除 y 轴刻度
            axis.ticks.y = element_blank(),  # 移除 y 轴刻度线
            
            #axis.title.y=element_text(size=zihao1),
            plot.title=element_text(size=23)) +
      labs(y = '',
           x = x_label) +
      theme(legend.title = element_text(size = zihao3),
            legend.text=element_text(size=zihao3))+  
      geom_text(mapping = aes(label = cover_ratio),vjust = -0.5,size=zihao4)
    
    # pyt_combine3 <- ggarrange(power_plot_v2,cover_plot_v2,
    #                           ncol=1,nrow=2,
    #                           #labels = c('(A)','(B)','(C)','(D)','(E)','(F)'),
    #                           common.legend = TRUE,
    #                           legend="none")
  }
  
  power_plot_v2 <- power_plot_v2 + geom_hline(aes(yintercept = 0.95), linetype = "dashed", color = "darkred")
  cover_plot_v2 <- cover_plot_v2 + geom_hline(aes(yintercept = 0.95), linetype = "dashed", color = "darkred")
  if(y_exist){
    pyt_combine3 <- ggarrange(power_plot_v2,cover_plot_v2,
                              ncol=1,nrow=2,
                              labels = lable_list,
                              hjust= -2.6,
                              vjust = 2,
                              common.legend = TRUE,
                              legend="none")
  }else{
    pyt_combine3 <- ggarrange(power_plot_v2,cover_plot_v2,
                              ncol=1,nrow=2,
                              labels = lable_list,
                              hjust= -1,
                              vjust = 2,
                              common.legend = TRUE,
                              legend="none")
  }
  
  #pyt_combine3
  
  return(pyt_combine3)
  #+ geom_hline(aes(yintercept = 1), linetype = "dashed", color = "darkred")
}

power_plot <- function(data,color_setting,x_exist=TRUE,y_exist=TRUE,x_label,y_label,axis,zihao4){  
  data_test <- data[,c(1,axis)]
  colnames(data_test) <- c('method','power')
  #x轴y轴都有
  if(x_exist & y_exist){                     
    power_plot_v2 <- ggplot(data_test,aes(x=method,y=power,fill=method))+
      geom_bar(stat="identity",width=0.5)+
      theme_bw()+
      theme(panel.background = element_blank()) +
      scale_fill_manual(values=color_setting)+
      scale_y_continuous(limits = c(0, 1), 
                         breaks = c(0,0.2,0.4,0.6,0.8,1))+
      #theme(panel.grid =element_blank()) + 
      theme(#axis.text= element_text(size = 18),
        axis.text.x= element_text(angle = 45, hjust = 1,size = 9),
        axis.title.x = element_text(size=19),
        
        axis.title.y = element_text(size=19),
        
        plot.title=element_text(size=23)) +
      labs(y = y_label
           ,x = x_label) +
      theme(legend.title = element_text(size = 15),
            legend.text=element_text(size=15))+  
      geom_text(mapping = aes(label = power),vjust = -0.5,size=zihao4)
    #只有x轴没有y轴
  }else if(x_exist & (!y_exist) ){  
    power_plot_v2 <- ggplot(data_test,aes(x=method,y=power,fill=method))+
      geom_bar(stat="identity",width=0.5)+
      theme_bw()+
      theme(panel.background = element_blank()) +
      scale_fill_manual(values=color_setting)+
      scale_y_continuous(limits = c(0, 1), 
                         breaks = c(0,0.2,0.4,0.6,0.8,1))+
      #theme(panel.grid =element_blank()) + 
      theme(#axis.text= element_text(size = 18),
        axis.text.x= element_text(angle = 45, hjust = 1,size = 9),
        axis.title.x = element_text(size=19),
        
        axis.text.y = element_blank(),   
        axis.ticks.y = element_blank(),  
        axis.title.y=element_blank(),
        
        plot.title=element_text(size=23)) +
      labs(y = y_label
           ,x = x_label) +
      theme(legend.title = element_text(size = 15),
            legend.text=element_text(size=15))+  
      geom_text(mapping = aes(label = power),vjust = -0.5,size=zihao4)
    #没有x轴只有y轴
  }else if( (!x_exist) & y_exist){  
    power_plot_v2 <- ggplot(data_test,aes(x=method,y=power,fill=method))+
      geom_bar(stat="identity",width=0.5)+
      theme_bw()+
      theme(panel.background = element_blank()) +
      scale_fill_manual(values=color_setting)+
      scale_y_continuous(limits = c(0, 1), 
                         breaks = c(0,0.2,0.4,0.6,0.8,1))+
      #theme(panel.grid =element_blank()) + 
      theme(#axis.text= element_text(size = 18),
        
        axis.text.x = element_blank(),  
        axis.ticks.x = element_blank(),  
        axis.title.x = element_text(size=19), 
        
        axis.title.y=element_text(size=19),
        
        plot.title=element_text(size=23)) +
      labs(y = y_label
           ,x = x_label) +
      theme(legend.title = element_text(size = 15),
            legend.text=element_text(size=15))+  
      geom_text(mapping = aes(label = power),vjust = -0.5,size=zihao4)
    # x轴没有y轴也没有
  }else if( (!x_exist) & (!y_exist)){  
    power_plot_v2 <- ggplot(data_test,aes(x=method,y=power,fill=method))+
      geom_bar(stat="identity",width=0.5)+
      theme_bw()+
      theme(panel.background = element_blank()) +
      scale_fill_manual(values=color_setting)+
      scale_y_continuous(limits = c(0, 1), 
                         breaks = c(0,0.2,0.4,0.6,0.8,1))+
      #theme(panel.grid =element_blank()) + 
      theme(#axis.text= element_text(size = 18),
        
        axis.title.x = element_text(size=19),
        axis.text.x = element_blank(),   # 移除 x 轴刻度
        axis.ticks.x = element_blank(),  # 移除 x 轴刻度线
        
        axis.title.y = element_blank(),  # 移除 y 轴标题
        axis.text.y = element_blank(),   # 移除 y 轴刻度
        axis.ticks.y = element_blank(),  # 移除 y 轴刻度线
        
        #axis.title.y=element_text(size=19),
        plot.title=element_text(size=23)) +
      labs(y = y_label
           ,x = x_label) +
      theme(legend.title = element_text(size = 15),
            legend.text=element_text(size=15))+  
      geom_text(mapping = aes(label = power),vjust = -0.5,size=zihao4)
  }
  return(power_plot_v2)
  
}
#------------------
# for (i_z in 1:length(save_power_list)) {
# 
#   df_2m <- save_power_list[[i_z]]
#   df_2m$method <- factor(df_2m$method, levels = c("IVW","Egger" ,"Weighted median" ,      "Weighted mode"  ,       "MR-RAPS"      ,
#                                                   "MR-conmix","MR-Conmix","MR-Robust"   ,          "MR-Lasso"  ,  "TLMR-OR", "TLMR-IPW"  ,"TLMR-AIPW"))
#   # df_2m$method[df_2m$method=='Trans-OR'] <- 'TLMR-OR'
#   # df_2m$method[df_2m$method=='Trans-IPW'] <- 'TLMR-IPW'
#   # df_2m$method[df_2m$method=='Trans-AIPW'] <- 'TLMR-AIPW'
#   df_2m$method[df_2m$method=='MR-conmix'] <- 'MR-Conmix'
#   df_2m$method <- droplevels(df_2m$method) 
# 
#   save_power_list[[i_z]] <- df_2m
# 
# }
# 
# for (i_z in 1:length(save_power_list2)) {
# 
#   df_2m <- save_power_list2[[i_z]]
#   df_2m$method <- factor(df_2m$method, levels = c("IVW","Egger" ,"Weighted median" ,      "Weighted mode"  ,       "MR-RAPS"      ,
#                                                   "MR-conmix","MR-Conmix","MR-Robust"   ,          "MR-Lasso"  ,  "TLMR-OR", "TLMR-IPW"  ,"TLMR-AIPW"))
#   # df_2m$method[df_2m$method=='Trans-OR'] <- 'TLMR-OR'
#   # df_2m$method[df_2m$method=='Trans-IPW'] <- 'TLMR-IPW'
#   # df_2m$method[df_2m$method=='Trans-AIPW'] <- 'TLMR-AIPW'
#   df_2m$method[df_2m$method=='MR-conmix'] <- 'MR-Conmix'
#   df_2m$method <- droplevels(df_2m$method) 
# 
#   save_power_list2[[i_z]] <- df_2m
# 
# }
#--------------------------------------------------------power图1
data_power_list <- vector('list',length = length(save_power_list))  #最终画图数据
for(i_p in 1:length(save_power_list)){
  test_df <- save_power_list[[i_p]]
  data_test <- data.frame(method=method_names,power=rep(NA,len_methods),cover_ratio=rep(NA,len_methods),n_0=rep(n_0_list2[i_p],len_methods))
  data_test$method <- factor(data_test$method,levels = method_names)
  
  for (i in 1:len_methods) {
    test_df_process <- test_df[test_df$method == method_names[i],]
    #------错误处理
    process_p_power <- test_df_process$p_power
    process_p_cover <- test_df_process$p_cover
    
    process_p_power_final <- subset(process_p_power,process_p_power < 1 )
    process_p_cover_final <- subset(process_p_cover,process_p_cover < 1 )
    
    test_times_power <- length(process_p_power_final)
    test_times_cover <- length(process_p_cover_final)
    #------
    data_test$power[i] <- round(sum(process_p_power_final < 0.05)/test_times_power,3)
    data_test$cover_ratio[i] <- round(sum(process_p_cover_final > 0.05)/test_times_cover,3)
  }
  data_power_list[[i_p]] <- data_test
}                                                                                                                                           

data_power_df <- do.call(rbind,data_power_list)
fwrite(data_power_df,file=paste(path_output_tables,'/data_power_hand_',date_mark,'_',ver,'.csv',sep = ''))
#------------
# path_read <- 'output/categorical/tables/2024-11-03-v5/data_power_hand_2024-11-03_v5.csv'
# path_output_tables <- 'output/categorical/tables/2024-11-03-v5'
# data_power_df <- fread(path_read)
#--------------------------------------------------------power图2
data_power_list2 <- vector('list',length = length(save_power_list2))  #最终画图数据
for(i_p in 1:length(save_power_list2)){
  test_df <- save_power_list2[[i_p]]
  data_test <- data.frame(method=method_names,power=rep(NA,len_methods),cover_ratio=rep(NA,len_methods),n_0=rep(n_0_list2[i_p],len_methods))
  data_test$method <- factor(data_test$method,levels = method_names)
  
  for (i in 1:len_methods) {
    test_df_process <- test_df[test_df$method == method_names[i],]
    #------错误处理
    process_p_power <- test_df_process$p_power
    process_p_cover <- test_df_process$p_cover
    
    process_p_power_final <- subset(process_p_power,process_p_power < 1 )
    process_p_cover_final <- subset(process_p_cover,process_p_cover < 1 )
    
    test_times_power <- length(process_p_power_final)
    test_times_cover <- length(process_p_cover_final)
    #------
    data_test$power[i] <- round(sum(process_p_power_final < 0.05)/test_times_power,3)
    data_test$cover_ratio[i] <- round(sum(process_p_cover_final > 0.05)/test_times_cover,3)
  }
  data_power_list2[[i_p]] <- data_test
}

data_power_df2 <- do.call(rbind,data_power_list2)
fwrite(data_power_df2,file=paste(path_output_tables,'/data_power_boot_',date_mark,'_',ver,'.csv',sep = ''))

#------------------------------------------------------
# #data,color_setting,x_exist=TRUE,y_exist=TRUE,x_label,y_label,axis
# power_plot_1 <- power_plot(data=data_power_list[[1]]
#                            ,color_setting = color_setting
#                            ,x_exist = FALSE
#                            ,y_exist = TRUE
#                            ,x_label = paste0('Sample size = ',n_0_list2[1])
#                            ,y_label = 'Power'
#                            ,axis = 2
#                            ,zihao4 = zihao4
# )
# power_plot_2 <- power_plot(data=data_power_list[[2]]
#                            ,color_setting = color_setting
#                            ,x_exist = FALSE
#                            ,y_exist = FALSE
#                            ,x_label = paste0('Sample size = ',n_0_list2[2])
#                            ,y_label = 'Power'
#                            ,axis = 2
#                            ,zihao4 = zihao4
# )
# power_plot_3 <- power_plot(data=data_power_list[[3]]
#                            ,color_setting = color_setting
#                            ,x_exist = FALSE
#                            ,y_exist = FALSE
#                            ,x_label = paste0('Sample size = ',n_0_list2[3])
#                            ,y_label = 'Power'
#                            ,axis = 2
#                            ,zihao4 = zihao4
# )
# power_plot_4 <- power_plot(data=data_power_list[[4]]
#                            ,color_setting = color_setting
#                            ,x_exist = FALSE
#                            ,y_exist = TRUE
#                            ,x_label = paste0('Sample size = ',n_0_list2[4])
#                            ,y_label = 'Power'
#                            ,axis = 2
#                            ,zihao4 = zihao4
# )
# power_plot_5 <- power_plot(data=data_power_list[[5]]
#                            ,color_setting = color_setting,
#                            x_exist = FALSE
#                            ,y_exist = FALSE
#                            ,x_label = paste0('Sample size = ',n_0_list2[5])
#                            ,y_label = 'Power'
#                            ,axis = 2
#                            ,zihao4 = zihao4
# )
# power_plot_6 <- power_plot(data=data_power_list[[6]]
#                            ,color_setting = color_setting
#                            ,x_exist = FALSE
#                            ,y_exist = FALSE
#                            ,x_label = paste0('Sample size = ',n_0_list2[6])
#                            ,y_label = 'Power'
#                            ,axis = 2
#                            ,zihao4 = zihao4
# )
# 
# power_plot_all <- ggarrange(power_plot_1,power_plot_2
#                             ,power_plot_3,power_plot_4
#                             ,power_plot_5,power_plot_6
#                             ,ncol=3,nrow=2,
#                             common.legend = TRUE,
#                             legend="right")
# #print(power_plot_all)
# #---------------------
# cover_plot_1 <- power_plot(data=data_power_list[[1]]
#                            ,color_setting = color_setting
#                            ,x_exist = FALSE
#                            ,y_exist = TRUE
#                            ,x_label = paste0('Sample size = ',n_0_list2[1])
#                            ,y_label = 'Cover Ratio'
#                            ,axis = 3
#                            ,zihao4 = zihao4
# )
# cover_plot_2 <- power_plot(data=data_power_list[[2]]
#                            ,color_setting = color_setting
#                            ,x_exist = FALSE
#                            ,y_exist = FALSE
#                            ,x_label = paste0('Sample size = ',n_0_list2[2])
#                            ,y_label = 'Cover Ratio'
#                            ,axis = 3
#                            ,zihao4 = zihao4
# )
# cover_plot_3 <- power_plot(data=data_power_list[[3]]
#                            ,color_setting = color_setting
#                            ,x_exist = FALSE
#                            ,y_exist = FALSE
#                            ,x_label = paste0('Sample size = ',n_0_list2[3])
#                            ,y_label = 'Cover Ratio'
#                            ,axis = 3
#                            ,zihao4 = zihao4
# )
# cover_plot_4 <- power_plot(data=data_power_list[[4]]
#                            ,color_setting = color_setting
#                            ,x_exist = FALSE
#                            ,y_exist = TRUE
#                            ,x_label = paste0('Sample size = ',n_0_list2[4])
#                            ,y_label = 'Cover Ratio'
#                            ,axis = 3
#                            ,zihao4 = zihao4
# )
# cover_plot_5 <- power_plot(data=data_power_list[[5]]
#                            ,color_setting = color_setting,
#                            x_exist = FALSE
#                            ,y_exist = FALSE
#                            ,x_label = paste0('Sample size = ',n_0_list2[5])
#                            ,y_label = 'Cover Ratio'
#                            ,axis = 3
#                            ,zihao4 = zihao4
# )
# cover_plot_6 <- power_plot(data=data_power_list[[6]]
#                            ,color_setting = color_setting
#                            ,x_exist = FALSE
#                            ,y_exist = FALSE
#                            ,x_label = paste0('Sample size = ',n_0_list2[6])
#                            ,y_label = 'Cover Ratio'
#                            ,axis = 3
#                            ,zihao4 = zihao4
# )
# 
# cover_plot_all <- ggarrange(cover_plot_1,cover_plot_2
#                             ,cover_plot_3,cover_plot_4
#                             ,cover_plot_5,cover_plot_6
#                             ,ncol=3,nrow=2,
#                             # labels = c('(A)','(B)','(C)','(D)','(E)','(F)'),
#                             common.legend = TRUE,
#                             legend="right")

#-----------------------------------------------------
#组合图形1
#########
power_cover_plot1 <- power_cover_plot(data_test = data_power_list[[1]]
                                      ,color_setting = color_setting
                                      ,x_exist = TRUE
                                      ,y_exist = TRUE
                                      ,x_label = paste0('Sample size = ',n_0_list2[1])
                                      ,zihao1 = zihao1
                                      ,zihao2 = zihao2
                                      ,zihao3 = zihao3
                                      ,zihao4 = zihao4
                                      ,lable_list = c("(A)","(D)")
) 
# + theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1, linetype = "dashed")
#           ,plot.margin = margin(t = 10, r = 10, b = 10, l = 10) )

power_cover_plot2 <- power_cover_plot(data_test = data_power_list[[2]]
                                      ,color_setting = color_setting
                                      ,x_exist = TRUE
                                      ,y_exist = FALSE
                                      ,x_label = paste0('Sample size = ',n_0_list2[2])
                                      ,zihao1 = zihao1
                                      ,zihao2 = zihao2
                                      ,zihao3 = zihao3
                                      ,zihao4 = zihao4
                                      ,lable_list = c("(B)","(E)")
)
# + theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1, linetype = "dashed")
#          ,plot.margin = margin(t = 10, r = 10, b = 10, l = 10))

power_cover_plot3 <- power_cover_plot(data_test = data_power_list[[3]]
                                      ,color_setting = color_setting
                                      ,x_exist = TRUE
                                      ,y_exist = FALSE
                                      ,x_label = paste0('Sample size = ',n_0_list2[3])
                                      ,zihao1 = zihao1
                                      ,zihao2 = zihao2
                                      ,zihao3 = zihao3
                                      ,zihao4 = zihao4
                                      ,lable_list = c("(C)","(F)")
)
# + theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1, linetype = "dashed")
#          ,plot.margin = margin(t = 10, r = 10, b = 10, l = 10))

# power_cover_plot4 <- power_cover_plot(data_test = data_power_list[[4]]
#                                       ,color_setting = color_setting
#                                       ,x_exist = TRUE
#                                       ,y_exist = TRUE
#                                       ,x_label = paste0('Sample size = ',n_0_list2[4])
#                                       ,zihao1 = zihao1
#                                       ,zihao3 = zihao3
#                                       ,zihao4 = zihao4
# )+ theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1, linetype = "dashed")
#          ,plot.margin = margin(t = 10, r = 10, b = 10, l = 10))
# power_cover_plot5 <- power_cover_plot(data_test = data_power_list[[5]]
#                                       ,color_setting = color_setting
#                                       ,x_exist = TRUE
#                                       ,y_exist = FALSE
#                                       ,x_label = paste0('Sample size = ',n_0_list2[5])
#                                       ,zihao1 = zihao1
#                                       ,zihao3 = zihao3
#                                       ,zihao4 = zihao4
# )+ theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1, linetype = "dashed")
#          ,plot.margin = margin(t = 10, r = 10, b = 10, l = 10))
# power_cover_plot6 <- power_cover_plot(data_test = data_power_list[[6]]
#                                       ,color_setting = color_setting
#                                       ,x_exist = TRUE
#                                       ,y_exist = FALSE
#                                       ,x_label = paste0('Sample size = ',n_0_list2[6])
#                                       ,zihao1 = zihao1
#                                       ,zihao3 = zihao3
#                                       ,zihao4 = zihao4
# )+ theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1, linetype = "dashed")
#          ,plot.margin = margin(t = 10, r = 10, b = 10, l = 10))

# 创建图例的空白图
p_with_legend <- ggplot(data_power_list[[1]], aes(x=method,y=power,fill=method) ) +
  geom_bar(stat="identity",width=0.5) +
  scale_fill_manual(values=color_setting)+
  theme(legend.title = element_text(size = 15),
        legend.text=element_text(size=15))+ 
  theme_void() +  # 不显示任何图形内容
  theme(legend.position = "right")  # 只保留图例
legend_only <- cowplot::get_plot_component(p_with_legend, 'guide-box-right', return_all = TRUE)
# #组合图形
# power_cover_combined_plot <- ggarrange(power_cover_plot1,power_cover_plot2
#                                        ,power_cover_plot3,power_cover_plot4
#                                        ,power_cover_plot5,power_cover_plot6
#                                        ,ncol=3,nrow=2,
#                                        common.legend = TRUE,
#                                        legend="none")
# power_cover_plot_with_legend <- cowplot::plot_grid(
#   power_cover_combined_plot,  # 主图
#   legend_only,  # 单独的图例
#   ncol = 2,  # 两列布局
#   rel_widths = c(4, 0.5)  # 调整图形与图例的宽度比例
# )

#组合图形
power_cover_combined_plot <- ggarrange(power_cover_plot1,power_cover_plot2
                                       ,power_cover_plot3
                                       ,ncol=3,nrow=1
                                       ,common.legend = TRUE,
                                       # labels = c('   (A)','(B)','(C)'),
                                       # # labels.x = 1,
                                       # hjust= -1.6,
                                       # vjust = 2,
                                       legend="none"
)
power_cover_combined_plot
power_cover_plot_with_legend <- cowplot::plot_grid(
  power_cover_combined_plot,  # 主图
  legend_only,  # 单独的图例
  ncol = 2,  # 两列布局
  rel_widths = c(4, 0.5)  # 调整图形与图例的宽度比例
)
power_cover_plot_with_legend
#-----------------------------------------------------
#组合图形2
#########
power_cover_plot1_2 <- power_cover_plot(data_test = data_power_list2[[1]]
                                        ,color_setting = color_setting
                                        ,x_exist = TRUE
                                        ,y_exist = TRUE
                                        ,x_label = paste0('Sample size of two samples = ',n_0_list2[1])
                                        ,zihao1 = zihao1
                                        ,zihao2 = zihao2
                                        ,zihao3 = zihao3
                                        ,zihao4 = zihao4
                                        ,lable_list = c("(A)","(D)")
) 

power_cover_plot2_2 <- power_cover_plot(data_test = data_power_list2[[2]]
                                        ,color_setting = color_setting
                                        ,x_exist = TRUE
                                        ,y_exist = FALSE
                                        ,x_label = paste0('Sample size of two samples = ',n_0_list2[2])
                                        ,zihao1 = zihao1
                                        ,zihao2 = zihao2
                                        ,zihao3 = zihao3
                                        ,zihao4 = zihao4
                                        ,lable_list = c("(D)","(E)")
)
power_cover_plot3_2 <- power_cover_plot(data_test = data_power_list2[[3]]
                                        ,color_setting = color_setting
                                        ,x_exist = TRUE
                                        ,y_exist = FALSE
                                        ,x_label = paste0('Sample size of two samples = ',n_0_list2[3])
                                        ,zihao1 = zihao1
                                        ,zihao2 = zihao2
                                        ,zihao3 = zihao3
                                        ,zihao4 = zihao4
                                        ,lable_list = c("(C)","(F)")
)

# power_cover_plot4_2 <- power_cover_plot(data_test = data_power_list2[[4]]
#                                         ,color_setting = color_setting
#                                         ,x_exist = TRUE
#                                         ,y_exist = TRUE
#                                         ,x_label = paste0('Sample size of two samples = ',n_0_list2[4])
#                                         ,zihao1 = zihao1
#                                         ,zihao3 = zihao3
#                                         ,zihao4 = zihao4
# )+ theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1, linetype = "dashed")
#          ,plot.margin = margin(t = 10, r = 10, b = 10, l = 10))
# power_cover_plot5_2 <- power_cover_plot(data_test = data_power_list2[[5]]
#                                         ,color_setting = color_setting
#                                         ,x_exist = TRUE
#                                         ,y_exist = FALSE
#                                         ,x_label = paste0('Sample size of two samples = ',n_0_list2[5])
#                                         ,zihao1 = zihao1
#                                         ,zihao3 = zihao3
#                                         ,zihao4 = zihao4
# )+ theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1, linetype = "dashed")
#          ,plot.margin = margin(t = 10, r = 10, b = 10, l = 10))
# power_cover_plot6_2 <- power_cover_plot(data_test = data_power_list2[[6]]
#                                         ,color_setting = color_setting
#                                         ,x_exist = TRUE
#                                         ,y_exist = FALSE
#                                         ,x_label = paste0('Sample size of two samples = ',n_0_list2[6])
#                                         ,zihao1 = zihao1
#                                         ,zihao3 = zihao3
#                                         ,zihao4 = zihao4
# )+ theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1, linetype = "dashed")
#          ,plot.margin = margin(t = 10, r = 10, b = 10, l = 10))

#组合图形
power_cover_combined_plot2 <- ggarrange(power_cover_plot1_2,power_cover_plot2_2
                                        ,power_cover_plot3_2
                                        ,ncol=3,nrow=1,
                                        common.legend = TRUE,
                                        legend="none")
power_cover_plot_with_legend2 <- cowplot::plot_grid(
  power_cover_combined_plot2,  # 主图
  legend_only,  # 单独的图例
  ncol = 2,  # 两列布局
  rel_widths = c(4, 0.5)  # 调整图形与图例的宽度比例
)

#power_cover_plot_with_legend
#------------------------- 
# ggsave(paste('power_',date_mark,'_',ver,'.pdf',sep = ''),plot = power_plot_all, bg="white", path = path_output_plots, height = 9,width = 20,dpi = 300)
# ggsave(paste('cover_ratio_',date_mark,'_',ver,'.pdf',sep = ''),plot = cover_plot_all, bg="white", path = path_output_plots, height = 9,width = 20,dpi = 300)
# ggsave(paste('power+cover_hand_',date_mark,'_',ver,'.pdf',sep = ''),plot = power_cover_plot_with_legend, bg="white", path = path_output_plots, height = 15,width = 20,dpi = 300)
# 
# ggsave(paste('power_',date_mark,'_',ver,'.png',sep = ''),plot = power_plot_all, bg="white", path = path_output_plots, height = 9,width = 20,dpi = 300)
# ggsave(paste('cover_ratio_',date_mark,'_',ver,'.png',sep = ''),plot = cover_plot_all, bg="white", path = path_output_plots, height = 9,width = 20,dpi = 300)
# ggsave(paste('power+cover_hand_',date_mark,'_',ver,'.png',sep = ''),plot = power_cover_plot_with_legend, bg="white", path = path_output_plots, height = 12,width = 18,dpi = 600)
# 

ggsave(paste('power+cover_boot_',date_mark,'_',ver,'.pdf',sep = ''),plot = power_cover_plot_with_legend2, bg="white", path = path_output_plots, height = 13.5,width = 18,dpi = 500)
ggsave(paste('power+cover_boot_',date_mark,'_',ver,'.png',sep = ''),plot = power_cover_plot_with_legend2, bg="white", path = path_output_plots, height = 13.5,width = 18,dpi = 500)
ggsave(paste('power+cover_hand_',date_mark,'_',ver,'.pdf',sep = ''),plot = power_cover_plot_with_legend, bg="white", path = path_output_plots, height = 13.5,width = 18,dpi = 500)
ggsave(paste('power+cover_hand_',date_mark,'_',ver,'.png',sep = ''),plot = power_cover_plot_with_legend, bg="white", path = path_output_plots, height = 13.5,width = 18,dpi = 500)
#-----------------------------------------------

