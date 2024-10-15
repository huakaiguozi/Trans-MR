shapes_1 <- c(1:len_methods)-1

colour_table1 <- c(  '#e6194B', '#3cb44b', '#ffe119', '#4363d8', '#f58231'
                     , '#42d4f4', '#f032e6', '#fabed4', '#469990', '#dcbeff'
                     , '#9A6324', '#fffac8', '#800000', '#aaffc3', '#000075'
                     , '#a9a9a9', '#ffffff', '#000000')

colour_table2 <- c(  '#e6194B', '#3cb44b', '#ffe119', '#4363d8', '#f58231'
                     , '#911eb4', '#42d4f4', '#f032e6', '#bfef45', '#fabed4'
                     , '#469990', '#dcbeff', '#9A6324', '#fffac8', '#800000'
                     , '#aaffc3', '#808000', '#ffd8b1', '#000075', '#a9a9a9'
                     , '#ffffff', '#000000')


colour_table_select <- c( '#e6194B', '#ffe119', '#4363d8', '#f58231'
                          , '#911eb4', '#42d4f4', '#f032e6'
                          , '#469990', '#9A6324', '#fffac8', '#800000'
                          , '#808000', '#000075', '#a9a9a9'
                          , '#000000','#fabed4','#dcbeff','#ffd8b1')

#coloo_1 <- c('#aaffc3','#bfef45','#3cb44b', colour_table_select[1:(len_methods-3)])
coloo_1 <- c(colour_table2[1:len_methods])
#---------
#& save_df2$n1 == n_1_list[k_min] 
data_gamma <- save_df2[save_df2$heterogeneity == heter_level_list[k_heter] & save_df2$n_0 == n_0_list[k_n0] & save_df2$inter_ratio == ratio_inter_list[k_inter] & save_df2$p_Z == p_Z_list[k_pz] & save_df2$inter_strength==inter_strength_list[k_inter_strength]& save_df2$alpha_modifier == alpha_list[k_alpha] & save_df2$eta_modifier == eta_list[k_eta], ]
data_heterogeneity <- save_df2[save_df2$gamma == gamma_list[k_gamma] & save_df2$n_0 == n_0_list[k_n0] & save_df2$inter_ratio == ratio_inter_list[k_inter] & save_df2$p_Z == p_Z_list[k_pz] & save_df2$inter_strength==inter_strength_list[k_inter_strength]& save_df2$alpha_modifier == alpha_list[k_alpha] & save_df2$eta_modifier == eta_list[k_eta], ]
data_n0 <- save_df2[save_df2$heterogeneity == heter_level_list[k_heter] & save_df2$gamma == gamma_list[k_gamma] & save_df2$inter_ratio == ratio_inter_list[k_inter] & save_df2$p_Z == p_Z_list[k_pz] & save_df2$inter_strength==inter_strength_list[k_inter_strength]& save_df2$alpha_modifier == alpha_list[k_alpha] & save_df2$eta_modifier == eta_list[k_eta], ]
data_inter_ratio <- save_df2[save_df2$heterogeneity == heter_level_list[k_heter] & save_df2$gamma == gamma_list[k_gamma]& save_df2$n_0 == n_0_list[k_n0] & save_df2$p_Z == p_Z_list[k_pz] & save_df2$inter_strength==inter_strength_list[k_inter_strength]& save_df2$alpha_modifier == alpha_list[k_alpha] & save_df2$eta_modifier == eta_list[k_eta], ]
data_pz <- save_df2[save_df2$heterogeneity == heter_level_list[k_heter] & save_df2$gamma == gamma_list[k_gamma]& save_df2$n_0 == n_0_list[k_n0] & save_df2$inter_ratio == ratio_inter_list[k_inter] & save_df2$inter_strength==inter_strength_list[k_inter_strength]& save_df2$alpha_modifier == alpha_list[k_alpha] & save_df2$eta_modifier == eta_list[k_eta], ]
data_inter_strength <- save_df2[save_df2$heterogeneity == heter_level_list[k_heter] & save_df2$gamma == gamma_list[k_gamma]& save_df2$n_0 == n_0_list[k_n0] & save_df2$p_Z == p_Z_list[k_pz] & save_df2$inter_ratio == ratio_inter_list[k_inter]& save_df2$alpha_modifier == alpha_list[k_alpha2] & save_df2$eta_modifier == eta_list[k_eta], ]
data_alpha <- save_df2[save_df2$heterogeneity == heter_level_list[k_heter] & save_df2$gamma == gamma_list[k_gamma]& save_df2$n_0 == n_0_list[k_n0] & save_df2$p_Z == p_Z_list[k_pz] & save_df2$inter_ratio == ratio_inter_list[k_inter] & save_df2$inter_strength==inter_strength_list[k_inter_strength] & save_df2$eta_modifier == eta_list[k_eta], ]
data_eta <- save_df2[save_df2$heterogeneity == heter_level_list[k_heter] & save_df2$gamma == gamma_list[k_gamma]& save_df2$n_0 == n_0_list[k_n0] & save_df2$p_Z == p_Z_list[k_pz] & save_df2$inter_ratio == ratio_inter_list[k_inter]& save_df2$inter_strength==inter_strength_list[k_inter_strength] & save_df2$alpha_modifier == alpha_list[k_alpha], ]
#----
data_pz <- data_pz[data_pz$p_Z != p_Z_fix,]
#------------------
ylim <- 1
pyt1 <- ggplot(data=data_gamma, aes(x=gamma,y=num, colour=method, group=method)) +
  geom_point(size=2,aes(shape=method,colour=method)) +
  labs(x='Relative magnitude of horizontal pleiotropy effect')+
  scale_shape_manual(values=shapes_1) +
  geom_line(linewidth=1.0,aes( colour=method)) +
  scale_colour_manual(values=coloo_1)+
  facet_grid(vars(),vars(index),scales="free_x")+
  theme_bw() +
  theme(plot.title = element_text(size = 15,hjust = 0.5,vjust = 0.5),
        axis.title.y = element_blank()) +
  ylim(-ylim,ylim)

pyt2 <- ggplot(data=data_heterogeneity, aes(x=heterogeneity, y=num, colour=method, group=method)) +
  geom_point(size=2,aes(shape=method,colour=method)) +
  labs(x='heterogeneity level between the two samples')+
  scale_shape_manual(values=shapes_1) +
  geom_line(size=1.0,aes(colour=method)) +
  scale_colour_manual(values=coloo_1)+
  facet_grid( vars(), vars(index),scales="free_x")+
  #ylab("heterogeneity of valid instruments") +
  theme_bw() +
  theme(plot.title = element_text(size = 15,hjust = 0.5,vjust = 0.5),
        axis.title.y = element_blank()) +
  ylim(-ylim,ylim)

pyt3 <- ggplot(data=data_n0, aes(x=n_0,y=num, colour=method, group=method)) +
  geom_point(size=2,aes(shape=method,colour=method)) +
  labs(x='sample size of two datasets')+
  #geom_point(size=2,aes(x=heterogeneity,shape=method,colour=method)) +
  scale_shape_manual(values=shapes_1) +
  geom_line(linewidth=1.0,aes( colour=method)) +
  #geom_line(linewidth=0.7,aes(colour=method)) +
  scale_colour_manual(values=coloo_1)+
  facet_grid(vars(),vars(index),scales="free_x")+
  theme_bw() +
  theme(plot.title = element_text(size = 15,hjust = 0.5,vjust = 0.5),
        axis.title.y = element_blank()) +
  ylim(-ylim,ylim)

# pyt4 <- ggplot(data=data_n1, aes(x=n1, y=num, colour=method, group=method)) +
#   geom_point(size=2,aes(shape=method,colour=method)) +
#   scale_shape_manual(values=shapes_1) +
#   geom_line(size=0.7,aes(colour=method)) +
#   scale_colour_manual(values=coloo_1)+
#   facet_grid( vars(), vars(index),scales="free_x")+
#   #ylab("heterogeneity of valid instruments") +
#   theme_bw() +
#   theme(plot.title = element_text(size = 15,hjust = 0.5,vjust = 0.5),
#         axis.title.y = element_blank()) +
#   ylim(-ylim,ylim)

pyt5 <- ggplot(data=data_inter_ratio, aes(x=inter_ratio,y=num, colour=method, group=method)) +
  geom_point(size=2,aes(shape=method,colour=method)) +
  labs(x='proportion of instruments with interactive effects')+
  #geom_point(size=2,aes(x=heterogeneity,shape=method,colour=method)) +
  scale_shape_manual(values=shapes_1) +
  geom_line(linewidth=1.0,aes( colour=method)) +
  #geom_line(linewidth=0.7,aes(colour=method)) +
  scale_colour_manual(values=coloo_1)+
  facet_grid(vars(),vars(index),scales="free_x")+
  theme_bw() +
  theme(plot.title = element_text(size = 15,hjust = 0.5,vjust = 0.5),
        axis.title.y = element_blank()) +
  ylim(-ylim,ylim)

pyt6 <- ggplot(data=data_pz, aes(x=p_Z,y=num, colour=method, group=method)) +
  geom_point(size=2,aes(shape=method,colour=method)) +
  labs(x='Number of instrumental variables')+
  #geom_point(size=2,aes(x=heterogeneity,shape=method,colour=method)) +
  scale_shape_manual(values=shapes_1) +
  geom_line(linewidth=1.0,aes( colour=method)) +
  #geom_line(linewidth=0.7,aes(colour=method)) +
  scale_colour_manual(values=coloo_1)+
  facet_grid(vars(),vars(index),scales="free_x")+
  theme_bw() +
  theme(plot.title = element_text(size = 15,hjust = 0.5,vjust = 0.5),
        axis.title.y = element_blank()) +
  ylim(-ylim,ylim)

pyt7 <- ggplot(data=data_inter_strength, aes(x=inter_strength,y=num, colour=method, group=method)) +
  geom_point(size=2,aes(shape=method,colour=method)) +
  labs(x='Number of instrumental variables')+
  #geom_point(size=2,aes(x=heterogeneity,shape=method,colour=method)) +
  scale_shape_manual(values=shapes_1) +
  geom_line(linewidth=1.0,aes( colour=method)) +
  #geom_line(linewidth=0.7,aes(colour=method)) +
  scale_colour_manual(values=coloo_1)+
  facet_grid(vars(),vars(index),scales="free_x")+
  theme_bw() +
  theme(plot.title = element_text(size = 15,hjust = 0.5,vjust = 0.5),
        axis.title.y = element_blank()) +
  ylim(-ylim,ylim)

pyt_combine <- ggarrange(pyt2,pyt6,pyt7,pyt3,
                         ncol=2,nrow=2,
                         #labels = c('A','B','C','D','E','F'),
                         common.legend = TRUE,
                         legend="right")

# pyt_combine <- ggarrange(pyt2,pyt5,pyt6,pyt3,
#                          ncol=2,nrow=2,
#                          #labels = c('A','B','C','D','E','F'),
#                          common.legend = TRUE,
#                          legend="right")

#pyt_combine
ggsave(paste('gamma_',date_mark,'_',ver,'.png',sep = ''),plot = pyt1, bg="white", path = path_output_plots)
ggsave(paste('heterogeneity_',date_mark,'_',ver,'.png',sep = ''),plot = pyt2, bg="white", path = path_output_plots)
ggsave(paste('n_0_',date_mark,'_',ver,'.png',sep = ''),plot = pyt3, bg="white", path = path_output_plots)
ggsave(paste('inter_ratio_',date_mark,'_',ver,'.png',sep = ''),plot = pyt5, bg="white", path = path_output_plots)
ggsave(paste('inter_strength_',date_mark,'_',ver,'.png',sep = ''),plot = pyt7, bg="white", path = path_output_plots)
ggsave(paste('pZ_',date_mark,'_',ver,'.png',sep = ''),plot = pyt6, bg="white", path = path_output_plots)
ggsave(paste('parameters_',date_mark,'_',ver,'.png',sep = ''),plot = pyt_combine, bg="white", path = path_output_plots,width = 20,height = 14.144,dpi = 300)
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

color_vec_1 <- c(Genshinpalette::Genshinpalette('KAMISATOAYAKA',4)
                 ,Genshinpalette::Genshinpalette('SANGONOMIIYAKOKOMI',4)
                 ,Genshinpalette::Genshinpalette('HUTAO',(len_methods- len_other_methods))
)
color_vec_2 <- c(RColorBrewer::brewer.pal(8,"Pastel2"),RColorBrewer::brewer.pal((len_methods- len_other_methods),"Accent"))
color_vec_3 <- c(RColorBrewer::brewer.pal(8,"Pastel2"),Genshinpalette::Genshinpalette('HUTAO',(len_methods- len_other_methods+1))[-1])
color_vec_4 <- c(RColorBrewer::brewer.pal(8,"Pastel2"),RColorBrewer::brewer.pal((len_methods- len_other_methods),"Pastel1"))
color_vec_5 <- c(RColorBrewer::brewer.pal(8,"Pastel2"),RColorBrewer::brewer.pal((len_methods- len_other_methods),"Set3"))

color_vec_6 <- c(RColorBrewer::brewer.pal(8,"Set3"),Genshinpalette::Genshinpalette('HUTAO',(len_methods- len_other_methods+1))[-1])
color_vec_7 <- c(RColorBrewer::brewer.pal(8,"Set3"),Genshinpalette::Genshinpalette('HUTAO',(len_methods- len_other_methods)))

color_vec_8 <- c(RColorBrewer::brewer.pal(8,"Set2"),Genshinpalette::Genshinpalette('HUTAO',(len_methods- len_other_methods+1))[-1])
color_vec_9 <- c(RColorBrewer::brewer.pal(8,"Set2"),Genshinpalette::Genshinpalette('HUTAO',(len_methods- len_other_methods)))

color_vec_10 <- c(RColorBrewer::brewer.pal(8,"Set1"),Genshinpalette::Genshinpalette('HUTAO',(len_methods- len_other_methods+1))[-1])
color_vec_11 <- c(RColorBrewer::brewer.pal(8,"Set1"),Genshinpalette::Genshinpalette('HUTAO',(len_methods- len_other_methods)))

color_vec_12 <- c(RColorBrewer::brewer.pal(8,"Accent"),Genshinpalette::Genshinpalette('HUTAO',(len_methods- len_other_methods+1))[-1])
color_vec_13 <- c(RColorBrewer::brewer.pal(8,"Accent"),Genshinpalette::Genshinpalette('HUTAO',(len_methods- len_other_methods)))

color_vec_14 <- c(RColorBrewer::brewer.pal(8,"Dark2"),Genshinpalette::Genshinpalette('HUTAO',(len_methods- len_other_methods+1))[-1])
color_vec_15 <- c(RColorBrewer::brewer.pal(8,"Dark2"),Genshinpalette::Genshinpalette('HUTAO',(len_methods- len_other_methods)))


color_setting <- color_vec_9
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
#     scale_y_continuous(limits = c(kedu_down-0.1, kedu_up+0.1), 
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
zihao1 <- 15 #x轴标题字号
zihao2 <- 12 #x轴刻度字号
zihao3 <- 18 #图例文本字号
#-----------
data_gamma_beta <- save_df[save_df$heterogeneity == heter_level_list[k_heter] & save_df$n_0 == n_0_list[k_n0] & save_df$inter_ratio == ratio_inter_list[k_inter] & save_df$p_Z == p_Z_list[k_pz] & save_df$inter_strength == inter_strength_list[k_inter_strength] & save_df$alpha_modifier == alpha_list[k_alpha] & save_df$eta_modifier == eta_list[k_eta], ]
data_heterogeneity_beta <- save_df[save_df$gamma == gamma_list[k_gamma] & save_df$n_0 == n_0_list[k_n0] & save_df$inter_ratio == ratio_inter_list[k_inter] & save_df$p_Z == p_Z_list[k_pz] & save_df$inter_strength == inter_strength_list[k_inter_strength] & save_df$alpha_modifier == alpha_list[k_alpha] & save_df$eta_modifier == eta_list[k_eta], ]
data_n0_beta <- save_df[save_df$heterogeneity == heter_level_list[k_heter] & save_df$gamma == gamma_list[k_gamma] & save_df$inter_ratio == ratio_inter_list[k_inter] & save_df$p_Z == p_Z_list[k_pz] & save_df$inter_strength == inter_strength_list[k_inter_strength] & save_df$alpha_modifier == alpha_list[k_alpha] & save_df$eta_modifier == eta_list[k_eta], ]
data_inter_ratio_beta <- save_df[save_df$heterogeneity == heter_level_list[k_heter] & save_df$gamma == gamma_list[k_gamma]& save_df$n_0 == n_0_list[k_n0] & save_df$p_Z == p_Z_list[k_pz] & save_df$inter_strength == inter_strength_list[k_inter_strength] & save_df$alpha_modifier == alpha_list[k_alpha] & save_df$eta_modifier == eta_list[k_eta], ]
data_pz_beta <- save_df[save_df$heterogeneity == heter_level_list[k_heter] & save_df$gamma == gamma_list[k_gamma]& save_df$n_0 == n_0_list[k_n0] & save_df$inter_ratio == ratio_inter_list[k_inter] & save_df$inter_strength == inter_strength_list[k_inter_strength] & save_df$alpha_modifier == alpha_list[k_alpha] & save_df$eta_modifier == eta_list[k_eta], ]
data_inter_strength_beta <- save_df[save_df$heterogeneity == heter_level_list[k_heter] & save_df$gamma == gamma_list[k_gamma]& save_df$n_0 == n_0_list[k_n0] & save_df$inter_ratio == ratio_inter_list[k_inter] & save_df$p_Z == p_Z_list[k_pz] & save_df$alpha_modifier == alpha_list[k_alpha2] & save_df$eta_modifier == eta_list[k_eta], ]
data_alpha_beta <- save_df[save_df$heterogeneity == heter_level_list[k_heter] & save_df$gamma == gamma_list[k_gamma]& save_df$n_0 == n_0_list[k_n0] & save_df$inter_ratio == ratio_inter_list[k_inter] & save_df$p_Z == p_Z_list[k_pz]& save_df$inter_strength == inter_strength_list[k_inter_strength]& save_df$eta_modifier == eta_list[k_eta], ]
data_eta_beta <- save_df[save_df$heterogeneity == heter_level_list[k_heter] & save_df$gamma == gamma_list[k_gamma]& save_df$n_0 == n_0_list[k_n0] & save_df$inter_ratio == ratio_inter_list[k_inter] & save_df$p_Z == p_Z_list[k_pz]& save_df$inter_strength == inter_strength_list[k_inter_strength]& save_df$alpha_modifier == alpha_list[k_alpha], ]
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


kedu_up <- 1.8
kedu_down <- 0.2
kedu_fen <- 0.2
#------------------------------
P1<-ggplot(data=data_heterogeneity_beta1,aes(x = heterogeneity, y = beta_hat)) +
  geom_boxplot(aes(fill = method),outlier.alpha=0.5) +
  scale_fill_manual(values=color_setting)+
  geom_hline(aes(yintercept = 1), linetype = "dashed", color = "darkred") +
  theme_bw() +
  scale_y_continuous(limits = c(kedu_down-0.1, kedu_up+0.1), 
                     breaks = seq(kedu_down,kedu_up,kedu_fen))+
  xlab("heterogeneity level between the two samples") +
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
  geom_hline(aes(yintercept = 1), linetype = "dashed", color = "darkred") +
  theme_bw() +
  scale_y_continuous(limits = c(kedu_down-0.1, kedu_up+0.1), 
                     breaks = seq(kedu_down,kedu_up,kedu_fen))+
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

P3<-ggplot(data=data_gamma_beta1,aes(x = gamma, y = beta_hat)) +
  geom_boxplot(aes(fill = method),outlier.alpha=0.5) +
  scale_fill_manual(values=color_setting)+
  geom_hline(aes(yintercept = 1), linetype = "dashed", color = "darkred") +
  theme_bw() +
  scale_y_continuous(limits = c(kedu_down-0.1, kedu_up+0.1), 
                     breaks = seq(kedu_down,kedu_up,kedu_fen))+
  xlab("Relative magnitude of horizontal pleiotropy effect") +
  #ylab("Estimation") +
  ylab("") +
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

P4<-ggplot(data=data_n0_beta1,aes(x = n_0, y = beta_hat)) +
  geom_boxplot(aes(fill = method),outlier.alpha=0.5) +
  scale_fill_manual(values=color_setting)+
  geom_hline(aes(yintercept = 1), linetype = "dashed", color = "darkred") +
  theme_bw() +
  scale_y_continuous(limits = c(kedu_down-0.1, kedu_up+0.1), 
                     breaks = seq(kedu_down,kedu_up,kedu_fen))+
  xlab("Sample Size") +
  ylab("") +
  #ylab("Estimation") +
  theme(axis.title.y = element_blank(),  # 去掉 y 轴标题
        axis.text.y = element_blank(),   # 去掉 y 轴文本
        axis.ticks.y = element_blank()) + # 去掉 y 轴刻度
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
  geom_hline(aes(yintercept = 1), linetype = "dashed", color = "darkred") +
  theme_bw() +
  scale_y_continuous(limits = c(kedu_down-0.1, kedu_up+0.1), 
                     breaks = seq(kedu_down,kedu_up,kedu_fen))+
  xlab("Relative magnitude of interaction effects") +
  ylab("") +
  #ylab("") +
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

P6<-ggplot(data=data_alpha_beta1,aes(x = alpha_modifier, y = beta_hat)) +
  geom_boxplot(aes(fill = method),outlier.alpha=0.5) +
  scale_fill_manual(values=color_setting)+
  geom_hline(aes(yintercept = 1), linetype = "dashed", color = "darkred") +
  theme_bw() +
  scale_y_continuous(limits = c(kedu_down-0.1, kedu_up+0.1), 
                     breaks = seq(kedu_down,kedu_up,kedu_fen))+
  xlab("Relative magnitude of the effect of Z on X") +
  ylab("Estimation") +
  #ylab("") +
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
  geom_hline(aes(yintercept = 1), linetype = "dashed", color = "darkred") +
  theme_bw() +
  scale_y_continuous(limits = c(kedu_down-0.1, kedu_up+0.1), 
                     breaks = seq(kedu_down,kedu_up,kedu_fen))+
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


pyt_combine2 <- ggarrange(P1,P2,P4,P6,P5,P3,
                          ncol=3,nrow=2,
                          #labels = c('A','B','C','D','E','F'),
                          common.legend = TRUE,
                          legend="right")
pyt_combine2

#--------------
# data_inter_strength_beta1<-data_inter_strength_beta1[data_inter_strength_beta1$inter_strength != 0,]
# 
# P51<-ggplot(data=data_inter_strength_beta1,aes(x = inter_strength, y = beta_hat)) +
#   geom_boxplot(aes(fill = method),outlier.alpha=0.5) +
#   scale_fill_manual(values=color_setting)+
#   geom_hline(aes(yintercept = 1), linetype = "dashed", color = "darkred") +
#   theme_bw() +
#   scale_y_continuous(limits = c(kedu_down-0.1, kedu_up+0.1), 
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
#                           #labels = c('A','B','C','D','E','F'),
#                           common.legend = TRUE,
#                           legend="right")
# print(pyt_combine2)
#--------------
ggsave(paste('gamma_2_',date_mark,'_',ver,'.png',sep = ''),plot = P3, bg="white", path = path_output_plots)
ggsave(paste('heterogeneity_2_',date_mark,'_',ver,'.png',sep = ''),plot = P1, bg="white", path = path_output_plots)
ggsave(paste('n_0_2_',date_mark,'_',ver,'.png',sep = ''),plot = P4, bg="white", path = path_output_plots)
#ggsave(paste('inter_ratio_2_',date_mark,'_',ver,'.png',sep = ''),plot = pyt5, bg="white", path = path_output_plots)
ggsave(paste('inter_strength_2_',date_mark,'_',ver,'.png',sep = ''),plot = P5, bg="white", path = path_output_plots)
ggsave(paste('pZ_2_',date_mark,'_',ver,'.png',sep = ''),plot = P2, bg="white", path = path_output_plots)
ggsave(paste('alpha_2_',date_mark,'_',ver,'.png',sep = ''),plot = P6, bg="white", path = path_output_plots)
ggsave(paste('eta_2_',date_mark,'_',ver,'.png',sep = ''),plot = P7, bg="white", path = path_output_plots)
ggsave(paste('parameters_2_',date_mark,'_',ver,'.png',sep = ''),plot = pyt_combine2, bg="white", path = path_output_plots,width = 20,height = 14.144,dpi = 300)
ggsave(paste('parameters_2_',date_mark,'_',ver,'.pdf',sep = ''),plot = pyt_combine2, bg="white", path = path_output_plots,width = 20,height = 14.144,dpi = 300)
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
power_cover_plot <- function(data_test,color_setting,x_exist=TRUE,y_exist=TRUE,x_label){  
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
      theme(axis.text= element_text(size = 0),
            axis.title.x = element_blank(),  # 移除 x 轴标题
            axis.text.x = element_blank(),   # 移除 x 轴刻度
            axis.ticks.x = element_blank(),  # 移除 x 轴刻度线
            axis.title.y=element_text(size=19),
            plot.title=element_text(size=23)) +
      labs(y = 'Power') +
      theme(legend.title = element_text(size = 15),
            legend.text=element_text(size=15))+  
      geom_text(mapping = aes(label = power),vjust = -0.5)
    cover_plot_v2 <- ggplot(data_test,aes(x=method,y=cover_ratio,fill=method))+
      geom_bar(stat="identity",width=0.5)+
      theme_bw()+
      theme(panel.background = element_blank()) +
      scale_fill_manual(values=color_setting)+
      scale_y_continuous(limits = c(0, 1), 
                         breaks = c(0,0.2,0.4,0.6,0.8,1))+
      #theme(panel.grid =element_blank()) + 
      theme(axis.text= element_text(size = 0),
            #axis.title.x = element_blank(),  # 移除 x 轴标题
            axis.title.x=element_text(size=19),
            axis.text.x = element_blank(),   # 移除 x 轴刻度
            axis.ticks.x = element_blank(),  # 移除 x 轴刻度线
            axis.title.y=element_text(size=19),
            plot.title=element_text(size=23)) +
      labs(y = 'Cover Ratio',
           x = x_label) +
      theme(legend.title = element_text(size = 15),
            legend.text=element_text(size=15))+  
      geom_text(mapping = aes(label = cover_ratio),vjust = -0.5)
    
    pyt_combine3 <- ggarrange(power_plot_v2,cover_plot_v2,
                              ncol=1,nrow=2,
                              #labels = c('A','B','C','D','E','F'),
                              common.legend = TRUE,
                              legend="none")
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
      theme(axis.text= element_text(size = 0),
            axis.title.x = element_blank(),  # 移除 x 轴标题
            axis.text.x = element_blank(),   # 移除 x 轴刻度
            axis.ticks.x = element_blank(),  # 移除 x 轴刻度线
            
            axis.title.y = element_blank(),  # 移除 y 轴标题
            axis.text.y = element_blank(),   # 移除 y 轴刻度
            axis.ticks.y = element_blank(),  # 移除 y 轴刻度线
            
            #axis.title.y=element_text(size=19),
            plot.title=element_text(size=23)) +
      labs(y = '') +
      theme(legend.title = element_text(size = 15),
            legend.text=element_text(size=15))+  
      geom_text(mapping = aes(label = power),vjust = -0.5)
    
    cover_plot_v2 <- ggplot(data_test,aes(x=method,y=cover_ratio,fill=method))+
      geom_bar(stat="identity",width=0.5)+
      theme_bw()+
      theme(panel.background = element_blank()) +
      scale_fill_manual(values=color_setting)+
      scale_y_continuous(limits = c(0, 1), 
                         breaks = c(0,0.2,0.4,0.6,0.8,1))+
      #theme(panel.grid =element_blank()) + 
      theme(axis.text= element_text(size = 0),
            #axis.title.x = element_blank(),  # 移除 x 轴标题
            axis.title.x=element_text(size=19),
            axis.text.x = element_blank(),   # 移除 x 轴刻度
            axis.ticks.x = element_blank(),  # 移除 x 轴刻度线
            
            axis.title.y = element_blank(),  # 移除 y 轴标题
            axis.text.y = element_blank(),   # 移除 y 轴刻度
            axis.ticks.y = element_blank(),  # 移除 y 轴刻度线
            
            #axis.title.y=element_text(size=19),
            plot.title=element_text(size=23)) +
      labs(y = '',
           x = x_label) +
      theme(legend.title = element_text(size = 15),
            legend.text=element_text(size=15))+  
      geom_text(mapping = aes(label = cover_ratio),vjust = -0.5)
    
    pyt_combine3 <- ggarrange(power_plot_v2,cover_plot_v2,
                              ncol=1,nrow=2,
                              #labels = c('A','B','C','D','E','F'),
                              common.legend = TRUE,
                              legend="none")
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
      theme(axis.text= element_text(size = 0),
            axis.title.x = element_blank(),  # 移除 x 轴标题
            axis.text.x = element_blank(),   # 移除 x 轴刻度
            axis.ticks.x = element_blank(),  # 移除 x 轴刻度线
            
            axis.title.y=element_text(size=19),
            plot.title=element_text(size=23)) +
      labs(y = 'Power') +
      theme(legend.title = element_text(size = 15),
            legend.text=element_text(size=15))+  
      geom_text(mapping = aes(label = power),vjust = -0.5)
    
    cover_plot_v2 <- ggplot(data_test,aes(x=method,y=cover_ratio,fill=method))+
      geom_bar(stat="identity",width=0.5)+
      theme_bw()+
      theme(panel.background = element_blank()) +
      scale_fill_manual(values=color_setting)+
      scale_y_continuous(limits = c(0, 1), 
                         breaks = c(0,0.2,0.4,0.6,0.8,1))+
      #theme(panel.grid =element_blank()) + 
      theme(axis.text= element_text(size = 0),
            #axis.title.x = element_blank(),  # 移除 x 轴标题
            axis.title.x=element_text(size=19),
            axis.text.x = element_blank(),   # 移除 x 轴刻度
            axis.ticks.x = element_blank(),  # 移除 x 轴刻度线
            
            axis.title.y=element_text(size=19),
            plot.title=element_text(size=23)) +
      labs(y = 'Cover Ratio',
           x = x_label) +
      theme(legend.title = element_text(size = 15),
            legend.text=element_text(size=15))+  
      geom_text(mapping = aes(label = cover_ratio),vjust = -0.5)
    
    pyt_combine3 <- ggarrange(power_plot_v2,cover_plot_v2,
                              ncol=1,nrow=2,
                              #labels = c('A','B','C','D','E','F'),
                              common.legend = TRUE,
                              legend="none")
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
      theme(axis.text= element_text(size = 0),
            axis.title.x = element_blank(),  # 移除 x 轴标题
            axis.text.x = element_blank(),   # 移除 x 轴刻度
            axis.ticks.x = element_blank(),  # 移除 x 轴刻度线
            
            axis.title.y = element_blank(),  # 移除 y 轴标题
            axis.text.y = element_blank(),   # 移除 y 轴刻度
            axis.ticks.y = element_blank(),  # 移除 y 轴刻度线
            
            #axis.title.y=element_text(size=19),
            plot.title=element_text(size=23)) +
      labs(y = '') +
      theme(legend.title = element_text(size = 15),
            legend.text=element_text(size=15))+  
      geom_text(mapping = aes(label = power),vjust = -0.5)
    
    cover_plot_v2 <- ggplot(data_test,aes(x=method,y=cover_ratio,fill=method))+
      geom_bar(stat="identity",width=0.5)+
      theme_bw()+
      theme(panel.background = element_blank()) +
      scale_fill_manual(values=color_setting)+
      scale_y_continuous(limits = c(0, 1), 
                         breaks = c(0,0.2,0.4,0.6,0.8,1))+
      #theme(panel.grid =element_blank()) + 
      theme(axis.text= element_text(size = 0),
            #axis.title.x = element_blank(),  # 移除 x 轴标题
            axis.title.x=element_text(size=19),
            axis.text.x = element_blank(),   # 移除 x 轴刻度
            axis.ticks.x = element_blank(),  # 移除 x 轴刻度线
            
            axis.title.y = element_blank(),  # 移除 y 轴标题
            axis.text.y = element_blank(),   # 移除 y 轴刻度
            axis.ticks.y = element_blank(),  # 移除 y 轴刻度线
            
            #axis.title.y=element_text(size=19),
            plot.title=element_text(size=23)) +
      labs(y = '',
           x = x_label) +
      theme(legend.title = element_text(size = 15),
            legend.text=element_text(size=15))+  
      geom_text(mapping = aes(label = cover_ratio),vjust = -0.5)
    
    pyt_combine3 <- ggarrange(power_plot_v2,cover_plot_v2,
                              ncol=1,nrow=2,
                              #labels = c('A','B','C','D','E','F'),
                              common.legend = TRUE,
                              legend="none")
  }
  return(pyt_combine3)
}

power_plot <- function(data,color_setting,x_exist=TRUE,y_exist=TRUE,x_label,y_label,axis){  
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
      geom_text(mapping = aes(label = power),vjust = -0.5)
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
      geom_text(mapping = aes(label = power),vjust = -0.5)
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
      geom_text(mapping = aes(label = power),vjust = -0.5)
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
      geom_text(mapping = aes(label = power),vjust = -0.5)
  }
  return(power_plot_v2)
}
#--------------------------------------------------------
data_power_list <- vector('list',length = length(save_power_list))  #最终画图数据
for(i_p in 1:length(save_power_list)){
  test_df <- save_power_list[[i_p]]
  data_test <- data.frame(method=method_names,power=rep(NA,len_methods),cover_ratio=rep(NA,len_methods),n_0=rep(n_0_list[i_p],len_methods))
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
write.csv(data_power_df,file=paste(path_output_tables,'/data_power_',date_mark,'_',ver,'.csv',sep = ''))
#------------------------------------------------------
#data,color_setting,x_exist=TRUE,y_exist=TRUE,x_label,y_label,axis
power_plot_1 <- power_plot(data=data_power_list[[1]]
                           ,color_setting = color_setting
                           ,x_exist = FALSE
                           ,y_exist = TRUE
                           ,x_label = paste0('Sample size = ',n_0_list[1])
                           ,y_label = 'Power'
                           ,axis = 2
)
power_plot_2 <- power_plot(data=data_power_list[[2]]
                           ,color_setting = color_setting
                           ,x_exist = FALSE
                           ,y_exist = FALSE
                           ,x_label = paste0('Sample size = ',n_0_list[2])
                           ,y_label = 'Power'
                           ,axis = 2
)
power_plot_3 <- power_plot(data=data_power_list[[3]]
                           ,color_setting = color_setting
                           ,x_exist = FALSE
                           ,y_exist = FALSE
                           ,x_label = paste0('Sample size = ',n_0_list[3])
                           ,y_label = 'Power'
                           ,axis = 2
)
power_plot_4 <- power_plot(data=data_power_list[[4]]
                           ,color_setting = color_setting
                           ,x_exist = FALSE
                           ,y_exist = TRUE
                           ,x_label = paste0('Sample size = ',n_0_list[4])
                           ,y_label = 'Power'
                           ,axis = 2
)
power_plot_5 <- power_plot(data=data_power_list[[5]]
                           ,color_setting = color_setting,
                           x_exist = FALSE
                           ,y_exist = FALSE
                           ,x_label = paste0('Sample size = ',n_0_list[5])
                           ,y_label = 'Power'
                           ,axis = 2
)
power_plot_6 <- power_plot(data=data_power_list[[6]]
                           ,color_setting = color_setting
                           ,x_exist = FALSE
                           ,y_exist = FALSE
                           ,x_label = paste0('Sample size = ',n_0_list[6])
                           ,y_label = 'Power'
                           ,axis = 2
)

power_plot_all <- ggarrange(power_plot_1,power_plot_2
                            ,power_plot_3,power_plot_4
                            ,power_plot_5,power_plot_6
                            ,ncol=3,nrow=2,
                            common.legend = TRUE,
                            legend="right")
#print(power_plot_all)
#---------------------
cover_plot_1 <- power_plot(data=data_power_list[[1]]
                           ,color_setting = color_setting
                           ,x_exist = FALSE
                           ,y_exist = TRUE
                           ,x_label = paste0('Sample size = ',n_0_list[1])
                           ,y_label = 'Cover Ratio'
                           ,axis = 3
)
cover_plot_2 <- power_plot(data=data_power_list[[2]]
                           ,color_setting = color_setting
                           ,x_exist = FALSE
                           ,y_exist = FALSE
                           ,x_label = paste0('Sample size = ',n_0_list[2])
                           ,y_label = 'Cover Ratio'
                           ,axis = 3
)
cover_plot_3 <- power_plot(data=data_power_list[[3]]
                           ,color_setting = color_setting
                           ,x_exist = FALSE
                           ,y_exist = FALSE
                           ,x_label = paste0('Sample size = ',n_0_list[3])
                           ,y_label = 'Cover Ratio'
                           ,axis = 3
)
cover_plot_4 <- power_plot(data=data_power_list[[4]]
                           ,color_setting = color_setting
                           ,x_exist = FALSE
                           ,y_exist = TRUE
                           ,x_label = paste0('Sample size = ',n_0_list[4])
                           ,y_label = 'Cover Ratio'
                           ,axis = 3
)
cover_plot_5 <- power_plot(data=data_power_list[[5]]
                           ,color_setting = color_setting,
                           x_exist = FALSE
                           ,y_exist = FALSE
                           ,x_label = paste0('Sample size = ',n_0_list[5])
                           ,y_label = 'Cover Ratio'
                           ,axis = 3
)
cover_plot_6 <- power_plot(data=data_power_list[[6]]
                           ,color_setting = color_setting
                           ,x_exist = FALSE
                           ,y_exist = FALSE
                           ,x_label = paste0('Sample size = ',n_0_list[6])
                           ,y_label = 'Cover Ratio'
                           ,axis = 3
)

cover_plot_all <- ggarrange(cover_plot_1,cover_plot_2
                            ,cover_plot_3,cover_plot_4
                            ,cover_plot_5,cover_plot_6
                            ,ncol=3,nrow=2,
                            # labels = c('A','B','C','D','E','F'),
                            common.legend = TRUE,
                            legend="right")

#-----------------------------------------------------
#组合图形
#########
power_cover_plot1 <- power_cover_plot(data_test = data_power_list[[1]]
                                      ,color_setting = color_setting
                                      ,x_exist = FALSE
                                      ,y_exist = TRUE
                                      ,x_label = paste0('Sample size = ',n_0_list[1])
)
power_cover_plot2 <- power_cover_plot(data_test = data_power_list[[2]]
                                      ,color_setting = color_setting
                                      ,x_exist = FALSE
                                      ,y_exist = FALSE
                                      ,x_label = paste0('Sample size = ',n_0_list[2])
)
power_cover_plot3 <- power_cover_plot(data_test = data_power_list[[3]]
                                      ,color_setting = color_setting
                                      ,x_exist = FALSE
                                      ,y_exist = FALSE
                                      ,x_label = paste0('Sample size = ',n_0_list[3])
)
power_cover_plot4 <- power_cover_plot(data_test = data_power_list[[4]]
                                      ,color_setting = color_setting
                                      ,x_exist = TRUE
                                      ,y_exist = TRUE
                                      ,x_label = paste0('Sample size = ',n_0_list[4])
)
power_cover_plot5 <- power_cover_plot(data_test = data_power_list[[5]]
                                      ,color_setting = color_setting
                                      ,x_exist = TRUE
                                      ,y_exist = FALSE
                                      ,x_label = paste0('Sample size = ',n_0_list[5])
)
power_cover_plot6 <- power_cover_plot(data_test = data_power_list[[6]]
                                      ,color_setting = color_setting
                                      ,x_exist = TRUE
                                      ,y_exist = FALSE
                                      ,x_label = paste0('Sample size = ',n_0_list[6])
)

# 创建图例的空白图
p_with_legend <- ggplot(data_power_list[[1]], aes(x=method,y=power,fill=method) ) +
  geom_bar(stat="identity",width=0.5) +
  scale_fill_manual(values=color_setting)+
  theme(legend.title = element_text(size = 15),
        legend.text=element_text(size=15))+ 
  theme_void() +  # 不显示任何图形内容
  theme(legend.position = "right")  # 只保留图例
legend_only <- cowplot::get_plot_component(p_with_legend, 'guide-box-right', return_all = TRUE)
#组合图形
power_cover_combined_plot <- ggarrange(power_cover_plot1,power_cover_plot2
                                       ,power_cover_plot3,power_cover_plot4
                                       ,power_cover_plot5,power_cover_plot6
                                       ,ncol=3,nrow=2,
                                       common.legend = TRUE,
                                       legend="none")
power_cover_plot_with_legend <- cowplot::plot_grid(
  power_cover_combined_plot,  # 主图
  legend_only,  # 单独的图例
  ncol = 2,  # 两列布局
  rel_widths = c(5, 0.5)  # 调整图形与图例的宽度比例
)
#power_cover_plot_with_legend
#------------------------- 
ggsave(paste('power_',date_mark,'_',ver,'.pdf',sep = ''),plot = power_plot_all, bg="white", path = path_output_plots, height = 9,width = 20,dpi = 300)
ggsave(paste('cover_ratio_',date_mark,'_',ver,'.pdf',sep = ''),plot = cover_plot_all, bg="white", path = path_output_plots, height = 9,width = 20,dpi = 300)
ggsave(paste('power+cover_',date_mark,'_',ver,'.pdf',sep = ''),plot = power_cover_plot_with_legend, bg="white", path = path_output_plots, height = 15,width = 20,dpi = 300)

ggsave(paste('power_',date_mark,'_',ver,'.png',sep = ''),plot = power_plot_all, bg="white", path = path_output_plots, height = 9,width = 20,dpi = 300)
ggsave(paste('cover_ratio_',date_mark,'_',ver,'.png',sep = ''),plot = cover_plot_all, bg="white", path = path_output_plots, height = 9,width = 20,dpi = 300)
ggsave(paste('power+cover_',date_mark,'_',ver,'.png',sep = ''),plot = power_cover_plot_with_legend, bg="white", path = path_output_plots, height = 15,width = 20,dpi = 300)
#-----------------------------------------------