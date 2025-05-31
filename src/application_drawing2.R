C_0_draw <- data.table(C_0_origin)
C_1_draw <- data.table(C_1_origin)

C_0_draw$group <- 'sample0'
C_1_draw$group <- 'sample1'

C_draw <- rbind(C_0_draw,C_1_draw)


C_draw_long <-  C_draw %>%
  pivot_longer(cols = -group, names_to = "Variable", values_to = "Value")


ggplot(C_draw_long, aes(x = Variable, y = Value, fill = group)) +
  geom_violin(trim = FALSE, alpha = 0.5, position = position_dodge(0.8), size = 1) +  # 小提琴图
  geom_boxplot(width = 0.15, position = position_dodge(0.8), outlier.shape = NA, alpha = 0.8) +  # 叠加箱线图
  stat_summary(fun = mean, geom = "point", position = position_dodge(0.8), 
               shape = 21, size = 3, fill = "white", color = "black", stroke = 1.2) +  # 显示均值点
  scale_fill_manual(values = c("#1f77b4", "#ff7f0e")) +  # 自定义填充颜色
  labs(title = "Violin Plot of Covariates", x = "Covariate", y = "Value") +
  theme_minimal(base_size = 16) +
  theme(legend.position = "top",  # 图例放在顶部
        legend.title = element_blank(),  # 去掉图例标题
        plot.title = element_text(hjust = 0.5, face = "bold"),  # 标题居中 & 加粗
        axis.text = element_text(size = 14),  # 坐标轴标签大小
        axis.title = element_text(size = 16)) 


C_0_draw1 <- C_0_draw[,-c(2,3)]
C_1_draw1 <- C_1_draw[,-c(2,3)]
C_draw1 <- rbind(C_0_draw1,C_1_draw1)

ggplot(C_draw1, aes(x = 腰围, fill = group, color = group)) +
  geom_density(alpha = 0.4, size = 1.2) +  # 透明度 & 线条粗细
  scale_fill_manual(values = c("#1f77b4", "#ff7f0e")) +  # 填充色（蓝 & 橙）
  scale_color_manual(values = c("#1f77b4", "#ff7f0e")) + # 线条颜色
  labs(title = "Density Plot of Covariate Distribution", 
       x = "Covariate Value", 
       y = "Density") +
  theme_minimal(base_size = 16) +  # 使用minimal主题，字体大小16
  theme(legend.position = "top",  # 图例放在上方
        legend.title = element_blank(),  # 去掉图例标题
        plot.title = element_text(hjust = 0.5, face = "bold"),  # 标题居中 & 加粗
        axis.text = element_text(size = 14),  # 坐标轴字体大小
        axis.title = element_text(size = 16)) 


C_0_draw2 <- C_0_draw[,-c(1,3)]
C_1_draw2 <- C_1_draw[,-c(1,3)]

C_0_draw2 <- C_0_draw2[C_0_draw2$C反应蛋白 <= quantile(C_0_draw2$C反应蛋白,0.95),]
C_1_draw2 <- C_1_draw2[C_1_draw2$C反应蛋白 <= quantile(C_1_draw2$C反应蛋白,0.95),]

C_draw2 <- rbind(C_0_draw2,C_1_draw2)

ggplot(C_draw2, aes(x = C反应蛋白, fill = group, color = group)) +
  geom_density(alpha = 0.4, size = 1.2) +  # 透明度 & 线条粗细
  scale_fill_manual(values = c("#1f77b4", "#ff7f0e")) +  # 填充色（蓝 & 橙）
  scale_color_manual(values = c("#1f77b4", "#ff7f0e")) + # 线条颜色
  labs(title = "Density Plot of Covariate Distribution", 
       x = "Covariate Value", 
       y = "Density") +
  theme_minimal(base_size = 16) +  # 使用minimal主题，字体大小16
  theme(legend.position = "top",  # 图例放在上方
        legend.title = element_blank(),  # 去掉图例标题
        plot.title = element_text(hjust = 0.5, face = "bold"),  # 标题居中 & 加粗
        axis.text = element_text(size = 14),  # 坐标轴字体大小
        axis.title = element_text(size = 16)) 


C_0_draw3 <- C_0_draw[,-c(1,2)]
C_1_draw3 <- C_1_draw[,-c(1,2)]
C_draw3 <- rbind(C_0_draw3,C_1_draw3)

ggplot(C_draw3, aes(x = 睡眠时长, fill = group, color = group)) +
  geom_density(alpha = 0.4, size = 1.2) +  # 透明度 & 线条粗细
  scale_fill_manual(values = c("#1f77b4", "#ff7f0e")) +  # 填充色（蓝 & 橙）
  scale_color_manual(values = c("#1f77b4", "#ff7f0e")) + # 线条颜色
  labs(title = "Density Plot of Covariate Distribution", 
       x = "Covariate Value", 
       y = "Density") +
  theme_minimal(base_size = 16) +  # 使用minimal主题，字体大小16
  theme(legend.position = "top",  # 图例放在上方
        legend.title = element_blank(),  # 去掉图例标题
        plot.title = element_text(hjust = 0.5, face = "bold"),  # 标题居中 & 加粗
        axis.text = element_text(size = 14),  # 坐标轴字体大小
        axis.title = element_text(size = 16)) 
