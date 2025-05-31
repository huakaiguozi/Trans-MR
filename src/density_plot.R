library(ggplot2)

# 合并数据并添加分组标签
df <- rbind(
  data.frame(value = X_0, Population = "Southem UK population"),
  data.frame(value = X_1, Population = "Northem UK population")
)

# 绘制密度图
p1 <- ggplot(df, aes(x = value, fill = group)) +
  geom_density(alpha = 0.5) +
  labs(title = "Density plot of BMI across two samples",
       x = "BMI", y = "Density") +
  theme_minimal()
p1 <- ggplot(df, aes(x = value, fill = Population)) +
  geom_density(alpha = 0.5) +
  labs(title = "Density plot of BMI across two samples",
       x = "BMI", y = "Density") +
  theme_minimal(base_size = 16) +  # 设置基础字体大小
  theme(
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5), # 图标题居中加粗
    axis.title = element_text(size = 18),  # 坐标轴标题
    axis.text = element_text(size = 16),   # 坐标轴刻度
    legend.title = element_text(size = 16),# 图例标题
    legend.text = element_text(size = 15)  # 图例内容
  )

ggsave(paste('density_',date_mark,'_',ver,'.pdf',sep = ''),plot = p1, bg="white", path = path_output_plots, height = 13.5,width = 18,dpi = 500)
ggsave(paste('density',date_mark,'_',ver,'.png',sep = ''),plot = p1, bg="white", path = path_output_plots, height = 13.5,width = 18,dpi = 500)



# 先确保提取的是数据框，且列名是字符
df_0 <- as.data.frame(C_0_origin_save[, c(11,12)])
colnames(df_0) <- c("腰围", "臀围")
df_0$group <- "C_0"

df_1 <- as.data.frame(C_1_origin_save[, c(11,12)])
colnames(df_1) <- c("腰围", "臀围")
df_1$group <- "C_1"

# 合并数据并转换为长格式
df_all <- bind_rows(df_0, df_1) %>%
  pivot_longer(cols = c("腰围", "臀围"), names_to = "变量", values_to = "数值")

# 绘制密度图
ggplot(df_all, aes(x = 数值, fill = group)) +
  geom_density(alpha = 0.5) +
  facet_wrap(~变量, scales = "free") +
  labs(title = "C_0 与 C_1 的腰围与臀围密度对比", x = "数值", y = "密度") +
  theme_minimal()
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# 假设 data_result_save_6 已经加载
df <- data_result_save_6 %>%
  filter(sample == 0) %>%
  mutate(
    method = factor(method, levels = unique(method)),
    outcome = factor(outcome, levels = c("FEV1", "FEV1/FVC", "FVC")),
    fill_color = case_when(
      d1 == 1 ~ "palegreen",
      d1 == 0 ~ "lightcoral",
      TRUE ~ "grey80"
    ),
    label_text = paste0(round(beta, 3), ifelse(p_star == "-", "", p_star)),
    triangle = ifelse(outcome == "FEV1/FVC", ifelse(d2 == 1, "▲", "▼"), NA)
  )

ggplot(df, aes(x = method, y = outcome)) +
  geom_tile(aes(fill = fill_color), color = "white") +
  scale_fill_identity() +
  geom_text(aes(label = label_text), size = 3.5, color = "black") +
  geom_text(
    data = df %>% filter(!is.na(triangle)),
    aes(label = triangle),
    size = 5,
    color = "black",
    nudge_x = 0.3,
    nudge_y = -0.3
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  ) +
  labs(x = "方法", y = "表型", title = "方法与表型的效应图")
