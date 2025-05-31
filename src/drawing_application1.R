dt_forest <- dt_forest %>% 
  mutate(empty_column = '                      ') %>%  # 创建空列
  #mutate(empty_column = NA) %>%  # 创建空列
  select(Exposure,Outcome,target,Method,empty_column,beta,betaCI2,sd,p,loci,upci) #%>%
  #rename(` ` = empty_column)  # 将列名改为 "空格"


target_vec <- dt_forest$target
target_vec[target_vec == 1] <- 'Sample 1'
target_vec[target_vec == 0] <- 'Sample 0'
dt_forest$target <- target_vec
# 格式化P值（添加星号标记）
dt_forest$p_formatted <- ifelse(dt_forest$p < 0.001, 
                                "<0.001*",
                                sprintf("%.3f%s", dt_forest$p, 
                                        ifelse(dt_forest$p < 0.05, "*", " ")))

# 生成置信区间文本
dt_forest$CI_text <- sprintf("[%.3f, %.3f]", dt_forest$loci, dt_forest$upci)



# 创建空行模板（与原数据框结构一致）
empty_row <- dt_forest[NA, ][1, ]

# 按两行一组拆分数据框
groups <- split(dt_forest, (seq_len(nrow(dt_forest)) - 1) %/% 2)

# 在每组后添加空行（最后一组不添加）
new_groups <- lapply(seq_along(groups), function(i) {
  if (i < length(groups)) {
    rbind(groups[[i]], empty_row)
  } else {
    groups[[i]]
  }
})

# 合并所有组并重置行名
new_dt <- do.call(rbind, new_groups)
rownames(new_dt) <- NULL
#=====================================================================
Methods1 <- as.character(new_dt$Method)
# 创建一个逻辑向量，标记每个元素是否是该元素第一次出现
is_first <- !duplicated(Methods1)
# 创建新的向量，保留第一次出现的元素，其他设为 NA
Methods1_new <- ifelse(is_first, Methods1, NA)

# 1. 从数据框中提取
Methods2 <- c('Method',NA,NA,Methods1_new)
Pop2 <- c('Target population',NA,NA,new_dt$target)
Empty2 <- c(NA,NA,NA,new_dt$empty_column)
Beta2   <- c('Beta',NA,NA,new_dt$beta)
Low2    <- c(NA,NA,NA,new_dt$loci)
Up2     <- c(NA,NA,NA,new_dt$upci)
CI2     <- c('95% CI',NA,NA,new_dt$CI_text)
Pval2   <- c('Pvalue',NA,NA,new_dt$p_formatted)
is_sum  <- c(TRUE,FALSE,FALSE,rep(FALSE,nrow(new_dt)))  # 如果没有汇总行，全部设为 FALSE

Beta2_2 <- as.numeric(Beta2)
Low2_2 <- as.numeric(Low2)
Up2_2 <- as.numeric(Up2)

# 2. 调用 forestplot
library(forestplot)

# （先定义那个交替上色的 fn，和你原来完全一样）
fn <- local({
  i <- 0
  function(..., clr.line, clr.marker){
    i <<- i + 1
    if(i%%2==0){
      fpDrawNormalCI(..., clr.line = "#000000", clr.marker = "#00B9BF")
    } else {
      fpDrawNormalCI(..., clr.line = "#000000", clr.marker = "#F9675C")
    }
  }
})
# label_gp <- vector("list", length(Beta2))
# for(i in seq_len(nrow(new_dt))){
#   if(i == 1){
#     # 第1行：3列都加粗
#     label_gp[[i]] <- rep(list(gpar(fontface="bold", cex=1.5)), ncol(new_dt))
#   } else {
#     # 其他行：正常体
#     label_gp[[i]] <- rep(list(gpar(fontface="plain", cex=1.5)), ncol(new_dt))
#   }
# }

Range1 <- c(Low2,Up2)
Upper1 <- quantile(c(Low2,Up2),0.97,na.rm = T)
Lower1 <- quantile(c(Low2,Up2),0.03,na.rm = T)
Range2 <- Range1[Range1 <= Upper1]
Range3 <- Range2[Range2 >= Lower1]
Upper2 <- max(pretty(Range3, n = 5))
Lower2 <- min(pretty(Range3, n = 5))

pdf(file = paste0(path_output_plots, '/forest_BMI_250519_2', date_mark, '_', ver, '.pdf'),
    height = 11.25, width = 15)
forestplot(
  #labeltext   = cbind(Methods2, Pop2, Empty2, Beta2, CI2, Pval2),
  labeltext   = cbind(Methods2, Pop2, Beta2, CI2, Pval2),
  mean        = Beta2_2,
  lower       = Low2_2,
  upper       = Up2_2,
  align       =c('l','l','l','c','c','c'),
  #xlab        = "Estimated regression coefficient for depression per 1‑unit increase in BMI",
  xlab        = 'Estimated change in Forced Vital Capacity(FVC) per 1-unit increase in BMI',
  #xlim        =c(-0.1,0.1),
  #xticks      = pretty(c(Low2, Up2), n = 5),
  xticks      = seq(-0.1,0.05,0.05),#pretty(Range3, n = 5),
  clip        = c(-0.1,0.05),#c(Lower2,Upper2),
  is.summary  = is_sum,#is.na(Beta2),#is_sum,
  #fn.ci_nrom  = 'fpDrawDiamondCI',
  hrzl_lines  = list(
    "3" = gpar(lwd=2, col="gray"),
    "5" = gpar(lwd=60, lineend="butt", columns=2:6, col="#99999920"),
    "8" = gpar(lwd=60, lineend="butt", columns=2:6, col="#99999920"),
    "11" = gpar(lwd=60, lineend="butt", columns=2:6, col="#99999920"),
    "14" = gpar(lwd=60, lineend="butt", columns=2:6, col="#99999920"),
    "17" = gpar(lwd=60, lineend="butt", columns=2:6, col="#99999920"),
    "20" = gpar(lwd=60, lineend="butt", columns=2:6, col="#99999920"),
    "23" = gpar(lwd=60, lineend="butt", columns=2:6, col="#99999920"),
    "26" = gpar(lwd=60, lineend="butt", columns=2:6, col="#99999920"),
    "29" = gpar(lwd=60, lineend="butt", columns=2:6, col="#99999920")
    # … 你可以根据行数继续添加
  ),
  #txt_gp=fpTxtGp(label=gpar(cex=2),ticks=gpar(cex=2),xlab=gpar(cex = 2), title=gpar(cex = 2)),
  txt_gp      = fpTxtGp(label=gpar(cex=1.5),
                        ticks=gpar(cex=1.5),
                        xlab =gpar(cex=1.5),
                        summary = gpar(fontface="bold") ),
  col         = fpColors(box='steelblue', lines="black", zero="gray50"),
  fn.ci_norm  = fn,
  zero        = 0,
  cex         = 1.5,
  lineheight  = "auto",
  #line.margin = unit(2,'mm'),
  boxsize     = 0.4,
  colgap      = unit(6,"mm"),
  lwd.ci      = 3,
  ci.vertices = TRUE,
  ci.vertices.height = 0.2,
  graph.pos   = 3,
  #graphwidth  = unit(12, 'cm'),
  new_page    = TRUE
)
dev.off()

png(filename = paste0(path_output_plots, '/forest_BMI_250519_2', date_mark, '_', ver, '.png'),
    height = 11.25, width = 15, units = "in", res = 500)
forestplot(
  #labeltext   = cbind(Methods2, Pop2, Empty2, Beta2, CI2, Pval2),
  labeltext   = cbind(Methods2, Pop2, Beta2, CI2, Pval2),
  mean        = Beta2_2,
  lower       = Low2_2,
  upper       = Up2_2,
  align       =c('l','l','l','c','c','c'),
  #xlab        = c("Estimated regression coefficient for depression per 1-unit increase in BMI"),
  xlab        = 'Estimated change in Forced Vital Capacity(FVC) per 1-unit increase in BMI',
  #xlim        =c(-0.1,0.1),
  #xticks      = pretty(c(Low2, Up2), n = 5),
  xticks      = seq(-0.1,0.05,0.05),#pretty(Range3, n = 5),
  clip        = c(-0.1,0.05),#c(Lower2,Upper2),
  is.summary  = is_sum,#is.na(Beta2),#is_sum,
  #fn.ci_nrom  = 'fpDrawDiamondCI',
  hrzl_lines  = list(
    "3" = gpar(lwd=2, col="gray"),
    "5" = gpar(lwd=60, lineend="butt", columns=2:6, col="#99999920"),
    "8" = gpar(lwd=60, lineend="butt", columns=2:6, col="#99999920"),
    "11" = gpar(lwd=60, lineend="butt", columns=2:6, col="#99999920"),
    "14" = gpar(lwd=60, lineend="butt", columns=2:6, col="#99999920"),
    "17" = gpar(lwd=60, lineend="butt", columns=2:6, col="#99999920"),
    "20" = gpar(lwd=60, lineend="butt", columns=2:6, col="#99999920"),
    "23" = gpar(lwd=60, lineend="butt", columns=2:6, col="#99999920"),
    "26" = gpar(lwd=60, lineend="butt", columns=2:6, col="#99999920"),
    "29" = gpar(lwd=60, lineend="butt", columns=2:6, col="#99999920")
    # … 你可以根据行数继续添加
  ),
  #txt_gp=fpTxtGp(label=gpar(cex=2),ticks=gpar(cex=2),xlab=gpar(cex = 2), title=gpar(cex = 2)),
  txt_gp      = fpTxtGp(label=gpar(cex=1.5),
                        ticks=gpar(cex=1.5),
                        xlab =gpar(cex=1.5),
                        summary = gpar(fontface="bold") ),
  col         = fpColors(box='steelblue', lines="black", zero="gray50"),
  fn.ci_norm  = fn,
  zero        = 0,
  cex         = 1.5,
  lineheight  = "auto",
  #line.margin = unit(2,'mm'),
  boxsize     = 0.4,
  colgap      = unit(6,"mm"),
  lwd.ci      = 3,
  ci.vertices = TRUE,
  ci.vertices.height = 0.2,
  graph.pos   = 3,
  #graphwidth  = unit(12, 'cm'),
  new_page    = TRUE
)
dev.off()

# 
# ggsave(paste('forest_BMI_depression',date_mark,'_',ver,'.pdf',sep = ''),plot = forest_plot1, bg="white", path = path_output_plots, height = 13.5,width = 18,dpi = 500)
# ggsave(paste('forest_BMI_depression',date_mark,'_',ver,'.png',sep = ''),plot = forest_plot1, bg="white", path = path_output_plots, height = 13.5,width = 18,dpi = 500)
