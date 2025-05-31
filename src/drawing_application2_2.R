# 2. 调用 forestplot
library(forestplot)


dt_forest <- dt_forest_COPD

dt_forest <- dt_forest %>% 
  mutate(empty_column = '                      ') %>%  # 创建空列
  #mutate(empty_column = NA) %>%  # 创建空列
  select(Exposure,Outcome,target,Method,empty_column,beta,betaCI2,sd,p,loci,upci) #%>%
#rename(` ` = empty_column)  # 将列名改为 "空格"


target_vec <- dt_forest$target
target_vec[target_vec == 1] <- 'Southern UK Population'
target_vec[target_vec == 0] <- 'Northern UK Population'
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
#=========================
Methods1 <- as.character(new_dt$Method)
# 创建一个逻辑向量，标记每个元素是否是该元素第一次出现
is_first <- !duplicated(Methods1)
# 创建新的向量，保留第一次出现的元素，其他设为 NA
Methods1_new <- ifelse(is_first, Methods1, NA)

# 1. 从数据框中提取
Methods2 <- c('Method',NA,NA,Methods1_new)
Pop2 <- c('Outcome population',NA,NA,new_dt$target)
Empty2 <- c(NA,NA,NA,new_dt$empty_column)

Beta2_3 <- ifelse(
  is.na(new_dt$beta),
  NA,  # 保持为逻辑缺失值 NA
  formatC(new_dt$beta, format = "f", digits = 3)
)
Beta2   <- c('Beta',NA,NA,Beta2_3)
Low2    <- c(NA,NA,NA,new_dt$loci)
Up2     <- c(NA,NA,NA,new_dt$upci)
CI2     <- c('95% CI',NA,NA,new_dt$CI_text)
Pval2   <- c('Pvalue',NA,NA,new_dt$p_formatted)
is_sum  <- c(TRUE,FALSE,FALSE,rep(FALSE,nrow(new_dt)))  # 如果没有汇总行，全部设为 FALSE

Beta2_2 <- as.numeric(Beta2)
Low2_2 <- as.numeric(Low2)
Up2_2 <- as.numeric(Up2)

Range1 <- c(Low2,Up2)
Upper1 <- quantile(c(Low2,Up2),0.97,na.rm = T)
Lower1 <- quantile(c(Low2,Up2),0.03,na.rm = T)
Range2 <- Range1[Range1 <= Upper1]
Range3 <- Range2[Range2 >= Lower1]
Upper2 <- max(pretty(Range3, n = 5))
Lower2 <- min(pretty(Range3, n = 5))
#22222::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::;:

dt_forest <- dt_forest_Asthma

dt_forest <- dt_forest %>% 
  mutate(empty_column = '                      ') %>%  # 创建空列
  #mutate(empty_column = NA) %>%  # 创建空列
  select(Exposure,Outcome,target,Method,empty_column,beta,betaCI2,sd,p,loci,upci) #%>%
#rename(` ` = empty_column)  # 将列名改为 "空格"


target_vec <- dt_forest$target
target_vec[target_vec == 1] <- 'Southern UK Population'
target_vec[target_vec == 0] <- 'Northern UK Population'
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
#=========================
Methods1 <- as.character(new_dt$Method)
# 创建一个逻辑向量，标记每个元素是否是该元素第一次出现
is_first <- !duplicated(Methods1)
# 创建新的向量，保留第一次出现的元素，其他设为 NA
Methods1_new <- ifelse(is_first, Methods1, NA)

# 1. 从数据框中提取
Methods2 <- c('Method',NA,NA,Methods1_new)
Pop2 <- c('Population',NA,NA,new_dt$target)
Empty2 <- c(NA,NA,NA,new_dt$empty_column)

Beta3_3 <- ifelse(
  is.na(new_dt$beta),
  NA,  # 保持为逻辑缺失值 NA
  formatC(new_dt$beta, format = "f", digits = 3)
)
Beta3   <- c('Beta',NA,NA,Beta3_3)
Low3    <- c(NA,NA,NA,new_dt$loci)
Up3     <- c(NA,NA,NA,new_dt$upci)
CI3     <- c('95% CI',NA,NA,new_dt$CI_text)
Pval3   <- c('Pvalue',NA,NA,new_dt$p_formatted)
is_sum  <- c(TRUE,FALSE,FALSE,rep(FALSE,nrow(new_dt)))  # 如果没有汇总行，全部设为 FALSE

Beta3_2 <- as.numeric(Beta3)
Low3_2 <- as.numeric(Low3)
Up3_2 <- as.numeric(Up3)

Range1 <- c(Low3,Up3)
Upper1 <- quantile(c(Low3,Up3),0.97,na.rm = T)
Lower1 <- quantile(c(Low3,Up3),0.03,na.rm = T)
Range2 <- Range1[Range1 <= Upper1]
Range3 <- Range2[Range2 >= Lower1]
Upper2 <- max(pretty(Range3, n = 5))
Lower2 <- min(pretty(Range3, n = 5))


#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::;:



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



pdf(file = paste0(path_output_plots, '/forest_BMI(2)', date_mark, '_', ver, '.pdf'),
    height = 12, width = 18.188)
grid.newpage()
pushViewport(viewport(layout = grid.layout(1, 2,widths=c(1.75,1))))
pushViewport(viewport(layout.pos.col = 1))
#11111
forestplot(
  #labeltext   = cbind(Methods2, Pop2, Empty2, Beta2, CI2, Pval2),
  labeltext   = cbind(Methods2, Pop2, Beta2, CI2, Pval2),
  mean        = Beta2_2,
  lower       = Low2_2,
  upper       = Up2_2,
  align       =c('l','l','l','c','c','c'),
  #xlab        = "Estimated regression coefficient for depression per 1‑unit increase in BMI",
  xlab        = 'COPD',
  #xlim        =c(-0.1,0.1),
  #xticks      = pretty(c(Low2, Up2), n = 5),
  xticks      = seq(-0.2,0.1,0.1),#pretty(Range3, n = 5),
  clip        = c(-0.2,0.1),#c(Lower2,Upper2),
  is.summary  = is_sum,#is.na(Beta2),#is_sum,
  #fn.ci_nrom  = 'fpDrawDiamondCI',
  hrzl_lines  = list(
    "3" = gpar(lwd=2, col="gray"),
    "5" = gpar(lwd=60, lineend="butt", columns=3:6, col="#99999920"),
    "8" = gpar(lwd=60, lineend="butt", columns=3:6, col="#99999920"),
    "11" = gpar(lwd=60, lineend="butt", columns=3:6, col="#99999920"),
    "14" = gpar(lwd=60, lineend="butt", columns=3:6, col="#99999920"),
    "17" = gpar(lwd=60, lineend="butt", columns=3:6, col="#99999920"),
    "20" = gpar(lwd=60, lineend="butt", columns=3:6, col="#99999920"),
    "23" = gpar(lwd=60, lineend="butt", columns=3:6, col="#99999920"),
    "26" = gpar(lwd=60, lineend="butt", columns=3:6, col="#99999920"),
    "29" = gpar(lwd=60, lineend="butt", columns=3:6, col="#99999920")
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
  new_page    = F
)

popViewport()
pushViewport(viewport(layout.pos.col = 2))
#22222
forestplot(
  #labeltext   = cbind(Methods2, Pop2, Empty2, Beta2, CI2, Pval2),
  labeltext   = cbind(Beta3, CI3, Pval3),
  mean        = Beta3_2,
  lower       = Low3_2,
  upper       = Up3_2,
  align       =c('l','c','c','c'),
  #xlab        = "Estimated regression coefficient for depression per 1‑unit increase in BMI",
  xlab        = 'Asthma',
  #xlim        =c(-0.1,0.1),
  #xticks      = pretty(c(Low2, Up2), n = 5),
  xticks      = seq(-0.3,0.2,0.2),#pretty(Range3, n = 5),
  clip        = c(-0.3,0.2),#c(Lower2,Upper2),
  is.summary  = is_sum,#is.na(Beta2),#is_sum,
  #fn.ci_nrom  = 'fpDrawDiamondCI',
  hrzl_lines  = list(
    "3" = gpar(lwd=2, col="gray"),
    "5" = gpar(lwd=60, lineend="butt", columns=1:4, col="#99999920"),
    "8" = gpar(lwd=60, lineend="butt", columns=1:4, col="#99999920"),
    "11" = gpar(lwd=60, lineend="butt", columns=1:4, col="#99999920"),
    "14" = gpar(lwd=60, lineend="butt", columns=1:4, col="#99999920"),
    "17" = gpar(lwd=60, lineend="butt", columns=1:4, col="#99999920"),
    "20" = gpar(lwd=60, lineend="butt", columns=1:4, col="#99999920"),
    "23" = gpar(lwd=60, lineend="butt", columns=1:4, col="#99999920"),
    "26" = gpar(lwd=60, lineend="butt", columns=1:4, col="#99999920"),
    "29" = gpar(lwd=60, lineend="butt", columns=1:4, col="#99999920")
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
  graph.pos   = 1,
  #graphwidth  = unit(12, 'cm'),
  new_page    = F
)


popViewport()
dev.off()



png(filename = paste0(path_output_plots, '/forest_BMI(2)', date_mark, '_', ver, '.png'),
    height = 12, width = 18, units = "in", res = 500)
# pdf(file = paste0(path_output_plots, '/forest_BMI(2)', date_mark, '_', ver, '.pdf'),
#     height = 12, width = 18.188)
grid.newpage()
pushViewport(viewport(layout = grid.layout(1, 2,widths=c(1.75,1))))
pushViewport(viewport(layout.pos.col = 1))
#11111
forestplot(
  #labeltext   = cbind(Methods2, Pop2, Empty2, Beta2, CI2, Pval2),
  labeltext   = cbind(Methods2, Pop2, Beta2, CI2, Pval2),
  mean        = Beta2_2,
  lower       = Low2_2,
  upper       = Up2_2,
  align       =c('l','l','l','c','c','c'),
  #xlab        = "Estimated regression coefficient for depression per 1‑unit increase in BMI",
  xlab        = 'COPD',
  #xlim        =c(-0.1,0.1),
  #xticks      = pretty(c(Low2, Up2), n = 5),
  xticks      = seq(-0.2,0.1,0.1),#pretty(Range3, n = 5),
  clip        = c(-0.2,0.1),#c(Lower2,Upper2),
  is.summary  = is_sum,#is.na(Beta2),#is_sum,
  #fn.ci_nrom  = 'fpDrawDiamondCI',
  hrzl_lines  = list(
    "3" = gpar(lwd=2, col="gray"),
    "5" = gpar(lwd=60, lineend="butt", columns=3:6, col="#99999920"),
    "8" = gpar(lwd=60, lineend="butt", columns=3:6, col="#99999920"),
    "11" = gpar(lwd=60, lineend="butt", columns=3:6, col="#99999920"),
    "14" = gpar(lwd=60, lineend="butt", columns=3:6, col="#99999920"),
    "17" = gpar(lwd=60, lineend="butt", columns=3:6, col="#99999920"),
    "20" = gpar(lwd=60, lineend="butt", columns=3:6, col="#99999920"),
    "23" = gpar(lwd=60, lineend="butt", columns=3:6, col="#99999920"),
    "26" = gpar(lwd=60, lineend="butt", columns=3:6, col="#99999920"),
    "29" = gpar(lwd=60, lineend="butt", columns=3:6, col="#99999920")
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
  new_page    = F
)

popViewport()
pushViewport(viewport(layout.pos.col = 2))
#22222
forestplot(
  #labeltext   = cbind(Methods2, Pop2, Empty2, Beta2, CI2, Pval2),
  labeltext   = cbind(Beta3, CI3, Pval3),
  mean        = Beta3_2,
  lower       = Low3_2,
  upper       = Up3_2,
  align       =c('l','c','c','c'),
  #xlab        = "Estimated regression coefficient for depression per 1‑unit increase in BMI",
  xlab        = 'Asthma',
  #xlim        =c(-0.1,0.1),
  #xticks      = pretty(c(Low2, Up2), n = 5),
  xticks      = seq(-0.3,0.2,0.2),#pretty(Range3, n = 5),
  clip        = c(-0.3,0.2),#c(Lower2,Upper2),
  is.summary  = is_sum,#is.na(Beta2),#is_sum,
  #fn.ci_nrom  = 'fpDrawDiamondCI',
  hrzl_lines  = list(
    "3" = gpar(lwd=2, col="gray"),
    "5" = gpar(lwd=60, lineend="butt", columns=1:4, col="#99999920"),
    "8" = gpar(lwd=60, lineend="butt", columns=1:4, col="#99999920"),
    "11" = gpar(lwd=60, lineend="butt", columns=1:4, col="#99999920"),
    "14" = gpar(lwd=60, lineend="butt", columns=1:4, col="#99999920"),
    "17" = gpar(lwd=60, lineend="butt", columns=1:4, col="#99999920"),
    "20" = gpar(lwd=60, lineend="butt", columns=1:4, col="#99999920"),
    "23" = gpar(lwd=60, lineend="butt", columns=1:4, col="#99999920"),
    "26" = gpar(lwd=60, lineend="butt", columns=1:4, col="#99999920"),
    "29" = gpar(lwd=60, lineend="butt", columns=1:4, col="#99999920")
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
  graph.pos   = 1,
  #graphwidth  = unit(12, 'cm'),
  new_page    = F
)


popViewport()
dev.off()
