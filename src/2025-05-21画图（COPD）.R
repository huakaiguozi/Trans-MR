# 2. 调用 forestplot
library(forestplot)


dt_forest <- dt_forest_COPD
dt_forest <- dt_forest %>%
  mutate(source = as.character(1 - as.numeric(target)), .after = Outcome)


dt_forest <- dt_forest %>% 
  mutate(empty_column = '                      ') %>%  # 创建空列
  #mutate(empty_column = NA) %>%  # 创建空列
  select(Exposure,Outcome,source,target,Method,empty_column,beta,betaCI2,sd,p,loci,upci) #%>%
#rename(` ` = empty_column)  # 将列名改为 "空格"


target_vec <- dt_forest$target
target_vec[target_vec == 1] <- 'Southern UK'
target_vec[target_vec == 0] <- 'Northern UK'
dt_forest$target <- target_vec

source_vec <- dt_forest$source
source_vec[source_vec == 1] <- 'Southern UK'
source_vec[source_vec == 0] <- 'Northern UK'
dt_forest$source <- source_vec
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
groups <- split(dt_forest, (seq_len(nrow(dt_forest)) - 1) %/% 1)

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
Pop2 <- c('Target',NA,NA,new_dt$target)
Pop2_2 <- c('Source',NA,NA,new_dt$source)
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

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::;:

#===========================
# 新增一个参数，传给 forestplot 的 `fn.ci_norm`
# fn <- local({
#   i <- 0
#   highlight_rows <- c(1:9)  # 你希望高亮的行号（可改）
#   
#   function(..., clr.line, clr.marker) {
#     i <<- i + 1
#     # 画背景色
#     if (i %in% highlight_rows) {
#       grid::grid.rect(gp = grid::gpar(fill = "#99999920", col = NA))
#     }
#     # 再画 CI（你原本想画的东西）
#     fpDrawNormalCI(..., clr.line = "#000000", clr.marker = "#F9675C")
#   }
# })

# 全局计数器
i_global <- 0
highlight_rows <- c(1:9)
fn <- local({
  i_local <- 0  # 用于交替颜色的局部计数器
  function(..., clr.line, clr.marker) {
    i_global <<- i_global + 1
    i_local <<- i_local + 1
    
    # 超过18行时，重置
    if (i_global > 9) {
      i_global <<- 1
      i_local <<- 1
    }
    #高亮背景
    if (i_global %in% highlight_rows) {
      grid::grid.rect(gp = grid::gpar(fill = "#99999920", col = NA))
    }
    # 交替颜色
    clr <- if(i_local == 9){'#FFA630'} else "#F9675C"
    
    fpDrawNormalCI(..., clr.line = "#000000", clr.marker = clr)
  }
  
})



pdf(file = paste0(path_output_plots, '/forest_BMI_COPD', date_mark, '_', ver, '.pdf'),
    height = 10, width = 15)

forestplot(
  #labeltext   = cbind(Methods2, Pop2, Empty2, Beta2, CI2, Pval2),
  labeltext   = cbind(Methods2, Pop2_2,Pop2, Beta2, CI2, Pval2),
  mean        = Beta2_2,
  lower       = Low2_2,
  upper       = Up2_2,
  align       =c('l','l','l','l','c','c','c'),
  #xlab        = "Estimated regression coefficient for depression per 1‑unit increase in BMI",
  xlab        = 'COPD',
  #xlim        =c(-0.1,0.1),
  #xticks      = pretty(c(Low2, Up2), n = 5),
  xticks      = seq(-0.2,0.1,0.1),#pretty(Range3, n = 5),
  clip        = c(-0.2,0.1),#c(Lower2,Upper2),
  is.summary  = is_sum,#is.na(Beta2),#is_sum,
  #fn.ci_norm = fn.ci_fillbg,  # 使用带背景的绘图函数
  #fn.ci_nrom  = 'fpDrawDiamondCI',
  hrzl_lines  = list(
    "3" = gpar(lwd=2, col="gray")
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
  graph.pos   = 4,
  #graphwidth  = unit(12, 'cm'),
  new_page    = F
)

dev.off()


#===========================
# 全局计数器
i_global <- 0
highlight_rows <- c(1:9)
fn <- local({
  i_local <- 0  # 用于交替颜色的局部计数器
  function(..., clr.line, clr.marker) {
    i_global <<- i_global + 1
    i_local <<- i_local + 1
    
    # 超过18行时，重置
    if (i_global > 9) {
      i_global <<- 1
      i_local <<- 1
    }
    #高亮背景
    if (i_global %in% highlight_rows) {
      grid::grid.rect(gp = grid::gpar(fill = "#99999920", col = NA))
    }
    # 交替颜色
    clr <- if(i_local == 9){'#FFA630'} else "#F9675C"
    
    fpDrawNormalCI(..., clr.line = "#000000", clr.marker = clr)
  }
  
})


png(filename = paste0(path_output_plots, '/forest_BMI_COPD', date_mark, '_', ver, '.png'),
    height = 10, width = 15, units = "in", res = 500)
forestplot(
  #labeltext   = cbind(Methods2, Pop2, Empty2, Beta2, CI2, Pval2),
  labeltext   = cbind(Methods2, Pop2_2,Pop2, Beta2, CI2, Pval2),
  mean        = Beta2_2,
  lower       = Low2_2,
  upper       = Up2_2,
  align       =c('l','l','l','l','c','c','c'),
  #xlab        = "Estimated regression coefficient for depression per 1‑unit increase in BMI",
  xlab        = 'COPD',
  #xlim        =c(-0.1,0.1),
  #xticks      = pretty(c(Low2, Up2), n = 5),
  xticks      = seq(-0.2,0.1,0.1),#pretty(Range3, n = 5),
  clip        = c(-0.2,0.1),#c(Lower2,Upper2),
  is.summary  = is_sum,#is.na(Beta2),#is_sum,
  #fn.ci_norm = fn.ci_fillbg,  # 使用带背景的绘图函数
  #fn.ci_nrom  = 'fpDrawDiamondCI',
  hrzl_lines  = list(
    "3" = gpar(lwd=2, col="gray")
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
  graph.pos   = 4,
  #graphwidth  = unit(12, 'cm'),
  new_page    = F
)

dev.off()
