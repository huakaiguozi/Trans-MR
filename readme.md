# Trans-MR Project



### 1. 代码组织结构

 `project/
│
├── data/
│   ├── raw/            # 原始数据
│   └── processed/      # 处理后的数据
│
├── src/
│   ├── 01_setup.R      # 环境设置和依赖加载
│   ├── 02_data_cleaning.R  # 数据整理
│   ├── 03_data_analysis.R  # 数据分析
│   ├── 04_visualization.R  # 绘图
│   └── 05_report.Rmd      # 报告生成
│
├── output/
│   ├── tables/          # 输出的表格
│   ├── plots/           # 生成的图表
│   └── reports/         # 最终报告
│
├── renv/                # `renv` 环境文件夹
├── renv.lock            # `renv` 依赖锁定文件
├── .gitignore           # Git 忽略文件
├── README.md            # 项目说明
└── run_all.R            # 主执行脚本`

### 2. 