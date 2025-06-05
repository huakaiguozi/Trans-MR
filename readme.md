# TLMR Simulation and Application Code

This repository contains all simulation code for **TLMR (Transfer Learning based Mendelian Randomization)**, a novel approach for causal effect estimation under population heterogeneity using individual-level data. The method integrates transfer learning techniques into the Mendelian Randomization (MR) framework to improve the validity of causal inference when populations differ between source and target datasets.

---

## 📁 Repository Structure

```
Trans-MR/
├── data/                             # Input data folder for simulations or applications
├── output/                           # Output results folder from simulations and applications
├── renv/                             # R environment management files
├── src/                              # All the R scripts that may be reused
│   ├── requirements                  # All the packages used in the simulations or applications
│
│   ├──  ...
│   ├── # Main computation process codes for simulations and applications
│   ├── computation_categorical.R
│   ├── computation_simulation2.R
│   ├── computation_simulation3.R
│   ├──  ...
│
│   ├── # All the functions used for computing process in simulations and applications
│   ├── computing_functions_categorical.R
│   ├── computing_functions_continuous.R
│   ├── computing_functions_application.R
│   ├──  ...
│
│   ├──  # All the functions used for data generation process in simulations
│   ├── data_generating_functions_categorical.R
│   ├── data_generating_functions_continuous.R
│   ├── ...
│
│   ├── # All drawing codes for simulations and applications
│   ├── drawing_categorical.R
│   ├── drawing_continuous.R
│   ├── drawing_application1.R
│   ├──  ...
│
├── .Rprofile                         # R configuration file
├── .gitignore                        # Git ignore configuration
├── LICENSE                           # License file
├── Trans_MR_project.Rproj            # RStudio project file
├── application1.R                    # Main R script: real data application
├── simulation1_continuousoutcome_categoricalcovariates.R   # Main R script: simulation 1-1: Continuous outcome, categorical covariates
├── simulation1_continuousoutcome_continuouscovariates.R    # Main R script: simulation 1-2: Continuous outcome, continuous covariates
├── simulation1_continuousoutcome_misspecifiedORmodel.R     # Main R script: simulation 1-3: Continuous outcome, misspecified outcome model
├── simulation2_binaryoutcome.R       # Main R script: simulation 2: Binary outcome
├── renv.lock                         # R environment lock file
└── readme.md                         # Project description file
```

---

## ▶️ Getting Started

### Prerequisites

You need R and all the packages listed in /src/requirements.R installed. 

- Unless otherwise noted, these packages can be installed directly via install.packages('[package name]'). 
- Some packages require installation through GitHub or other sources; we have provided the corresponding installation commands as comments in requirements.R.

### Run Example Simulation

Simply run the main R script to obtain the result:

```r
source("simulation1_continuousoutcome_categoricalcovariates.R")
```
You can also open the main R scripts to modify the names of output folders and other settings before running them.

---

## 📬 Contact

**Author**: [huakaiguozi](https://github.com/huakaiguozi)  
**E-mail**: yw86023@outlook.com