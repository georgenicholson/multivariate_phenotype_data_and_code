# multivariate_phenotype_data_and_code

This repository contains the data and R scripts to reproduce the results reported 
in 'Illuminating the mammalian genome with multivariate phenotype analysis'. 

## Installation

To run these scripts, you will need R version 3.6.3 or later, [available](https://www.r-project.org/) for Unix-like, Windows and Mac families of operating systems.

### Windows users please install Rtools
If you have R version 4.0.0 or above on Windows, you will also need to install 
[Rtools](https://cran.r-project.org/bin/windows/Rtools/). 

### Clone the repository
To get started, choose the local directory where you want the repository (let's call this <your_path>), and [clone](https://git-scm.com/book/en/v2/Git-Basics-Getting-a-Git-Repository)
this repository onto your local machine in folder "<your_path>/multivariate_phenotype_data_and_code"
```
cd <your_path>
git clone https://github.com/georgenicholson/multivariate_phenotype_data_and_code.git
```
### Install renv

Next, open an R console, and install the 
[renv](https://rstudio.github.io/renv/index.html) R package if you don't have it 
already
```
install.packages("renv")
```
###  Install the required packages
Run the following, remembering to place <your_path> in the code below
```
path_to_dir <- "<your_path>/multivariate_phenotype_data_and_code"
setwd(path_to_dir)
renv::activate()
install.packages("BiocManager")
BiocManager::install("AnnotationDbi", force = TRUE)
renv::restore()
```


## Data

Accompanying data are provided in "data/Data_all.RDS" which can be loaded in R via
```
Data_all <- readRDS("data/Data_all.RDS")
```

## Demo run

You can run a demonstration of the software on a smaller data set via
```
source("scripts/01_model_fitting_wrapper.R")
```
## Main analysis

To reproduce the main analysis, please re-define run_type at the top of "scripts/01_model_fitting_wrapper.R as follows
```
run_type <- "main"
```
The total CPU time of the main analysis is of the order of days, but can be parallelized over the 50 cross-validation folds (the subsamseed loop in "scripts/01_model_fitting_wrapper.R" and shown below)
```
for (scen in 1:nrow(analysis_table)) {
  # Parallelize the subsamseed (folds) loop for the "main" analysis
  for (subsamseed in 1:analysis_table$n_subsamples[scen]) {

  }
}
```
## Benchmarking analysis

To reproduce the benchmarking analysis, please re-define run_type at the top of "scripts/01_model_fitting_wrapper.R as follows
```
run_type <- "benchmark"
```
The total CPU time of the benchmarking analysis is of the order of months, but can be parallelized over the models and cross-validation folds:
```
# Parallelize the scen (models) and subsamseed (folds) loops for the "benchmark" analysis
for (scen in 1:nrow(analysis_table)) {
  for (subsamseed in 1:analysis_table$n_subsamples[scen]) {

  }
}
```

