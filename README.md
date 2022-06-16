# multivariate_phenotype_data_and_code

This repository contains the data and R scripts to reproduce the results reported 
in 'Multivariate phenotype analysis enables genome-wide inference of mammalian gene function'. 

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
renv::init(bioconductor = "3.13")
install.packages("BiocManager")
BiocManager::install(version = '3.13')
BiocManager::install(pkgs = c("AnnotationDbi", "Mus.musculus"), force = TRUE, version = "3.13")
renv::restore()
```

## Data

Accompanying data are provided in "data/Data_all.RDS" which can be loaded in R via
```
Data_all <- readRDS("data/Data_all.RDS")
```

## Generate figures from paper

To generate the Figures and Tables in the paper, download and unzip the analysis outputs file "output.zip", an asset in our [GitHub release](https://github.com/georgenicholson/multivariate_phenotype_data_and_code/releases/tag/v1.1). Unzip it directly in the root of the repository folder, creating a folder structure "output/global_results/" containing the precomputed output files. Check the files are there:
```
list.files("output/global_results")
```
Once the files are in place please run
```
source("scripts/05_generate_results.R")
```
The Figures will be generated in the [/figures](figures) and [/tables](tables) folders. The script takes as input the raw data as well as the outputs of model fits, stored in [/output/global_results](output/global_results)

## Demo run of software

You can run a demonstration of the software on a smaller data set via
```
run_type <- "demo"
source("scripts/01_model_fitting_wrapper.R")
```
## Main analysis

To reproduce the main analysis (without benchmarking against a variety of alternative methods), the total CPU time is of the order of days, but can be parallelized over the 50 cross-validation folds (the subsamseed loop in "scripts/01_model_fitting_wrapper.R" and shown below)
```
for (scen in 1:nrow(analysis_table)) {
  # Parallelize the subsamseed (folds) loop for the "main" analysis
  for (subsamseed in 1:analysis_table$n_subsamples[scen]) {
    < Loop content omitted, full code in "scripts/01_model_fitting_wrapper.R" >
  }
}
```
After parallelizing the CV folds above, you can run 
```
run_type <- "main"
source("scripts/01_model_fitting_wrapper.R")
```

## Full benchmarking analysis

To reproduce the full benchmarking analysis, the total CPU time is of the order of months, but can be parallelized over the models and cross-validation folds at this point in "scripts/01_model_fitting_wrapper.R":
```
# Parallelize the scen (models) and subsamseed (folds) loops for the "benchmark" analysis
for (scen in 1:nrow(analysis_table)) {
  for (subsamseed in 1:analysis_table$n_subsamples[scen]) {
    < Loop content omitted, full code in "scripts/01_model_fitting_wrapper.R" >
  }
}
```
After parallelizing the CV folds above, you can run 
```
run_type <- "benchmark"
source("scripts/01_model_fitting_wrapper.R")
```

## Post model fitting analyses

To perform the post-model-fitting analyses and generate figures and tables, please run
```
source("scripts/02_collect_results.R")
source("scripts/03_estimate_global_factors.R")
source("scripts/04_calculate_hit_rates.R")
source("scripts/05_generate_results.R")

```


## 
