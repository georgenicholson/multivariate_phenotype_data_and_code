# multivariate_phenotype_data_and_code

This repository contains the data and R scripts to reproduce the results reported 
in 'Illuminating the mammalian genome with multivariate phenotype analysis'. 

## Installation

To run these scripts, you will need R version 3.6.3 or later, [available](https://www.r-project.org/) for Unix-like, Windows and Mac families of operating systems.

### Windows users install Rtools
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
Run the following, changing <your_path> to get to the path of your local version of this repository
```
path_to_dir <- "<your_path>/multivariate_phenotype_data_and_code"
setwd(path_to_dir)
renv::activate(path_to_dir)
renv::restore(path_to_dir)
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

## Full analysis


