install.packages("renv")
renv::activate()
renv::init(bioconductor = "3.13")
install.packages("BiocManager")
BiocManager::install(version = '3.13')
BiocManager::install(pkgs = c("AnnotationDbi", "Mus.musculus"), force = TRUE, version = "3.13")
renv::restore()
Data_all <- readRDS("data/Data_all.RDS")
list.files("output/global_results")
source("scripts/05_generate_results.R")
run_type <- "demo"
source("scripts/01_model_fitting_wrapper.R")

