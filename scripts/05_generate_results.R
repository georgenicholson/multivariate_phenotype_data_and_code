##########################################
# Source function files
fns_to_source <- list.files("scripts/functions", full.names = TRUE)
for (file_curr in fns_to_source) {
  source(file_curr)
}

##########################################
# control contains parameter settings
control <- get_control_parameters_mv()

###################################
# Load data
Data_all <- readRDS(control$Data_all_file)
linemap <- Data_all$impc$linemap
phmap <- Data_all$impc$phmap

resl.comp <- readRDS(file = control$file.resl.comp)
# compl <- readRDS(file = control$file.compl)
objl <- readRDS(file = control$file.objl)
resll <- readRDS(file = control$file.resll)
restabl <- readRDS(file = file.path(control$global_res_dir, paste0("restabl_comb.RDS")))

# This script generates Figures 2 and 3
source("scripts/05a_heatmap_of_hits.R")
# This script generates Figures S3, 5, and S4
source("scripts/05b_power_comparisons.R")
# This script generates Figures 6, S5, S6
source("scripts/05c_reference_lines.R")
# Please note: 05d_gene_ontology.R will crash if the Bioconductor packages "AnnotationDbi" and "Mus.musculus" 
# are not installed (see README for installation instructions)
try ({
  # This script generates Figures 7, S7
  source("scripts/05d_gene_ontology.R")
})
# This script generates Figures S8, 8
source("scripts/05e_factor_model.R")
# Please note: 05f_benchmarking.R will crash if you have generated your own set of results 
# using run_type == "main" instead of run_type == "benchmark"
try ({
  # This script generates Figure S9, and Tables 9, 10, 11 
  source("scripts/05f_benchmarking.R")
})
# This script generates Figure 4, and Tables 1a and 1b
source("scripts/05g_ebi_results_comparison.R")
# This script generates Figure S2
source("scripts/05h_create_figure_of_uv_qc.R")
# This script generates Figure 1
source("scripts/05i_create_figure_1_cartoon.R")
# This script generates Figure S10
source("scripts/06a_sensitivity_analysis.R")
# This script generates Table 8
source("scripts/06b_data_subsampling.R")
# This script generates Figure S11
source("scripts/06c_data_masking.R")
