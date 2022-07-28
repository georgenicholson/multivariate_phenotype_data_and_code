##########################################
# Source function files
fns_to_source <- list.files("scripts/functions", full.names = TRUE)
for (file_curr in fns_to_source) {
  source(file_curr)
}

##########################################
# control contains parameter settings
control <- get_control_parameters_mv(output_to_dropbox = F)

##########################################
# Create directory structure
dirs_to_create <- c("output_dir", "methods_comp_dir", "global_res_dir", "data_dir", "train_test_samples_dir",
                    "dropbox_figure_dir", "dropbox_text_numbers_dir", "dropbox_table_dir")
for (dirc in c(control[dirs_to_create])) {
  dir.create(dirc, recursive = TRUE, showWarnings = FALSE)
}

###################################
# Load data
Data_all <- readRDS(control$Data_all_file)
linemap <- Data_all$impc$linemap
phmap <- Data_all$impc$phmap

resl.comp <- readRDS(file = control$file.resl.comp)
objl <- readRDS(file = control$file.objl)
resll <- readRDS(control$file_ref_lines_raw_outputs)
restabl <- readRDS(file = file.path(control$global_res_dir, paste0("restabl_comb.RDS")))


# This script generates Figures 2 and 3
source("scripts/05a_heatmap_of_hits.R")
# This script generates Figures S3, 5, and S4
source("scripts/05b_power_comparisons.R")
# This script generates Figures 6, S5, S6
source("scripts/05c_reference_lines.R")
try ({
  # This script generates Figures 8, S7
  source("scripts/05d_gene_ontology.R")
})
# This script generates Figures S8, 10
source("scripts/05e_factor_model.R")
# Please note: 05f_benchmarking.R will crash if you have 
# generated a set of incomplete benchmarking results 
# that overwrite the files in output/global_results/ .
# Re-downloading and unzipping output.zip fixes this
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
