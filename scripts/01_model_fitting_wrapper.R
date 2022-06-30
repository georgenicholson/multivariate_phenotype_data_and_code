
##################################################################
# This section imports arguments from the command line if present
arguments <- commandArgs()
if("--args" %in% arguments){
  argnam.in <- data.frame(nam = c("run_type", "analysis_table_row", "subsamseed"), 
                          coersion.fn = c("as.character", "as.numeric", "as.numeric"), 
                          stringsAsFactors = F)
  for(i in 1:nrow(argnam.in))
    assign(argnam.in$nam[i], eval(call(argnam.in$coersion.fn[i], arguments[grep("--args", arguments) + i])))
}

if (!all(c("run_type", "analysis_table_row", "subsamseed") %in% ls())) {
  stop("All of 'run_type', 'analysis_table_row', 'subsamseed' must be specified as inputs to 01_model_fitting_wrapper.R. 
       This can be done via running at the command line with 'run_type', 'analysis_table_row', 'subsamseed' as arguments in that order, 
       or by sourcing 01_model_fitting_wrapper.R in an R session with 'run_type', 'analysis_table_row', 'subsamseed' as objects. See README.md for details")
}

if (!"run_type" %in% ls()) {
  warning("run_type not specified, setting run_type <- 'demo' by default")
  ########################################################## 
  # Choose analysis type (see README)
  run_type <- c("demo", "main", "benchmark", "test_benchmark")[1]
}

library(foreach)
##########################################
# Source function files
fns_to_source <- list.files("scripts/functions", full.names = TRUE)
for (file_curr in fns_to_source) {
  source(file_curr)
}

##########################################
# control to contain parameters
control <- get_control_parameters_mv()

##########################################
# Create directory structure
dirs_to_create <- c("output_dir", "methods_comp_dir", "global_res_dir", "data_dir", "train_test_samples_dir",
                              "dropbox_figure_dir", "dropbox_text_numbers_dir", "dropbox_table_dir")
for (dirc in c(control[dirs_to_create])) {
  dir.create(dirc, recursive = TRUE, showWarnings = FALSE)
}

##########################################
# Load data
Data_all <- readRDS(file = control$Data_all_file)

##########################################
# Get table of analyses
analysis_table <- create_table_of_analyses(control = control, check_status = F, run_type = run_type)

##########################################
# Loop through analysis_table
# NOTE: You will need to parallelise this for the "main" or "benchmark" analyses 
# as each run is computationally intensive for the full sample size
# for (analysis_table_row in 1:nrow(analysis_table)) {
#   for (subsamseed in 1:analysis_table$n_subsamples[analysis_table_row]) {
    Data <- analysis_table$Data[analysis_table_row]
    Meth <- analysis_table$Meth[analysis_table_row]
    N <- analysis_table$N[analysis_table_row]
    P <- analysis_table$P[analysis_table_row]
    nSig <- analysis_table$nSig[analysis_table_row]
    Data <- analysis_table$Data[analysis_table_row]
    MVphen_K <- analysis_table$MVphen_K[analysis_table_row]
    n_subsamples <- analysis_table$n_subsamples[analysis_table_row]
    run_masking <- analysis_table$loocv[analysis_table_row]
    XDmeth <- Meth
    
    file_list <- get_file_list(control = control, 
                               file_core_name = analysis_table[analysis_table_row, "file_core_name"], 
                               subsamseed = subsamseed)
    
    Y_raw <- Data_all[[Data]]$Y_raw
    Y_zeroed <- Data_all[[Data]]$Y_zeroed
    S_raw <- Data_all[[Data]]$S_raw
    S_zeroed <- Data_all[[Data]]$S_zeroed
    Y.eml <- S.eml <- list()
    dat.typev <- switch(Data, 
                        impc = c("raw", "zero"), 
                        eqtl = "raw")
    Y.eml$raw <- Data_all[[Data]]$Y_raw
    Y.eml$zero <- Data_all[[Data]]$Y_zeroed
    S.eml$raw <- Data_all[[Data]]$S_raw
    S.eml$zero <- Data_all[[Data]]$S_zeroed
    
    ###################################
    # Create and Load train-test samples information 
    force_overwrite_train_test_splits <- FALSE
    train_test_list_file_curr <- file.path("output", "train_test_splits", paste0(Data, "_N_", N, "_P_", P, ".RDS"))
    if (!file.exists(train_test_list_file_curr) | force_overwrite_train_test_splits) {
      train_test_list <- get_train_test_split(control = control, Data_all = Data_all, N = N, P = P, Data = Data, n_subsamples = n_subsamples)
      saveRDS(object = train_test_list, file = train_test_list_file_curr)
    } else {
      train_test_list <- readRDS(file = train_test_list_file_curr)
    }
    phens_to_use <- train_test_list$phens_to_use
    
    ###################################
    # Assign samples
    for (var_assign_subsample in c("sams_for_lik_cross_val", 
                                   "sams_for_model_testing", 
                                   "sams_for_model_training", 
                                   "sams_for_cor_est",
                                   "sams_for_strong_cov_est")) {
      subsam_mat_curr <- train_test_list[[var_assign_subsample]]
      assign(x = var_assign_subsample, value = rownames(subsam_mat_curr)[subsam_mat_curr[, subsamseed]])
    }
    
    ##############################################
    # Estimate R
    if(Data == "impc"){
      R.init <- cor((Y_zeroed / S_zeroed)[sams_for_cor_est, phens_to_use], use = "p", meth = "p")
      R.init[is.na(R.init)] <- 0
      diag(R.init) <- 1
      if(qr(R.init)$rank < P)
        R.init <- (1 - control$rank_deficient_R_eps) * R.init + control$rank_deficient_R_eps * diag(rep(1, P))
    }
    if(Data == "eqtl"){
      data.in <- list(Bhat = Y_zeroed[sams_for_cor_est, phens_to_use], Shat = S_zeroed[sams_for_cor_est, phens_to_use])
      R.init <- mashr::estimate_null_correlation_simple(data = data.in, z_thresh = 2)
    }
    dimnames(R.init) <- list(phens_to_use, phens_to_use)
    
    ##############################################
    # Initilize Sig list
    Sigl.em.init <- initialize_Sig_list(control = control, 
                                        Y = Y_zeroed[sams_for_model_training, phens_to_use], 
                                        nSig = nSig, 
                                        random_init = grepl("rand", Meth))
    
    ##############################################
    # Run MVphen
    if (grepl("MVphen", Meth)) {
      source(file = "scripts/01a_fit_MVphen.R")
    }  
    
    ##################################################
    # Fit models using XD and/or MASH
    if(Meth %in% c("mash", "XD")){
      source(file = "scripts/01b_fit_XD_and_mash.R")
    }
    print("Done!")
#   }
# }
    
