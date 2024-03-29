create_table_of_analyses <- function(control, check_status = F, run_type = c("demo", "main", "benchmark", "test_benchmark")[1]) {
  
  if(run_type == "demo") {
    runtab <- data.frame(N = 100, P = 10, Data = "impc", 
                          Meth = c("MVphen", "XD", "mash"), 
                          nSig = c(2, 2, 1), 
                          MVphen_K = c(5, NA, NA), 
                          mem = 2000,
                          n_subsamples = 1,
                          stringsAsFactors = F)[1, ]
  }  
  
  if(run_type == "main") {
    runtab <- data.frame(N = 2000, P = 148, Data = "impc", 
                         Meth = "MVphen", 
                         nSig = 1, 
                         MVphen_K = 20, 
                         mem = 2000,
                         n_subsamples = control$n_subsamples_main,
                         stringsAsFactors = F)
  }  
  
  
  if (run_type %in% c("benchmark", "test_benchmark")) {
    runtab <- expand.grid(N = NA, P = NA, Data = c("impc", "eqtl"), 
                          Meth = c("MVphen", "XD", "mash", "MVphen_rand", "MVphen_N_500"), 
                          nSig = 1:2, 
                          MVphen_K = c(15, 20, 30, 40),
                          n_subsamples = control$n_subsamples_benchmark,
                          stringsAsFactors = F)
    runtab[which(!grepl("MVphen", runtab$Meth)), "MVphen_K"] <- NA
    #############################
    # Filter and tweak runtab
    #############################
    runtab[runtab$Meth == "mash", c("nSig")] <- 1
    runtab[runtab$Data == "impc", "P"] <- control$default_parameters$impc$P
    runtab[runtab$Data == "impc", "N"] <- control$default_parameters$impc$N
    runtab[runtab$Data == "eqtl", "P"] <- control$default_parameters$eqtl$P
    runtab[runtab$Data == "eqtl", "N"] <- control$default_parameters$eqtl$N
    runtab[runtab$Meth == "MVphen_N_500", "N"] <- 500
    runtab <- runtab[which(!(runtab$Meth == "MVphen_N_500" & (runtab$Data == "eqtl" | runtab$nSig > 1 | runtab$MVphen_K != 20))), ]
    runtab <- runtab[which(!(runtab$Meth == "MVphen_rand" & (runtab$Data == "eqtl" | runtab$MVphen_K != 20))), ]
    runtab <- runtab[which(!(runtab$Meth == "MVphen_rand" & runtab$nSig > 1)), ]
  }  
  
  ###############################################################
  #Specify memory requirements for each type of job
  #################################################################
  runtab$mem <- 2300
  runtab[runtab$Data %in% c("impc", "eqtl") & runtab$Meth == "MVphen", "mem"] <- 2000
  runtab[runtab$Data  == "impc" & runtab$Meth == "mash", "mem"] <- 6000
  runtab[runtab$Data == "eqtl" & runtab$Meth == "mash", "mem"] <- 4600
  runtab <- runtab[order(runtab$Data, runtab$Meth, runtab$nSig), ]
  rownames(runtab) <- 1:nrow(runtab)
  if (run_type == "test_benchmark") {
    runtab$n_subsamples <- 2  
    runtab$N <- 200  
    runtab$P <- 20
    runtab$MVphen_K <- 5
  }
  runtab <- unique(runtab)
  
  # Increased n_subsamples for main analysis
  runtab[runtab$N == 2000 & runtab$P == 148 & runtab$nSig == 1 & 
           runtab$Meth == "MVphen" & runtab$MVphen_K == control$nfac, "n_subsamples"] <- control$n_subsamples_main
  runtab$rand <- ifelse(grepl("rand", runtab$Meth), T, F)
  runtab$loocv <- ifelse(runtab$Meth == "MVphen" & runtab$nSig == 1 & runtab$N == 2000 & runtab$MVphen_K == 20, T, F)
  
  runtab[, c("file_core_name")] <- NA
  for (analysis_table_row in 1:nrow(runtab)) {
    for (j in 1:ncol(runtab)) {
      assign(colnames(runtab)[j], runtab[analysis_table_row, j], pos = sys.frame(which = 0))
    }
    variables_in_filename_use <- control$variables_in_filename
    XDmeth <- Meth
    file_core_name <- c()
    for (subsamseed in 1:n_subsamples) {
      file_core_name <- c(file_core_name, paste0(paste(paste(variables_in_filename_use, 
                                         sapply(variables_in_filename_use, function(x) get(x, envir = environment())), sep = "_"), 
                                         collapse = "_")))
    }
    loocv_results_filename <- paste0(control$methods_comp_dir, "/", file_core_name, "_loocv_res.RData")
    factor_results_filename <- paste0(control$methods_comp_dir, "/", file_core_name, "_facres.RData")
    basic_results_filename <- paste0(control$methods_comp_dir, "/", file_core_name, 
                             switch(Meth, eb = "_res.RData", mash = "_mash_resl.RData", XD = "_bovy_resl.RData"))
    emout_file_name <- paste0(control$methods_comp_dir, "/", file_core_name, 
                              switch(Meth, eb = "_emout.RData", mash = NA, XD = NA))
    runtab$file_core_name[analysis_table_row] <- gsub("seed\\_1", "seed\\_XXX", file_core_name[1])
  }
  rownames(runtab) <- 1:nrow(runtab)
  runtab <- runtab[order(match(runtab$Data, c("impc", "eqtl")),
                         match(runtab$Meth, c("XD", "mash", "MVphen", "MVphen_N_500", "MVphen_rand"))), ]
  runtab[runtab$Meth == "mash" & runtab$Data == "impc", "nSig"] <- 158
  runtab[runtab$Meth == "mash" & runtab$Data == "eqtl", "nSig"] <- 54
  runtab
  
  
  return(runtab)  
}


