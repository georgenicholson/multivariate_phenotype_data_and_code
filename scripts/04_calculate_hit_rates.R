Data <- "impc"

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

resl.comp <- readRDS(file = control$file.resl.comp)
compl <- readRDS(file = control$file.compl)
objl <- readRDS(file = control$file.objl)
resl.comp.fac <- readRDS(file = control$file_raw_factor_results)

resl.err.rates.comb <- c(resl.comp[grepl(Data, names(resl.comp)) | names(resl.comp) == "uv"],
                          list(varimax = resl.comp.fac[[control$mv_meth_nam_use]]))
split.use <- 1
resl.err.rates.one.split <- c(list(uv = resl.comp$uv), lapply(compl[grepl(Data, names(compl))], 
                                   function(x) list(mn = x$mnarr[, , split.use], sd = x$sdarr[, , split.use], lfsr = x$lfsrarr[, , split.use])))

for(err.data.type in c("comb", "single")[1]){
  resl.err.rates <- switch(err.data.type, 
                           comb = resl.err.rates.comb, 
                           single = resl.err.rates.one.split)
  
  out.perm <- err.rate.control(control = control, 
                               resl = resl.err.rates,
                               err.rate.meth = "perm", 
                               test.stat = "z", 
                               linemap = Data_all$impc$linemap, 
                               reflinemap = Data_all$impc$reflinemap, 
                               phmap = Data_all$impc$phmap, 
                               cenmap = Data_all$impc$cenmap,
                               Yhat = Data_all$impc$Y_raw,
                               control.level = c("line.fdr", "line.fwer", "phcen.fdr")[1],
                               err.thresh = control$fdr.th, 
                               p.complete.null.true = 1, 
                               p.test.null.true = 1)
  print(out.perm$restab)
  out.perm.lfsr <- err.rate.control(control = control, 
                               resl = resl.err.rates[!names(resl.err.rates) %in% c("uv", "uv.ss")],
                               err.rate.meth = "perm", 
                               test.stat = "lfsr", 
                               linemap = Data_all$impc$linemap, 
                               reflinemap = Data_all$impc$reflinemap, 
                               phmap = Data_all$impc$phmap, 
                               cenmap = Data_all$impc$cenmap,
                               Yhat = Data_all$impc$Y_raw,
                               control.level = c("line.fdr", "line.fwer", "phcen.fdr")[1],
                               err.thresh = control$fdr.th, 
                               p.complete.null.true = 1, 
                               p.test.null.true = 1)
  
  print(out.perm.lfsr$restab)
  out.lfsr <- err.rate.control(control = control, 
                                    resl = resl.err.rates[!names(resl.err.rates) %in% c("uv", "uv.ss")],
                                    err.rate.meth = "lfsr", 
                                    test.stat = "lfsr", 
                                    linemap = Data_all$impc$linemap, 
                                    reflinemap = Data_all$impc$reflinemap, 
                                    phmap = Data_all$impc$phmap, 
                                    cenmap = Data_all$impc$cenmap,
                                    Yhat = Data_all$impc$Y_raw,
                                    control.level = c("line.fdr", "line.fwer", "phcen.fdr")[1],
                                    err.thresh = control$fdr.th, 
                                    p.complete.null.true = 1, 
                                    p.test.null.true = 1)
  print(out.lfsr$restab)
  restabl <- list(perm = out.perm$restab, perm.lfsr = out.perm.lfsr$restab, lfsr = out.lfsr$restab)
  resll <- list(perm = out.perm$resl, perm.lfsr = out.perm.lfsr$resl, lfsr = out.lfsr$resl)
  resimpl <- list(perm = out.perm$resimp, perm.lfsr = out.perm.lfsr$resimp, lfsr = out.lfsr$resimp)
  
  saveRDS(resll, file = file.path(control$global_res_dir, paste0("resll_", err.data.type, ".RDS")))
  saveRDS(resimpl, file = file.path(control$global_res_dir, paste0("resimpl_", err.data.type, ".RDS")))
  saveRDS(restabl, file = file.path(control$global_res_dir, paste0("restabl_", err.data.type, ".RDS")))
}  

resimpl <- readRDS(file = file.path(control$global_res_dir, "resimpl_comb.RDS"))
resimp <- resimpl$perm
names(resimp)
resimp$eb.signsig <- resimp[, paste0(control$mv_meth_nam_use, ".perm.signsig")]
resimp$eb.sig <- abs(resimp$eb.signsig)
resimp$eb.t <- resimp[, paste0(control$mv_meth_nam_use, ".t")]
resimp$eb.th.final <- resimp[, paste0(control$mv_meth_nam_use, ".th.final")]
resimp <- resimp[resimp$line.type == control$nam.truemut, ]
resimp$uv.sig <- abs(resimp$uv.perm.signsig)
resimp$uv.signsig <- resimp$uv.perm.signsig

resimp_nonfac_rows <- which(!grepl("fac", resimp$ph))
resimp$eb.mn.loo.proc<-resimp$eb.sd.loo.proc <- NA
resimp$eb.mn.loo.proc[resimp_nonfac_rows] <- compl[[control$mv_meth_nam_use]]$loocv.mn.comb[cbind(resimp$geno, resimp$ph)[resimp_nonfac_rows, ]]
resimp$eb.sd.loo.proc[resimp_nonfac_rows] <- compl[[control$mv_meth_nam_use]]$loocv.sd.comb[cbind(resimp$geno, resimp$ph)[resimp_nonfac_rows, ]]
resimp$eb.t.loo.proc <- resimp$eb.mn.loo.proc / resimp$eb.sd.loo.proc


saveRDS(resimp, file = control$file.resimp)  


######################################################################
# Save smaller version of resll for reproducing figures
resll <- readRDS(file = file.path(control$global_res_dir, paste0("resll_comb.RDS")))
resll_sub <- list()
resll_sub$perm <- resll$perm[c("uv", control$mv_meth_nam_use)]
saveRDS(resll_sub, file = control$file_ref_lines_raw_outputs)
























