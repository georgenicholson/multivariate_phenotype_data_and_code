################################################
# This script generates Figure S10
################################################

#####################################################################
# Calculate KL divergence between split models and combined model for Factor Sensitivity Analysis
kl1 <- kl2 <- c()
Sig.comb <- resl.comp[[control$mv_meth_nam_use]]$Sig.comb
P <- control$default_parameters$impc$P
for(seed in 1:control$n_subsamples_main){
  Sig_for_curr_data_split <- objl[[control$mv_meth_nam_use]][[seed]]$Sigl[[1]]
  if(!is.null(Sig_for_curr_data_split)){
    kl1[seed] <- .5 * (sum(diag(solve(Sig.comb) %*% Sig_for_curr_data_split)) - P + 
                         determinant(Sig.comb, logarithm = T)$modulus - determinant(Sig_for_curr_data_split, logarithm = T)$modulus)
    kl2[seed] <- .5 * (sum(diag(solve(Sig_for_curr_data_split) %*% Sig.comb)) - P + 
                         determinant(Sig_for_curr_data_split, logarithm = T)$modulus - determinant(Sig.comb, logarithm = T)$modulus)
  }
}

seed.largest.kl <- which.max(kl1 + kl2)
Sig_main <- Sig.comb
Sig_compare <- objl[[control$mv_meth_nam_use]][[seed.largest.kl]]$Sigl[[1]]
Sig_corr_compare <- t(Sig_compare / sqrt(diag(Sig_compare))) / sqrt(diag(Sig_compare))
eigc_compare <- eigen(Sig_corr_compare)
Sig_corr_main <- t(Sig_main / sqrt(diag(Sig_main))) / sqrt(diag(Sig_main))
eigc_main <- eigen(Sig_corr_main)
loadings_main <- varimax(eigc_main$vectors[, 1:control$nfac])$loadings
loadings_compare <- varimax(eigc_compare$vectors[, 1:control$nfac])$loadings
loadings_main <- sweep(loadings_main, 2, apply(loadings_main, 2, function(v) v[which.max(abs(v))]), "/")
loadings_compare <- sweep(loadings_compare, 2, apply(loadings_compare, 2, function(v) v[which.max(abs(v))]), "/")
rownames(loadings_compare) <- rownames(Sig_corr_compare)
rownames(loadings_main) <- rownames(Sig_corr_main)
loadings_compare <- loadings_compare[rownames(loadings_main), ]
################################################
# Reorder compare loadings to match main
prox_mat <- abs(t(loadings_main) %*% loadings_compare)
map_to <- rep(0, control$nfac)
for (j in 1:control$nfac) {
    map_to[j] <- which.max(ifelse(1:control$nfac %in% map_to, -1, abs(t(loadings_main) %*% loadings_compare[, j])))
}

loadings_compare[, map_to] <- loadings_compare
loadings_compare <- sweep(loadings_compare, 2, sign(colMeans(loadings_main * loadings_compare)), "*")


jpegc <- T
devwid <- 9
devhei <- 11
if(jpegc){
  fnamc <- "sensitivity_analysis_plot.jpg"
  jpeg(paste(control$figure_dir, "/", fnamc, sep = ""), devwid, 11, units = "in", res = 1000)
} else {
  graphics.off()
}
wc <- .6
hc <- .8
laymat <- matrix(c(1, 2, 4, 3), 2, 2)
layl <- list(p1 = c(2, 2), p2 = c(3, 2), p3 = c(3, 3))
npl <- length(layl)
laymat <- matrix(npl + 1, 4, 4)
for(i in 1:length(layl))
  laymat[t(layl[[i]])] <- i
laymat
vp <- c(.025, .6, .2, .2)
hp <- c(.3, .25, .2, .3)
layout(laymat, hei = vp, wid = hp)
oma.use <- rep(0, 4)
par(mar = rep(.5, 4), oma = oma.use, xpd = F)
cexax <- .9
cexax2 <- .9
cexax3 <- .95
cexax4 <- .8
cexlet <- 1.15
cexleg <- 1
cexax.labrt <- .6
cexax.lablt <- 1
cex.axis <- .9
cexax.mtext.big <- .8
nfac <- control$nfac
nph <- nrow(loadings_main)
phord.fac <- Data_all$impc$phord

par(mfrow = c(1, 2), oma = c(4, 18, 2, 2), mar = c(1, 1, 1, 1))
zpl <- t(loadings_main)
for (j in 1:2) {
  if (j == 1) {
    zpl <- t(loadings_main[phord.fac, ])
  }
  if (j == 2) {
    zpl <- t(loadings_compare[phord.fac, ])
  }
  image(x = 1:nfac, y = 1:nph, z = zpl, col = control$heat_col_palette, yaxt = "n", zlim = c(-1, 1),
        ylab = "", cex.axis = cexax, xaxt = "n", xlab = "")
  ypl1 <- grconvertY(c(0, 1), from = "nfc", "ndc")
  abline(v = 1:nfac - .5)
  procv <- Data_all$impc$phmap[match(phord.fac, Data_all$impc$phmap$ph), "procnam"]
  procats <- sapply(Data_all$impc$procord, function(procc) mean(which(procv == procc)))
  if (j == 1) {
    axis(side = 2, labels = Data_all$impc$procord, at = procats, las = 2, cex.axis = cexax.lablt)
  }
  title_curr <- c(expression("(a) ", paste(Sigma[pool], " loadings")), 
                  expression("(b) ", paste(Sigma[maxKL], " loadings")))[j]
  mtext(side = 3, text = title_curr, line = .5, cex = 1.2)
  abline(h = match(Data_all$impc$procord, procv) - .5, lwd = 2)
}

if(jpegc){
  dev.off()
  if (control$output_to_dropbox) {
    file.copy(from = paste(control$figure_dir, "/", fnamc, sep = ""),
              to = paste(control$dropbox_figure_dir, "/", fnamc, sep = ""), overwrite = TRUE)
  } else {
    file.rename(from = paste(control$figure_dir, "/", fnamc, sep = ""), 
                to = file.path(control$figure_dir, "Figure_S10.jpg"))
  }
}



