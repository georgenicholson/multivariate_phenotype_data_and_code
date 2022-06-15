################################################
# This script generates Figures 6, S5, S6
################################################

resimp <- readRDS(file = control$file.resimp)

####################################################
# Reference lines numbers for text
resimp0.imp <- resimp[resimp$line.type == control$control$nam.negcon & resimp$imputed, ]
resimp1.imp <- resimp[resimp$line.type == control$control$nam.truemut & resimp$imputed, ]
glob.fdr.imp <- mean(abs(resimp0.imp$eb.t) > resimp0.imp$eb.th.final) / 
  mean(abs(resimp1.imp$eb.t) > resimp1.imp$eb.th.final)
resimp0.non <- resimp[resimp$line.type == control$nam.negcon & !resimp$imputed, ]
resimp1.non <- resimp[resimp$line.type == control$nam.truemut & !resimp$imputed, ]
glob.fdr.non <- mean(abs(resimp0.non$eb.t) > resimp0.non$eb.th.final) / 
  mean(abs(resimp1.non$eb.t) > resimp1.non$eb.th.final)

matl <- resll$perm[[control$mv_meth_nam_use]]$ref.lines$matl
matl$pre <- resll$perm$uv$ref.lines$matl$post
matl.sig <- lapply(matl, function(M) {M[abs(M) < 1] <- 0; sign(M)})
matl.sig.ord <- lapply(matl.sig, function(M)t(apply(M, 1, function(v) sort(v, decreasing = T))))
tabl <- resll$perm[[control$mv_meth_nam_use]]$ref.lines$tabl

fdr.ests <- resll$perm[[control$mv_meth_nam_use]]$ref.lines$fdr.ests
fdr.ests$pre <- resll$perm$uv$ref.lines$fdr.ests$post

#FDR ests and CIs
for (i in 1:length(fdr.ests)) {
  fdrc <- fdr.ests[[i]]
  loomv.uv.ci.numv <- formatC(unlist(fdrc) * 100, format = "f", digits = 1)
  loomv.uv.ci <- paste(loomv.uv.ci.numv[1], "\\% (95\\% CI: ", loomv.uv.ci.numv[2], "\\% - ", loomv.uv.ci.numv[3], "\\%)", sep = "")
  save.text <- paste("reflines.fdr.ci.", names(fdr.ests)[i], sep = "")
  for (numc in save.text)
    write.table(loomv.uv.ci, file = paste(control$dropbox_text_numbers_dir, "/", save.text, ".txt", sep = ""),
                col.names = F, row.names = F, quote = F)
}

#global concordance/discordance
ref.lines.concordant.n <- tabl$post["-1", "-1"] + tabl$post["1", "1"]
ref.lines.discordant.n <- tabl$post["-1", "1"] + tabl$post["1", "-1"]
save.num <- c("ref.lines.concordant.n", "ref.lines.discordant.n")
for (numc in save.num)
  write.table(prettyNum(eval(as.name(numc)), big.mark = ","), file = paste(control$dropbox_text_numbers_dir, "/", numc, ".txt", sep = ""),
              col.names = F, row.names = F, quote = F)


####################################################
#Reference lines concordance scatter (Figure 6)
####################################################
fnamc <- "ref_lines_z_scatter.jpg"
jpeg(paste(control$figure_dir, "/", fnamc, sep = ""), 9, 9, units = "in", res = 500)
par(mfrow = c(2, 2), oma = c(3, 2, 1, 1), mar = c(3, 3, 3, 3))
namv <- c("pre", "postcomp", "postimp", "postimpcomp")
for (namc in namv) {
  ntab <- table(matl.sig[[namc]][, 1], matl.sig[[namc]][, 2])
  ptab <- ntab / sum(ntab)
  limn <- 3
  legcex <- .65
  cexax <- .9
  plot(matl[[namc]][, 1], matl[[namc]][, 2], pch = 1, cex = .2, 
       xlim = c(-1, 1) * limn, ylim = c(-1, 1) * limn, ylab = "", xlab = "", xaxs = "i", yaxs = "i")
  abline(0, 1)
  labc1 <- switch(namc, pre = "UV replicate 1 (measured)",
                  postcomp = "MV replicate 1 (measured)",
                  postimp = "MV replicate 1 (missing)",
                  postimpcomp = "MV replicate 1 (missing)")
  labc2 <- switch(namc, pre = "UV replicate 2 (measured), ",
                  postcomp = "MV replicate 2 (measured)",
                  postimp = "MV replicate 2 (missing)",
                  postimpcomp = "MV replicate 2 (measured)")
  mtext(side = 1, line = 2.5, text = substitute(paste(labc1, ", ", italic(tilde(z))), list(labc1 = labc1)), cex = cexax)
  mtext(side = 2, line = 2.5, text = substitute(paste(labc2, ", ", italic(tilde(z))), list(labc2 = labc2)), cex = cexax)
  abline(h = c(-1, 1), col = 2, lty = 2)
  abline(v = c(-1, 1), col = 2, lty = 2)
  atv <- c((-limn - 1) / 2, 0, (limn + 1) / 2)
  xl <- list(c(1, 1, 3, 3), c(1, 1, 3, 3) * -1, c(1, 1, 3, 3), c(1, 1, 3, 3) * -1)
  yl <- list(c(1, 3, 3, 1), c(1, 3, 3, 1) * -1, c(1, 3, 3, 1) * -1, c(1, 3, 3, 1))
  for (quad in 1:4) {
    if(quad %in% 1:2)
      colc <- rgb(0, 0, 1, .35)
    if(quad %in% 3:4)
      colc <- rgb(1, 0, 0, .35)
    polygon(x = xl[[quad]], y = yl[[quad]], col = colc)
  }
  for (i in 1:3) {
    for (j in 1:3) {
      labc <- paste(prettyNum(ntab[i, j], big.mark = ","), " (", formatC(100 * ptab[i, j], format = "f", digits = 1), "%)", sep = "")
      legend(x = atv[i], y = atv[j], legend = labc, text.col = 1, bg = adjustcolor("white", .75), xjust = .5, cex = legcex)
    }
  }
  mtext(side = 3, text = paste("(", letters[match(namc, namv)], ")", sep = ""), at = -limn - .25, line = 1, cex = 1.2)
  numv <- formatC(c(fdr.ests[[namc]][[1]], fdr.ests[[namc]][[2]]) * 100, format = "f", digits = 1)
  methlab <- switch(namc, pre = "UV (measured)", postcomp = "MV (measured)", postimp = "MV (missing)", postimpcomp = "MV (meas. vs miss.)")
  fdrlab <- paste0(methlab, " Fsr = ", numv[1], "% (", numv[2], "% - ", numv[3], "%)")
  mtext(side = 3, text = fdrlab, line = .5, cex = .9)
}
dev.off()
if (control$output_to_dropbox) {
  file.copy(from = paste(control$figure_dir, "/", fnamc, sep = ""),
            to = paste(control$dropbox_figure_dir, "/", fnamc, sep = ""), overwrite = TRUE)
} else {
  file.rename(from = paste(control$figure_dir, "/", fnamc, sep = ""), 
              to = file.path(control$figure_dir, "Figure_6.jpg"))
}



#######################################################
# Heatmap of reference lines results (Figure S5)
#######################################################
graphics.off()
tl.pre <- resll$perm$uv$ref.lines$tl
thl.pre <- resll$perm$uv$ref.lines$thl
tl.post <- resll$perm[[control$mv_meth_nam_use]]$ref.lines$tl
thl.post <- resll$perm[[control$mv_meth_nam_use]]$ref.lines$thl

nlin <- length(tl.pre)
nlinrep <- sapply(tl.pre, ncol)
zth = 1
names(tl.pre)
big.mar <- 1.5
small.mar <- .5
fnamc <- "ref_lines_heatmap.jpg"
jpeg(paste(control$figure_dir, "/", fnamc, sep = ""), 12, 18, units = "in", res = 600)
layout(matrix(1:(2 * nlin), 1), width = rep(nlinrep, each = 2) + big.mar)
par(oma = c(6, 22, 0, 20))
martop <- 10
for (i in 1:length(tl.pre)) {
  for (ty in c("pre", "pos")) {
    cenv <- resimp[match(colnames(tl.pre[[i]]), resimp$geno), "cenlong"]
    nlinrepc <- nlinrep[i]
    if(ty == "pre") {
      plc = (tl.pre[[i]] / thl.pre[[i]])[Data_all$impc$phord, ]
      par(mar = c(1, big.mar / 2, martop, small.mar / 2))
    }
    if(ty == "pos") {
      plc = (tl.post[[i]] / thl.post[[i]])[Data_all$impc$phord, ]
      par(mar = c(1, small.mar / 2, martop, big.mar / 2))
    }
    plc[which(abs(plc) > zth)] = zth * sign(plc[which(abs(plc) > zth)])
    image(x = 1:ncol(plc), y = 1:nrow(plc), z = t(plc), col = control$heat_col_palette, zlim = c(-1, 1) * zth, xaxt = "n", yaxt = "n", xlab = "", ylab = "")
    for (j in 1:nrow(plc)) {
      for (k in 1:ncol(plc)) {
        if(any(na.omit(abs(plc[j, k])) == 1)) {
          lines(x = k + c(-.5, .5), y = j + c(-.5, .5))
          lines(x = k + c(.5, -.5), y = j + c(-.5, .5))
        }
      }
    }
    procv <- Data_all$impc$phmap[match(Data_all$impc$phord, Data_all$impc$phmap$ph), "procnam"]
    if(i == 1 & ty == "pre") {
      procats <- sapply(Data_all$impc$procord, function(procc) mean(which(procv == procc)))
      axis(side = 2, labels = Data_all$impc$procord, at = procats, las = 2, cex.axis = 1.2)
    }
    if(i == nlin & ty == "pos") {
      phnam <- Data_all$impc$phmap[match(Data_all$impc$phord, Data_all$impc$phmap$ph), "nam"]
      axis(side = 4, labels = phnam, at = 1:length(Data_all$impc$phord), las = 2, cex.axis = 1)
    }
    abline(h = match(Data_all$impc$procord, procv) - .5, lwd = 2)
    polygon(x = rep(c(0, nlinrepc), each = 2) + .5, y = c(0, nrow(plc))[c(1, 2, 2, 1)] + .5, lwd = 2)
    mtext(side = 3, text = ifelse(ty == "pre", "UV", "MV"), cex = .8)
    axis(side = 1, labels = cenv, at = 1:ncol(plc), cex.axis = .9, las = 2)
    splc <- strsplit(names(tl.pre)[i], spl = "_")[[1]]
    gnam1 <- quote(substitute(italic(a), list(a = splc)))
    gnam2 <- quote(substitute(paste("(", b, ")", sep = ""), list(b = ifelse(splc[2] == 1, "hom", "het"))))
    if(ty == "pre") {
      mtext(side = 3, at = nlinrepc + 1, line = 2.5, text = eval(gnam1), cex = 1)
      mtext(side = 3, at = nlinrepc + 1, line = 1, text = eval(gnam2), cex = 1)
    }
  }
}
par(fig = c(.3, .7, .96, .97), new = T)
cexax <- 1.1
par(mar = c(0, 0, 0, 0))
image(z = as.matrix(1:1000), x = seq(-1, 1, len = 1000), y = 1, col = control$heat_col_palette, xaxt = "n", yaxt = "n", 
      ylab = "")
axis(side = 3, las = 2, cex.axis = cexax, labels = c("< -1.0", -0.5, 0.0, 0.5, "> 1.0"), at = seq(-1, 1, by = .5), las = 0)
mtext(side = 3, text = expression(italic(tilde(z))), line = 2, cex = cexax, las = 0)
dev.off()
if (control$output_to_dropbox) {
  file.copy(from = paste(control$figure_dir, "/", fnamc, sep = ""),
            to = paste(control$dropbox_figure_dir, "/", fnamc, sep = ""), overwrite = TRUE)
} else {
  file.rename(from = paste(control$figure_dir, "/", fnamc, sep = ""), 
              to = file.path(control$figure_dir, "Figure_S5.jpg"))
}


hom.z <- (resll$perm[[control$mv_meth_nam_use]]$het.hom$tl$het.hom / resll$perm[[control$mv_meth_nam_use]]$het.hom$thl$het.hom)[, 1]
het.z <- (resll$perm[[control$mv_meth_nam_use]]$het.hom$tl$het.hom / resll$perm[[control$mv_meth_nam_use]]$het.hom$thl$het.hom)[, 2]
hom.z.sign <- sign(hom.z) * (abs(hom.z) > 1)
het.z.sign <- sign(het.z) * (abs(het.z) > 1)
tabhethom <- table(hom.z.sign, het.z.sign)

###############################################
# Het hom FDR ests and CIs (Figure S6)
###############################################
tabhethom <- table(hom.z.sign, het.z.sign)
fdrci <- fdr.est.tab(tabhethom)
fdrc <- fdrci
loomv.uv.ci.numv <- formatC(unlist(fdrc) * 100, format = "f", digits = 1)
loomv.uv.ci <- paste(loomv.uv.ci.numv[1], "\\% (95\\% CI: ", loomv.uv.ci.numv[2], "\\% - ", loomv.uv.ci.numv[3], "\\%)", sep = "")
save.text <- "hethom.fdr.ci"
for (numc in save.text)
  write.table(loomv.uv.ci, file = paste(control$dropbox_text_numbers_dir, "/", save.text, ".txt", sep = ""), 
              col.names = F, row.names = F, quote = F)



fnamc <- "hethom_scatter.jpg"
jpeg(paste(control$figure_dir, "/", fnamc, sep = ""), 9, 9, units = "in", res = 500)
par(mfrow = c(1, 1), oma = c(3, 2, 1, 1), mar = c(3, 3, 3, 3))
ntab <- t(tabhethom)
ptab <- ntab / sum(ntab)
limn <- 3
legcex <- .65
plot(het.z, hom.z, pch = 1, cex = .2, 
     xlim = c(-1, 1) * limn, ylim = c(-1, 1) * limn, ylab = "", xlab = "", xaxs = "i", yaxs = "i")
abline(0, 1)
labc1 <- "Heterozygote MV "
labc2 <- "Homozygote MV"
mtext(side = 1, line = 2.5, text = substitute(paste(labc1, " ", italic(tilde(z))), list(labc1 = labc1)))
mtext(side = 2, line = 2.5, text = substitute(paste(labc2, " ", italic(tilde(z))), list(labc2 = labc2)))
abline(h = c(-1, 1), col = 2, lty = 2)
abline(v = c(-1, 1), col = 2, lty = 2)
atv <- c((-limn - 1) / 2, 0, (limn + 1) / 2)
for (i in 1:3) {
  for (j in 1:3) {
    labc <- paste(prettyNum(ntab[i, j], big.mark = ","), " (", formatC(100 * ptab[i, j], format = "f", digits = 1), "%)", sep = "")
    legend(x = atv[i], y = atv[j], legend = labc, text.col = 1, bg = adjustcolor("white", .75), xjust = .5, cex = legcex)
  }
}
numv <- formatC(c(fdrci[[1]], fdrci[[2]]) * 100, format = "f", digits = 1)
fdrlab <- paste("Estimated Fsr = ", numv[1], "% (", numv[2], "% - ", numv[3], "%)", sep = "")
mtext(side = 3, text = fdrlab, line = .5, cex = 1)
dev.off()
if (control$output_to_dropbox) {
  file.copy(from = paste(control$figure_dir, "/", fnamc, sep = ""),
            to = paste(control$dropbox_figure_dir, "/", fnamc, sep = ""), overwrite = TRUE)
} else {
  file.rename(from = paste(control$figure_dir, "/", fnamc, sep = ""), 
              to = file.path(control$figure_dir, "Figure_S6.jpg"))
}


#global concordance/discordance
hethom.concordant.n <- tabhethom["-1", "-1"] + tabhethom["1", "1"]
hethom.discordant.n <- tabhethom["-1", "1"] + tabhethom["1", "-1"]
save.num <- c("hethom.concordant.n", "hethom.discordant.n")
for (numc in save.num) {
  write.table(prettyNum(eval(as.name(numc)), big.mark = ","), file = paste(control$dropbox_text_numbers_dir, "/", numc, ".txt", sep = ""), 
              col.names = F, row.names = F, quote = F)
}

#het/hom hit rate
hom.hitrate <- mean(hom.z.sign != 0, na.rm = T)
het.hitrate <- mean(het.z.sign != 0, na.rm = T)
save.num1 <- c("hom.hitrate", "het.hitrate")
for (numc in save.num1) {
  write.table(formatC(eval(as.name(numc)) * 100, format = "f", digits = 1), 
              file = paste(control$dropbox_text_numbers_dir, "/", numc, ".txt", sep = ""), 
              col.names = F, row.names = F, quote = F)
}


