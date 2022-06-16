################################################
# This script generates Figure 1
################################################

Data_all <- readRDS(control$Data_all_file)
d1 <- Data_all$impc$figure_1_data$meas[[1]]
d2 <- Data_all$impc$figure_1_data$meas[[2]]
gendo <-  Data_all$impc$figure_1_data$geno["do"]
genup <-  Data_all$impc$figure_1_data$geno["up"]
phnam1use <-  Data_all$impc$figure_1_data$ph_names[1]
phnam2use <-  Data_all$impc$figure_1_data$ph_names[2]

resimp <- readRDS(file = control$file.resimp)
ph_use <- Data_all$impc$phmap[match(c(phnam1use, phnam2use), Data_all$impc$phmap$nam), "ph"]

mean(abs(resl.comp$impc_MVphen_nSig_1_K_20$mn) < abs(resl.comp$uv$mn), na.rm = T)
mean(abs(resl.comp$impc_MVphen_nSig_1_K_20$sd) <= abs(resl.comp$uv$sd), na.rm = T)

#######################################################
#Version for paper
graphics.off()
namc <- paste0("Figure_1.jpg")
jpeg(file = paste(control$figure_dir, "/", namc, sep = ""), 12, 7, units = "in", res = 500)
par(mfrow = c(1, 2), oma = c(1, 2, 1, 2), mar = c(3, 5, 3, 5))
for(i in 1:2){
  dc <- list(d1, d2)[[i]]
  dc$yuse <- dc$yuse / dc$scale_factor[1]
  phnamc <- switch(i, phnam1use, phnam2use)
  plot(dc$day, dc$yuse, cex = .5, ylab = "", xlab = "", xaxt = "n", yaxt = "n", col = "grey")
  labcex <- 1.2
  titcex <- 1.25
  ltyuse <- 3
  if(i == 1)
    mtext(side = 3, line = 1, cex = titcex, text = expression(paste("(a) Triglycerides, phenotype p")))
  if(i == 2)
    mtext(side = 3, line = 1, cex = titcex, text = expression(paste("(b) Body fat percentage, phenotype ", tilde(p))))
  mtext(side = 1, line = 1, cex = labcex, text = "Date")
  mtext(side = 2, line = 1, cex = labcex, text = "y", las = 2)
  points(x = dc[dc$geno == genup, "day"], y = dc[dc$geno == genup, "yuse"], col = "red", pch = 19)
  points(x = dc[dc$geno == gendo, "day"], y = dc[dc$geno == gendo, "yuse"], col = "blue", pch = 19)
  mnv <- list(wt = median(as.numeric(dc[dc$geno == "0", "yuse"]), na.rm = T),
              ko1 = mean(as.numeric(dc[dc$geno == genup, "yuse"]), na.rm = T),
              ko2 = mean(as.numeric(dc[dc$geno == gendo, "yuse"]), na.rm = T))
  spl_fit <- smooth.spline(x = dc[dc$geno == "0", "day"], y = dc[dc$geno == "0", "yuse"], spar = .75)
  lines(spl_fit, lwd = 1)
  dayun <- sort(unique(dc[, "day"]))
  
  mean_day_genup <- dayun[match(T, max(dc[dc$geno == genup, "day"]) < dayun)]
  mean_day_gendo <- dayun[match(T, max(dc[dc$geno == gendo, "day"]) < dayun)]
  do_uv_mn <- predict(spl_fit, x = mean_day_gendo)$y + resl.comp$uv$mn[gendo, ph_use[i]]
  do_uv_sd <- resl.comp$uv$sd[gendo, ph_use[i]]
  do_mv_mn <- predict(spl_fit, x = mean_day_gendo)$y + resl.comp$impc_MVphen_nSig_1_K_20$mn[gendo, ph_use[i]]
  do_mv_sd <- resl.comp$impc_MVphen_nSig_1_K_20$sd[gendo, ph_use[i]]
  up_uv_mn <- predict(spl_fit, x = mean_day_genup)$y + resl.comp$uv$mn[genup, ph_use[i]]
  up_uv_sd <- resl.comp$uv$sd[genup, ph_use[i]]
  up_mv_mn <- predict(spl_fit, x = mean_day_genup)$y + resl.comp$impc_MVphen_nSig_1_K_20$mn[genup, ph_use[i]]
  up_mv_sd <- resl.comp$impc_MVphen_nSig_1_K_20$sd[genup, ph_use[i]]
  names(resl.comp)
  eps = diff(range(dayun)) * .04
  eps2 = diff(range(dayun)) * .01
  lwd_whisk <- 2
  lines(x = rep(mean_day_gendo + eps, 2), y = do_uv_mn + 2 * do_uv_sd * c(-1, 1), col = "black", lwd = lwd_whisk)
  for (yadd in c(-1, 1)) {
    lines(x = mean_day_gendo + eps + eps2 * c(-1, 1), y = rep(do_uv_mn + 2 * do_uv_sd * yadd, 2), col = "black", lwd = lwd_whisk)
  }
  points(x = mean_day_gendo + eps, y = do_uv_mn, col = "black", pch = 22, bg = "white")
  lines(x = rep(mean_day_gendo + 2 * eps, 2), y = do_mv_mn + 2 * do_mv_sd * c(-1, 1), col = "black", lwd = lwd_whisk)
  for (yadd in c(-1, 1)) {
    lines(x = mean_day_gendo + 2 * eps + eps2 * c(-1, 1), y = rep(do_mv_mn + 2 * do_mv_sd * yadd, 2), col = "black", lwd = lwd_whisk)
  }
  points(x = mean_day_gendo + 2 * eps, y = do_mv_mn, col = "black", pch = 15)
  lines(x = rep(mean_day_genup + eps, 2), y = up_uv_mn + 2 * up_uv_sd * c(-1, 1), col = "black", lwd = lwd_whisk)
  for (yadd in c(-1, 1)) {
    lines(x = mean_day_genup + eps + eps2 * c(-1, 1), y = rep(up_uv_mn + 2 * up_uv_sd * yadd, 2), col = "black", lwd = lwd_whisk)
  }
  points(x = mean_day_genup + eps, y = up_uv_mn, col = "black", pch = 22, bg = "white")
  lines(x = rep(mean_day_genup + 2 * eps, 2), y = up_mv_mn + 2 * up_mv_sd * c(-1, 1), col = "black", lwd = lwd_whisk)
  for (yadd in c(-1, 1)) {
    lines(x = mean_day_genup + 2 * eps + eps2 * c(-1, 1), y = rep(up_mv_mn + 2 * up_mv_sd * yadd, 2), col = "black", lwd = lwd_whisk)
  }
  points(x = mean_day_genup + 2 * eps, y = up_mv_mn, col = "black", pch = 15)
  ltyuse_uv <- 2
  ltyuse_mv <- 3
  linat <- .5
  newlabcex <- 1
}
yat <- .7
par(fig = c(.4, .6, yat, yat + .25), 
    mar = c(0, 0, 0, 0),
    new = T)
plot(0, axes = F)
legend(x = "center", legend = c("WT mice", "KO gene g mice", 
                                expression(paste("KO gene ", tilde(g), " mice")),
                                "UV model estimate",
                                "MV model estimate"), 
         pch = c(1, 19, 19, 0, 15), col = c("grey", "red", "blue", 1, 1))
dev.off()


