################################################
# This script generates Figure S2
################################################

tpl <- Data_all$impc$figure_S2_data

resimp <- readRDS(file = control$file.resimp)
resimp <- resimp[!grepl("fac", resimp$ph), ]

omitl <- list(list(cen = c("Bcm"), proc = c("Auditory Brain Stem Response"), omit = T),
              list(cen = c("Gmc"), proc = c("Eye Morphology"), omit = T),
              list(cen = c("Gmc"), proc = c("Echo"), omit = T),
              list(cen = c("Gmc"), proc = c("FACS"), omit = T),
              list(cen = c("Ics"), proc = c("Hematology"), omit = T),
              list(cen = c("Ics"), proc = c("Open Field"), omit = T),
              list(cen = c("J"), proc = c("Clinical Chemistry"), omit = T),
              list(cen = c("J"), proc = c("FACS"), omit = T),
              list(cen = c("J"), proc = c("Hematology"), omit = T),
              list(cen = c("J"), proc = c("Heart Weight"), omit = T),
              list(cen = c("J"), proc = c("Electroconvulsive Threshold Testing"), omit = T),
              list(cen = c("Tcp"), proc = c("FACS"), omit = T),
              list(cen = c("Tcp"), proc = c("Combined SHIRPA and Dysmorphology"), omit = T),
              list(cen = c("Tcp"), proc = c("Acoustic Startle and Pre-pulse Inhibition (PPI)"), omit = T),
              list(cen = c("Tcp"), proc = c("Hematology"), omit = T),
              list(cen = c("Tcp"), proc = c("Open Field"), omit = T),
              list(cen = c("Wtsi"), proc = c("Grip Strength"), omit = T))


##########################################################
# Producing Figure S2 from MV paper 
##############################################
fnamc <- "uv_qc_plot_for_paper.jpg"
jpeg(paste(control$figure_dir, "/", fnamc, sep = ""), 12, 11, units = "in", res = 800)
par(mfrow = c(1, 1), mar = c(6, 14, 3, 1))

image(x = 1:nrow(tpl), y = 1:ncol(tpl), z = tpl, col = control$heat_col_palette, zlim = c(-1, 1), xlab = "", ylab = "", xaxt = "n", yaxt = "n")

procv <- Data_all$impc$phmap[match(colnames(tpl), Data_all$impc$phmap$ph), "procnam"] 
procun = unique(procv)
ats = sapply(procun, function(x) mean(which(procv == x)))
linats = sapply(procun, function(x) match(x, procv) - .5)
axis(side = 2, at = ats, labels = procun, las = 2, cex.axis = .7)
abline(h = linats)
cenv = resimp[match(rownames(tpl), resimp$geno), "cenlong"]
cennamun <- unique(cenv)
cenats = sapply(cennamun, function(x) mean(which(cenv == x)))
cenlinats = sapply(cennamun, function(x) match(x, cenv) - .5)
axis(side = 1, at = cenats, labels = cennamun, las = 2, cex.axis = .7)
abline(v = cenlinats)
for(j in 1:length(omitl)){
  cenc <- omitl[[j]]$cen
  procc <- omitl[[j]]$proc
  omitc <- omitl[[j]]$omit
  colc <- ifelse(omitc, "red", NA)
  cenran <- range(which(cenv == cenc)) + c(-1, 1) * .5
  procran <- range(which(procv == procc)) + c(-1, 1) * .5
  polygon(x = cenran[c(1, 1, 2, 2, 1)], y = procran[c(1, 2, 2, 1, 1)], border = colc, lwd = 2, col = NA)
}
dev.off()
if (control$output_to_dropbox) {
  file.copy(from = paste(control$figure_dir, "/", fnamc, sep = ""),
            to = paste(control$dropbox_figure_dir, "/", fnamc, sep = ""), overwrite = TRUE)
} else {
  file.rename(from = paste(control$figure_dir, "/", fnamc, sep = ""), 
              to = file.path(control$figure_dir, "Figure_S2.jpg"))
}


