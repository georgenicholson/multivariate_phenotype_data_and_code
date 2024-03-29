# control = control
# resl = resl.err.rates[!names(resl.err.rates) %in% c("uv", "uv.ss")]
# err.rate.meth = "lfsr"
# test.stat = "lfsr"
# linemap = Data_all$impc$linemap
# reflinemap = Data_all$impc$reflinemap
# phmap = Data_all$impc$phmap
# cenmap = Data_all$impc$cenmap
# Yhat = Data_all$impc$Y_raw
# use.upper.fp.est = F
# control.level = c("line.fdr", "line.fwer", "phcen.fdr")[1]
# err.thresh = .05
# p.complete.null.true = 1
# p.test.null.true = 1




########################################
#Need to fix this function to have arguments for all variables used
# err.rate.meth = "perm"; calibrate = F; sep.imp.thresh = T; test.stat = c("z", "lfsr")[1]
err.rate.control <- function(control, 
                              resl, 
                              err.rate.meth = "perm", 
                              test.stat = c("z", "lfsr")[1], 
                              linemap, 
                              reflinemap,
                              phmap = NULL, 
                              cenmap,
                              Yhat, 
                              control.level = c("line.fdr", "line.fwer", "phcen.fdr")[1],
                              err.thresh = .05, 
                              p.complete.null.true = 1, 
                              p.test.null.true = 1){
  methv <- names(resl)
  require(Hmisc)
  fdr.th.seq.finemesh <- seq(0, 10, by = .01)
  ####################################################################################
  #Create resimp, update fields in resmat, copy fields from resmat to resimp
  ph.use.inc.fac <- unique(unlist(sapply(resl, function(x) colnames(x$mn), simplify = F)))
  lines.all <- linemap$geno#unique(unlist(sapply(resl, function(x) rownames(x$mn), simplify = F)))
  resimp <- data.frame(ph = rep(ph.use.inc.fac, length(lines.all)), geno = rep(lines.all, each = length(ph.use.inc.fac)))
  resimp[, "cen"] <- linemap[match(resimp$geno, linemap$geno), "cen"]
  resimp$line.type <- linemap[match(resimp$geno, linemap$geno), "line.type"]
  resimp$testid <- paste(resimp$cen, resimp$ph, resimp$geno, sep = "_")
  resimp$cenlong <- cenmap[match(resimp$cen, cenmap$cen), "nam"]
  resimp$cenphen <- paste(resimp$cenlong, resimp$ph, sep = "_")
  resimp$imputed <- NA
  faclog <- !(resimp$ph %in% colnames(Yhat))
  resimp$imputed[faclog] <- F
  resimp$imputed[!faclog] <- is.na(Yhat[cbind(resimp$geno, resimp$ph)[!faclog, ]])
  genospl <- strsplit(resimp$geno, spl = "_")
  resimp$geno.sh <- sapply(genospl, function(x) x[1])
  resimp$zyg <- sapply(genospl, function(x) x[2])
  controls.for.fitting <- controls.for.th.select <- "negCon"
  truekos.for.fitting <- truekos.for.th.select <- "trueMut"
  linemap$cenlong <- cenmap[match(linemap$cen, cenmap$cen), "nam"]
  if(err.rate.meth == "perm"){
    ###############################################
    #Calibration and Threshold selection
    cenphenun <- unique(resimp$cenphen[resimp$line.type == truekos.for.th.select])
    cennamun.sub <- unique(resimp$cenlong[which(resimp$line.type == truekos.for.th.select)])
    thdat <- resimp[match(cenphenun, resimp$cenphen), c("cenphen", "cenlong", "ph")]
    for(methc in methv){#methc <- "mash"#
      ###########################
      #Thresh select
      cpnam <- paste(methc, "cenphen.th", sep = ".")
      cnam <- paste(methc, "cen.th", sep = ".")
      gnam <- paste(methc, "glob.th", sep = ".")
      fnam <- paste(methc, "th.final", sep = ".")
      tnam <- paste(methc, "t", sep = ".")
      signsignam <- paste(methc, "perm.signsig", sep = ".")
      if(test.stat == "z"){
        test.stat.mat <- resl[[methc]]$mn / resl[[methc]]$sd
      }
      if(test.stat == "lfsr"){
        test.stat.mat <- qnorm(.5 + (1 - resl[[methc]]$lfsr) / 2) * sign(resl[[methc]]$mn)
        # test.stat.mat <- (1 / resl[[methc]]$lfsr) * sign(resl[[methc]]$mn)
      }
      resimp[, tnam] <- NA
      ph.availc <- which(resimp$ph %in% colnames(test.stat.mat))
      resimp[ph.availc, tnam] <- test.stat.mat[cbind(resimp$geno[ph.availc], resimp$ph[ph.availc])]
      if(control.level %in% c("line.fdr", "line.fwer")){
        zmat.null <- test.stat.mat[linemap[linemap$line.type %in% controls.for.th.select, "geno"], ]
        zmat.true <- test.stat.mat[linemap[linemap$line.type %in% truekos.for.th.select, "geno"], ]
        zmat.null <- zmat.null[rowMeans(is.na(zmat.null)) < 1, ]
        zmat.true <- zmat.true[rowMeans(is.na(zmat.true)) < 1, ]
        tmax0 <- apply(abs(zmat.null), 1, function(v) max(v, na.rm = T))
        tmax1 <- apply(abs(zmat.true), 1, function(v) max(v, na.rm = T))
        prob.disc.false <- colMeans(outer(tmax0, fdr.th.seq.finemesh, '>'), na.rm = T)
        prob.disc <- colMeans(outer(tmax1, fdr.th.seq.finemesh, '>'), na.rm = T)
        line.fdr.out <- switch(control.level, 
                               line.fdr = p.complete.null.true * prob.disc.false / prob.disc, 
                               line.fwer = prob.disc.false)
        resimp[, gnam] <- fdr.th.seq.finemesh[match(T, line.fdr.out < err.thresh)]
        resimp[, cnam] <- NA
        for(cenc in cennamun.sub){
          zmat.null <- test.stat.mat[linemap[linemap$line.type %in% controls.for.th.select & linemap$cenlong == cenc, "geno"], ]
          zmat.true <- test.stat.mat[linemap[linemap$line.type %in% truekos.for.th.select & linemap$cenlong == cenc, "geno"], ]
          zmat.null <- zmat.null[rowMeans(is.na(zmat.null)) < 1, ]
          zmat.true <- zmat.true[rowMeans(is.na(zmat.true)) < 1, ]
          tmax0 <- apply(abs(zmat.null), 1, function(v) max(v, na.rm = T))
          tmax1 <- apply(abs(zmat.true), 1, function(v) max(v, na.rm = T))
          prob.disc <- colMeans(outer(tmax1, fdr.th.seq.finemesh, '>'), na.rm = T)
          mat.neg <- outer(tmax0, fdr.th.seq.finemesh, '>')
          prob.disc.false <- colMeans(mat.neg, na.rm = T)
          line.fdr.out <- switch(control.level, 
                                 line.fdr = p.complete.null.true * prob.disc.false / prob.disc, 
                                 line.fwer = prob.disc.false)
          resimp[which(resimp$cenlong == cenc), cnam] <- fdr.th.seq.finemesh[match(T, line.fdr.out < err.thresh)]
        }  
        # resimp[, fnam] <- pmax(resimp[, gnam], resimp[, cnam])
        resimp[, fnam] <- resimp[, gnam]
        resimp[, fnam][is.na(resimp[, fnam]) & !is.na(resimp[, tnam])] <- Inf
      }
      resimp[, signsignam] <- (abs(resimp[, tnam]) > resimp[, fnam]) * sign(resimp[, tnam])
      resl[[methc]]$signsig <- resl[[methc]]$t <- resl[[methc]]$th <- resl[[methc]]$mn
      resl[[methc]]$signsig[] <- resl[[methc]]$t[] <- resl[[methc]]$th[] <- NA
      ph.availc <- which(resimp$ph %in% colnames(resl[[methc]]$signsig))
      resl[[methc]]$signsig[cbind(resimp$geno[ph.availc], resimp$ph[ph.availc])] <- resimp[ph.availc, signsignam]
      resl[[methc]]$t[cbind(resimp$geno[ph.availc], resimp$ph[ph.availc])] <- resimp[ph.availc, tnam]
      resl[[methc]]$th[cbind(resimp$geno[ph.availc], resimp$ph[ph.availc])] <- resimp[ph.availc, fnam]
    }
  }
  if(err.rate.meth %in% c("lfdr", "lfsr")){
    for(methc in methv){
      fnam <- paste(methc, "th.final", sep = ".")
      tnam <- paste(methc, "t", sep = ".")
      signsignam <- paste0(methc, ".", err.rate.meth, ".signsig")
      # following line req'd as factors available only if methc == varimax
      which_resimp_ph_in_resl <- which(resimp$ph %in% colnames(resl[[methc]][[err.rate.meth]]))
      resimp[which_resimp_ph_in_resl, tnam] <- 1 / resl[[methc]][[err.rate.meth]][cbind(resimp$geno, resimp$ph)[which_resimp_ph_in_resl, ]] * 
                                                sign(resl[[methc]]$mn[cbind(resimp$geno, resimp$ph)[which_resimp_ph_in_resl, ]])
      resimp[which_resimp_ph_in_resl, fnam] <- 1 / err.thresh
      resimp[which_resimp_ph_in_resl, signsignam] <- (resimp[which_resimp_ph_in_resl, tnam] > resimp[which_resimp_ph_in_resl, fnam]) * 
                                              sign(resl[[methc]]$mn[cbind(resimp$geno, resimp$ph)[which_resimp_ph_in_resl, ]])
      resl[[methc]]$signsig <- resl[[methc]]$t <- resl[[methc]]$th <- resl[[methc]]$mn
      resl[[methc]]$signsig[] <- resl[[methc]]$t[] <- resl[[methc]]$th[] <- NA
      ph.availc <- which(resimp$ph %in% colnames(resl[[methc]]$signsig))
      resl[[methc]]$signsig[cbind(resimp$geno[ph.availc], resimp$ph[ph.availc])] <- resimp[ph.availc, signsignam]
      resl[[methc]]$t[cbind(resimp$geno[ph.availc], resimp$ph[ph.availc])] <- resimp[ph.availc, tnam]
      resl[[methc]]$th[cbind(resimp$geno[ph.availc], resimp$ph[ph.availc])] <- resimp[ph.availc, fnam]
    }
  }

  
  resimp$geno[!resimp$geno %in% rownames(resl[[methc]][[err.rate.meth]])]
  resimp$ph[!resimp$ph %in% colnames(resl[[methc]][[err.rate.meth]])]
  
  ########################################################
  #Extract reference lines results
  # reflinemap <- read.csv(paste0(base.dir, "/data_in/reference_line_genotypeIds.csv"))
  linun <- unique(reflinemap$gene_symbol)
  phun <- unique(resimp$ph)
  for(methc in methv){
    ssl <- tl <- thl <- list()
    for(linc in linun){#linc <- linun[2]#
      for(zyg in 0:1){#zyg <- 0#
        namc <- paste(linc, zyg, sep = "_")
        gen.shv <- reflinemap[reflinemap$gene_symbol == linc, "genotype_id"]
        resc <- resimp[resimp$zyg == zyg & resimp$line.type %in% truekos.for.fitting, ]
        genv <- c(na.omit(resc[match(gen.shv, resc$geno.sh), "geno"]))
        if(length(genv) <= 1)
          next
        for(genc in genv){
          resimpc <- resimp[resimp$geno == genc, ]
          ssl[[namc]] <- cbind(ssl[[namc]], resimpc[match(phun, resimpc$ph), paste(methc, err.rate.meth, "signsig", sep = ".")])
          tl[[namc]] <- cbind(tl[[namc]], resimpc[match(phun, resimpc$ph), paste(methc, "t", sep = ".")])
          thl[[namc]] <- cbind(thl[[namc]], resimpc[match(phun, resimpc$ph), paste(methc, "th.final", sep = ".")])
        }
        dimnames(ssl[[namc]]) <- dimnames(tl[[namc]]) <- dimnames(thl[[namc]]) <- list(phun, genv)
      }
    }
    resl[[methc]]$ref.lines <- list(ssl = ssl, tl = tl, thl = thl)
  }
  
  ########################################################
  #Extract het/hom matched results
  resimp1 <- resimp[resimp$line.type %in% truekos.for.fitting, ]
  resimp1$geno.sh.ph = paste(resimp1$geno.sh, resimp1$ph, sep = "_")
  tabc = table(resimp1$geno.sh.ph)
  geno.sh.ph.hethom = names(tabc)[tabc == 2]
  resimp.homfirst = resimp1[order(-as.numeric(resimp1$zyg)), ]
  resimp.hetfirst = resimp1[order(as.numeric(resimp1$zyg)), ]
  resimp.hom = resimp.homfirst[match(geno.sh.ph.hethom, resimp.homfirst$geno.sh.ph), ]
  resimp.het = resimp.hetfirst[match(geno.sh.ph.hethom, resimp.hetfirst$geno.sh.ph), ]
  for(methc in methv){
    ssl <- tl <- thl <- ph_genl <- list()
    ssl$het.hom <- cbind(resimp.hom[, paste(methc, err.rate.meth, "signsig", sep = ".")], 
                         resimp.het[, paste(methc, err.rate.meth, "signsig", sep = ".")])
    tl$het.hom <- cbind(resimp.hom[, paste(methc, "t", sep = ".")], 
                         resimp.het[, paste(methc, "t", sep = ".")])
    thl$het.hom <- cbind(resimp.hom[, paste(methc, "th.final", sep = ".")], 
                         resimp.het[, paste(methc, "th.final", sep = ".")])
    geno.sh.ph.hethom.just.geno.sh <- sapply(strsplit(geno.sh.ph.hethom, split = "_"), function(v) v[1])
    geno.sh.ph.hethom.just.phen <- sapply(strsplit(geno.sh.ph.hethom, split = "_"), function(v) paste(v[2:length(v)], collapse = "_"))
    ph_genl$het.hom <- cbind(paste0(geno.sh.ph.hethom.just.phen, "_", geno.sh.ph.hethom.just.geno.sh, "_1"), 
                             paste0(geno.sh.ph.hethom.just.phen, "_", geno.sh.ph.hethom.just.geno.sh, "_0"))
    dimnames(ssl$het.hom) <- dimnames(tl$het.hom) <- dimnames(thl$het.hom) <- dimnames(ph_genl$het.hom) <- list(geno.sh.ph.hethom, 1:0)
    # dimnames(ssl[[namc]]) <- dimnames(tl[[namc]]) <- dimnames(thl[[namc]]) <- list(phun, genv)
    resl[[methc]]$het.hom <- list(ssl = ssl, tl = tl, thl = thl, ph_genl = ph_genl)
  }

  ##################################################################
  #Analyse reference lines and het/homs results for hit concordance
  resimp$ph_geno <- paste(resimp$ph, resimp$geno, sep = "_")
  comp.typev <- c("ref.lines", "het.hom")
  for(comp.type in comp.typev){#comp.type <- "ref.lines"#
    for(methc in methv){#methc <- "impc_mash_nSig_1"#
      matl <- tabl <- list()
      for(namc in names(resl[[methc]][[comp.type]]$ssl)){
        # tabpost <- tabpostcomp <- tabpostimp <- resl[[methc]][[comp.type]]$ssl[[namc]]
        # tabpre <- resl$uv[[comp.type]]$ssl[[namc]]
        tabpost <- tabpostcomp <- tabpostimp <- resl[[methc]][[comp.type]]$tl[[namc]] / resl[[methc]][[comp.type]]$thl[[namc]]
        # tabpre <- resl$uv[[comp.type]]$tl[[namc]] / resl$uv[[comp.type]]$thl[[namc]]
        # is.na(Yhat[cbind(resimp$geno, resimp$ph)[!faclog, ]])
        if(comp.type == "ref.lines")
          ph_geno.tabpost <- outer(rownames(tabpost), colnames(tabpost), function(x, y) paste(x, y, sep = "_"))
        if(comp.type == "het.hom")
          ph_geno.tabpost <- resl[[methc]][[comp.type]]$ph_genl$het.hom
        implog <- resimp[match(ph_geno.tabpost, resimp$ph_geno), "imputed"]
        tabpostcomp[implog] <- NA
        tabpostimp[!implog] <- NA
        for(i in 1:(ncol(tabpost) - 1)){
          for(j in (i + 1):ncol(tabpost)){
            matl$post <- rbind(matl$post, tabpost[, c(i, j)])
            matl$postcomp <- rbind(matl$postcomp, tabpostcomp[, c(i, j)])
            matl$postimp <- rbind(matl$postimp, tabpostimp[, c(i, j)])
            matl$postimpcomp <- rbind(matl$postimpcomp, cbind(tabpostimp[, i], tabpostcomp[, j]), 
                                      cbind(tabpostimp[, j], tabpostcomp[, i]))
          }
        }
      }
      matl <- lapply(matl, function(M) M[rowSums(is.na(M)) == 0, , drop = F])
      matl.sig <- lapply(matl, function(M){M[abs(M) < 1] <- 0; return(sign(M))})
      for(ty in names(matl)){
        if(NROW(matl[[ty]]) != 0){
          M <- matl.sig[[ty]]
          Mtab <- table(M[, 1], M[, 2])
          namv <- seq(-1, 1, by = 1); 
          Mout <- Mtab[match(namv, rownames(Mtab)), match(namv, colnames(Mtab))]; 
          Mout[is.na(Mout)] <- 0; 
          dimnames(Mout) <- list(namv, namv)
          tabl[[ty]] <- Mout
        } else {
          tabl[[ty]] <- matrix(0, 3, 3, dimnames = rep(list((-1):1), 2))
        }
      }
      fdr.ests <- lapply(tabl, fdr.est.tab)
      resl[[methc]][[comp.type]][c("matl", "tabl", "fdr.ests")] <- list(matl, tabl, fdr.ests)
    }
  }  
  

  resimp1 <- resimp[resimp$line.type %in% truekos.for.fitting, ]
  resimp0 <- resimp[resimp$line.type %in% controls.for.fitting, ]
  all.null.lines <- unique(resimp0$geno)
  all.true.lines <- unique(resimp1$geno)
  B <- 1000
  bootmat.n.null.test <- bootmat.n.null.hit <- bootmat <- matrix(0, length(all.null.lines), B)
  bootmat.n.true.test.imp <- bootmat.n.true.hit.imp <- 
    bootmat.n.true.test.nonimp <- bootmat.n.true.hit.nonimp <- bootmat.true <- matrix(0, length(all.true.lines), B)
  set.seed(1)
  bootmat[] <- sample(all.null.lines, length(all.null.lines) * B, replace = T)
  bootmat.true[] <- sample(all.true.lines, length(all.true.lines) * B, replace = T)
  namv <- paste(methv, err.rate.meth, "signsig", sep = ".")


  null.line.nhit.mat <- as.matrix(sapply(all.null.lines, function(genoc) 
    colSums(resimp0[which(resimp0$geno == genoc), namv, drop = F] != 0, na.rm = T)))
  null.line.ntest.mat <- as.matrix(sapply(all.null.lines, function(genoc) 
    colSums(!is.na(resimp0[which(resimp0$geno == genoc), namv, drop = F]))))
  true.line.nhit.mat.imp <- as.matrix(sapply(all.true.lines, function(genoc) 
    colSums(resimp1[which(resimp1$geno == genoc & resimp1$imputed), namv, drop = F] != 0, na.rm = T)))
  true.line.ntest.mat.imp <- as.matrix(sapply(all.true.lines, function(genoc) 
    colSums(!is.na(resimp1[which(resimp1$geno == genoc & resimp1$imputed), namv, drop = F]))))
  true.line.nhit.mat.nonimp <- as.matrix(sapply(all.true.lines, function(genoc) 
    colSums(resimp1[which(resimp1$geno == genoc & !resimp1$imputed), namv, drop = F] != 0, na.rm = T)))
  true.line.ntest.mat.nonimp <- as.matrix(sapply(all.true.lines, function(genoc) 
    colSums(!is.na(resimp1[which(resimp1$geno == genoc & !resimp1$imputed), namv, drop = F]))))
  if (length(namv) > 1) {
    null.line.nhit.mat <- t(null.line.nhit.mat)
    null.line.ntest.mat <- t(null.line.ntest.mat)
    true.line.nhit.mat.imp <- t(true.line.nhit.mat.imp)
    true.line.ntest.mat.imp <- t(true.line.ntest.mat.imp)
    true.line.nhit.mat.nonimp <- t(true.line.nhit.mat.nonimp)
    true.line.ntest.mat.nonimp <- t(true.line.ntest.mat.nonimp)
  }
  dimnames(true.line.nhit.mat.imp) <- dimnames(true.line.ntest.mat.imp) <- 
    dimnames(true.line.nhit.mat.nonimp) <- dimnames(true.line.ntest.mat.nonimp) <- list(all.true.lines, methv)
  dimnames(null.line.nhit.mat) <- dimnames(null.line.ntest.mat) <- list(all.null.lines, methv)
  restab <- data.frame()
  ci.al <- .05
  for(methc in methv){#methc = methv[2]#
    namc <- paste(methc, err.rate.meth, "signsig", sep = ".")
    hit.rate.imp <- mean(resimp1[resimp1$imputed, namc] != 0, na.rm = T)
    hit.rate.nonimp <- mean(resimp1[!resimp1$imputed, namc] != 0, na.rm = T)
    bootmat.n.true.hit.imp[] <- true.line.nhit.mat.imp[match(c(bootmat.true), rownames(true.line.nhit.mat.imp)), methc]
    bootmat.n.true.test.imp[] <- true.line.ntest.mat.imp[match(c(bootmat.true), rownames(true.line.ntest.mat.imp)), methc]
    bootv.imp <- colSums(bootmat.n.true.hit.imp, na.rm = T) / colSums(bootmat.n.true.test.imp, na.rm = T)
    bootmat.n.true.hit.nonimp[] <- true.line.nhit.mat.nonimp[match(c(bootmat.true), rownames(true.line.nhit.mat.nonimp)), methc]
    bootmat.n.true.test.nonimp[] <- true.line.ntest.mat.nonimp[match(c(bootmat.true), rownames(true.line.ntest.mat.nonimp)), methc]
    bootv.nonimp <- colSums(bootmat.n.true.hit.nonimp, na.rm = T) / colSums(bootmat.n.true.test.nonimp, na.rm = T)
    hit.rate.imp.ci <- c(hit.rate.imp, quantile(bootv.imp, c(ci.al, 1 - ci.al), na.rm = T))
    hit.rate.nonimp.ci <- c(hit.rate.nonimp, quantile(bootv.nonimp, c(ci.al, 1 - ci.al), na.rm = T))
    
    line.fpr.x <- length(unique(resimp0$geno[resimp0[, namc] != 0]))
    line.fpr.n <- length(unique(resimp0$geno))
    line.fpr.ci <- Hmisc::binconf(x = line.fpr.x, n = line.fpr.n, method = "exact")
    line.fpr.est <- line.fpr.ci[, "PointEst"]
    line.pr.est <- length(unique(resimp1$geno[resimp1[, namc] != 0])) / length(unique(resimp1$geno))
    line.fdr.ci <- pmin(1, p.complete.null.true * line.fpr.ci / line.pr.est)
    line.fdr.est <- line.fdr.ci[1]
    line.fdr.est[line.fpr.est == 0 & line.pr.est == 0] <- 0
    fpr.est <- mean(resimp0[, namc] != 0, na.rm = T)
    bootmat.n.null.hit[] <- null.line.nhit.mat[match(c(bootmat), rownames(null.line.nhit.mat)), methc]
    bootmat.n.null.test[] <- null.line.ntest.mat[match(c(bootmat), rownames(null.line.ntest.mat)), methc]
    bootv <- colSums(bootmat.n.null.hit, na.rm = T) / colSums(bootmat.n.null.test, na.rm = T)
    pr.est <- mean(resimp1[, namc] != 0, na.rm = T)
    fpr.ci <- c(fpr.est, quantile(bootv, c(ci.al, 1 - ci.al), na.rm = T))
    fdr.est <- min(c(1, p.test.null.true * fpr.est / pr.est))
    fdr.est[fpr.est == 0 & pr.est == 0] <- 0
    fpr.est.nonimp <- mean(resimp0[!resimp0$imputed, namc] != 0, na.rm = T)
    pr.est.nonimp <- mean(resimp1[!resimp1$imputed, namc] != 0, na.rm = T)
    fdr.est.nonimp <- min(c(1, p.test.null.true * fpr.est.nonimp / pr.est.nonimp))
    fdr.ci <- pmin(1, p.test.null.true * fpr.ci / pr.est)
    fdr.est.nonimp[fpr.est.nonimp == 0 & pr.est.nonimp == 0] <- 0
    fpr.est.imp <- mean(resimp0[resimp0$imputed, namc] != 0, na.rm = T)
    pr.est.imp <- mean(resimp1[resimp1$imputed, namc] != 0, na.rm = T)
    fdr.est.imp <- min(c(1, p.test.null.true * fpr.est.imp / pr.est.imp))
    fdr.est.imp[fpr.est.imp == 0 & pr.est.imp == 0] <- 0
    restab.add <- data.frame(meth = methc, 
                             err.rate.meth = err.rate.meth,
                             test.stat = test.stat, 
                             mvhitimp = hit.rate.imp, 
                             hit.rate.imp.ci.l = hit.rate.imp.ci[2], 
                             hit.rate.imp.ci.u = hit.rate.imp.ci[3], 
                             mvhitnonimp = hit.rate.nonimp, 
                             hit.rate.nonimp.ci.l = hit.rate.nonimp.ci[2], 
                             hit.rate.nonimp.ci.u = hit.rate.nonimp.ci[3], 
                             line.fdr.est = line.fdr.est, 
                             line.fdr.ci.l = line.fdr.ci[2], 
                             line.fdr.ci.u = line.fdr.ci[3], 
                             fdr.est = fdr.est, 
                             fdr.ci.l = fdr.ci[2], 
                             fdr.ci.u = fdr.ci[3], 
                             fdr.est.imp = fdr.est.imp, 
                             fdr.est.nonimp = fdr.est.nonimp)
    for(comp.type in comp.typev){
      fdr.ests.c <- resl[[methc]][[comp.type]]$fdr.ests
      restab.add[, paste0(comp.type, ".", names(fdr.ests.c))] <- sapply(fdr.ests.c, function(x) x$est)
      restab.add[, paste0(comp.type, ".", names(fdr.ests.c), ".l")] <- sapply(fdr.ests.c, function(x) x$ci[1])
      restab.add[, paste0(comp.type, ".", names(fdr.ests.c), ".u")] <- sapply(fdr.ests.c, function(x) x$ci[2])
      restab.add[, paste0(comp.type, ".", names(fdr.ests.c), ".x")] <- sapply(fdr.ests.c, function(x) x$x)
      restab.add[, paste0(comp.type, ".", names(fdr.ests.c), ".n")] <- sapply(fdr.ests.c, function(x) x$n)
    }
    restab <- rbind(restab, restab.add)
  }
  return(list(resl = resl, resimp = resimp, restab = restab))
}





















































