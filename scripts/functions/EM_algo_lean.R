
EM_algo_mixture_multi_Sig <- function(control, Y.em, S.em, MVphen_K, Sigl.em.init, R.em.init,
                            update.Sig = T, update.K = F, pi.init = NULL,
                            fac.model = c("fa", "pca")[1], bic.pen.mult = 0.5, wish.pri = F,
                            verbose = F){
  M <- length(control$omegaseq)
  R <- R.em.init
  N <- nrow(Y.em)
  ph.use <- colnames(Y.em)
  P <- length(ph.use)
  
  #################################
  # Initialise
  Sigl <- Sigl.em.init
  nSig <- length(Sigl)
  if(length(update.Sig) == 1)
    update.Sig <- rep(update.Sig, nSig)
  if(length(update.Sig) > 1){
    if(length(update.Sig) != nSig)
      stop(paste0("You have specified length(update.Sig) = ", length(update.Sig), ", but update.Sig should be logical vector of length either 1 or ",
                  length(Sigl.em.init), " (= length(Sigl.em.init))"))
  }
  if(length(update.K) == 1)
    update.K <- rep(update.K, nSig)
  if(length(update.K) > 1){
    if(length(update.K) != nSig)
      stop(paste0("You have specified length(update.K) = ", length(update.K), ", but update.K should be logical vector of length either 1 or ",
                  length(Sigl.em.init), " (= length(Sigl.em.init))"))
  }

  if(is.null(pi.init)){
    pimat <- matrix(1, M, nSig)
    pimat <- pimat / sum(pimat)
  } else {
    pimat <- pi.init
    # TODO: input check that pimat is correctly specified here
  }
  Ksig <- rep(MVphen_K, nSig)

  ##############################################
  #Calc initial likelihood matrix
  llikmat.out <- em.update.function(Y.em = Y.em, S.em = S.em, Sigl = Sigl, R = R, omegaseq = control$omegaseq, 
                                  pimat = pimat, meth = "just.obj", prior.in.obj = F, update.Sig = update.Sig, 
                                  Ksig = Ksig, update.K = update.K, fac.model = fac.model, bic.pen.mult = bic.pen.mult,
                                  recalc.llmat = F, wish.pri = wish.pri)
  llmat <- llikmat.out$llmat
  
  ###################################
  # Start EM algorithm
  objv <- c()
  itnum <- 1
  converged <- F
  while(!converged){
    
    #############################################################
    # Compute updated Sigl, rmat
    em.up.out <- em.update.function(Y.em = Y.em, S.em = S.em, Sigl = Sigl, R = R, omegaseq = control$omegaseq,
                               pimat = pimat, meth = "update.Sigl", update.Sig = update.Sig, Ksig = Ksig, update.K = update.K, 
                               fac.model = fac.model, bic.pen.mult = bic.pen.mult, llmat = llmat, wish.pri = wish.pri)
    Sigl <- em.up.out$Sigl
    Ksig <- em.up.out$Ksig
    rmat <- em.up.out$rmat
    llmat <- em.up.out$llmat

    ############################################
    # Optimize wrt pi (note: implemented for a uniform Dirichlet here)
    diralpha.new <- apply(rmat, 2:3, function(v) sum(v, na.rm = T))
    pimat[] <- diralpha.new / sum(diralpha.new)
    
    #############################
    # Evaluate objective function
    calc.obj.out <- em.update.function(Y.em = Y.em, S.em = S.em, Sigl = Sigl, R = R, omegaseq = control$omegaseq,
                                       pimat = pimat, Ksig = Ksig, meth = "just.obj", update.Sig = update.Sig, 
                                       fac.model = fac.model, bic.pen.mult = bic.pen.mult, llmat = llmat, recalc.llmat = F,
                                       wish.pri = wish.pri)
    calc.obj.out$obj
    objv[itnum] <- calc.obj.out$obj
    llikv <- calc.obj.out$llikv
    
    ##################################
    # Check convergence and print trace
    if(length(objv) > 1){
      most.recent.change <- objv[itnum] - objv[itnum - 1]
      llik.scale <- mad(llikv, na.rm = T)
      tolerance.eps <- llik.scale * N * control$MVphen_conv_tolerance
      if(most.recent.change < 0 & most.recent.change > -tolerance.eps){
        converged <- TRUE
      }
      dob <- objv[length(objv)] - objv[length(objv) - 1]
      if (verbose) {
        print(paste("Iteration = ", itnum))
        print(paste0("Ksig = ", Ksig))
        print(paste0("Obj = ", objv[length(objv)]))
        print(paste0("Change in obj = ", dob))
        if (ncol(pimat) > 1) {
          print("Pi = ")
          print(round(colSums(pimat) * 100, 2))
        }
        print(paste0("llik.scale = ", llik.scale))
        print(paste0("tolerance.eps = ", tolerance.eps))
        print(paste0("most.recent.change = ", most.recent.change))
        plot(objv, ty = "l")
      }
    }
    itnum <- itnum + 1
  }
  Sig.mn <- 0
  for(sc in 1:nSig){
    for(m in 1:M)
      Sig.mn <- Sig.mn + pimat[m, sc] * Sigl[[sc]] * control$omegaseq[m]
  }
  Sigl <- lapply(Sigl, function(M){ dimnames(M) <- list(ph.use, ph.use); M})
  dimnames(Sig.mn) <- list(ph.use, ph.use)
  return(list(Sigl = Sigl, Sig.mn = Sig.mn, Ksig = Ksig, R = R, pi = pimat, omegaseq = control$omegaseq, objv = objv))
}






























































