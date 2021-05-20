################################
## efficient influence functions - based on 2014 paper (same representation in binary and continuous exposures)
################################

.DcIPW <- function(n,
                   X,
                   Y,
                   Acol,
                   delta,
                   qfun,
                   gfun,
                   qfit,
                   gfits,
                   estimand,
                   wt,
                   isbin=FALSE,
                   ...){
  # define shifts
  Xa <- .shift(X,Acol, -delta)
  Xb <- .shift(X,Acol,  delta)
  #
  gn <- gfun(X,Acol,gfits=gfits)
  ga <- gfun(Xa,Acol,gfits=gfits)
  gb <- gfun(Xb,Acol,gfits=gfits)
  #
  #
  Haw = .Haw(gn, ga, gb)  # evaluated at A_i
  psi <- mean((Haw*(Y) - Y*(estimand != "mean"))*wt)
  dc1 <- Haw*(Y - 0)
  dc2 <- 0
  dc3 <- - psi - Y*(estimand != "mean")                # Y doesn't show up in Diaz,vdl 2012 b/c they are estimating mean Y|A+delta
  list(eif = as.vector(dc1 + dc2 + dc3)*wt, psi=psi)
}


.DbIPW <- function(n,
                   X,
                   Y,
                   Acol,
                   delta,
                   qfun,
                   gfun,
                   qfit,
                   gfits,
                   estimand,
                   wt,
                   isbin=FALSE,
                   ...){
  # define shifts
  #X0 <- vibr:::.shift(X,Acol, -X[,Acol])
  #X1 <- .shift(X,Acol,  (1-X[,Acol]))
  g0 <- 1-gfun(NULL,Acol,gfits=gfits)


  Haw <- .Hawb(g0, delta, X, Acol, retcols=1)
  psi <- mean(Haw*Y*wt) - mean(Y*(estimand != "mean")*wt)
  dc1 <- Haw*(Y - 0)
  dc2 <- 0
  dc3 <- - psi - Y*(estimand != "mean")                # Y doesn't show up in Diaz,vdl 2012 b/c they are estimating mean Y|A+delta
  list(eif = as.vector(dc1 + dc2 + dc3)*wt, psi=psi)
}


################################
#: estimating equations
################################

#.MakeiIpwEst <- .MakeTmleEst

.EstEqIPW <- function(n,
                      X,
                      Y,
                      delta,
                      qfun=NULL,
                      gfun,
                      qfit=NULL,
                      gfits,
                      estimand,
                      bounded=FALSE,
                      wt=rep(1,n),
                      isbin=FALSE,
                      ...){
  isbin_vec <- apply(X, 2, function(x) length(unique(x))==2)
  resmat <- matrix(NA, nrow=length(isbin_vec), ncol=3)
  for(Acol in seq_len(length(isbin_vec))){
    if(isbin_vec[Acol]){
      dphi <- .DbIPW(n,X,Y,Acol,delta,qfun=NULL,gfun,qfit=NULL,gfits,estimand,wt,isbin=isbin, ...)
    } else{
      dphi <- .DcIPW( n,X,Y,Acol,delta,qfun=NULL,gfun,qfit=NULL,gfits,estimand,wt,isbin=isbin, ...)
    }
    tm <- .MakeTmleEst(dphi)
    resmat[Acol,] <- tm
  }
  colnames(resmat) <- names(tm)
  rownames(resmat) <- names(X)
  resmat <- data.frame(resmat)
  resmat$p <- pnorm(-abs(resmat$z))*2
  resmat
}

################################
# expert wrappers
################################
.trained_ipw <- function(obj,
                          X,
                          Y,
                          delta,
                          qfun,
                          gfun,
                          estimand,
                          bounded,
                          updatetype){
  fittable <- .EstEqIPW(obj$n,X,Y,delta,qfun=NULL,gfun=.gfunction,qfit=NULL,gfits=obj$sl.gfits,estimand, bounded=FALSE,wt=obj$weights,isbin=obj$isbin)
  res <- list(
    res = fittable,
    qfit = NULL,
    gfits = obj$sl.gfits,
    binomial = obj$isbin,
    type = "IPW",
    weights=obj$weights
  )
  class(res) <- c("vibr.fit", class(res))
  res
}

#' @export
.varimp_ipw <- function(X,
                        Y,
                        V=NULL,
                        delta=0.1,
                        Y_learners=NULL,
                        Xdensity_learners=NULL,
                        Xbinary_learners=NULL,
                        verbose=TRUE,
                        estimand,
                        bounded=FALSE,
                        isbin=FALSE,
                        ...){
  obj = .prelims(X, Y, V, delta, Y_learners=NULL, Xbinary_learners, Xdensity_learners, verbose=verbose, isbin=isbin, ...)
  res = .trained_ipw(obj,X,Y,delta,qfun,gfun,estimand,bounded,updatetype)
  res
}


#' @importFrom future future value
#' @export
.varimp_ipw_boot <- function(X,
                             Y,
                             V=NULL,
                             delta=0.1,
                             Y_learners=NULL,
                             Xdensity_learners=NULL,
                             Xbinary_learners=NULL,
                             verbose=TRUE,
                             estimand="diff",
                             bounded=FALSE,
                             isbin=FALSE,
                             B=100,
                             showProgress=TRUE,
                             ...){
  if(is.null(isbin)) isbin <- as.logical((length(unique(Y))==2))
  est <- .varimp_ipw(X,Y,V,delta,Y_learners,Xdensity_learners,Xbinary_learners,verbose,estimand,bounded, isbin=isbin,...)
  rn <- rownames(est$res)
  n = length(Y)
  ee <- new.env()
  for(b in 1:B){
    #ridx <- sample(seq_len(n), n, replace=TRUE)
    ridx <- .bootsample(n)
    ee[[paste0("iter",b)]] <- future::future( {
      if(showProgress) cat(".")
      Xi = X[ridx,,drop=FALSE]
      Yi = Y[ridx]
      Vi = V[ridx,,drop=FALSE]
      obj = .prelims(X=Xi, Y=Yi, V=Vi, delta, Y_learners, Xbinary_learners, Xdensity_learners, verbose=verbose, isbin=isbin, ...)
      fittable <- .EstEqIPW(obj$n,Xi,Yi,delta,qfun=NULL,gfun=.gfunction,qfit=NULL,gfits=obj$sl.gfits, estimand,bounded,wt=obj$weights, isbin=obj$isbin)
      fittable$est
    }, seed=TRUE, lazy=TRUE)
  }
  bootests = do.call(rbind, as.list(future::value(ee)))
  if(showProgress) cat("\n")
  colnames(bootests) <- rn
  if(verbose) cat("\n")
  res <- list(
    est = est,
    boots = bootests,
    binomial = isbin,
    type = "IPW"
  )
  class(res) <- c("vibr.bootfit", class(res))
  res
}
