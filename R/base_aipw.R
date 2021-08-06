# estimators obtained by using updated TMLE algorithm (Diaz and vdl 2018, and setting all Q=0)
################################
## efficient influence functions
################################

.DcAIPW <- function(
  n,
  X,
  Y,
  Acol,
  delta,
  qfun,
  gfun,
  qfit,
  gfits,
  estimand,
  bounded,
  wt,
  isbin=FALSE,
  ...
){
  # define shifts
  Xa <- .shift(X,Acol, -delta)
  Xb <- .shift(X,Acol,  delta)
  #
  gn <- gfun(X,Acol,gfits=gfits)
  ga <- gfun(Xa,Acol,gfits=gfits)
  gb <- gfun(Xb,Acol,gfits=gfits)
  #
  qinit = qfun(X, Acol,qfit=qfit)
  qbinit = qfun(Xb, Acol,qfit=qfit)

  #
  #ga = .enforce_min_dens(ga,eps=1e-8)
  Haw = .Haw(gn, ga, gb)  # evaluated at A_i

  #eqfb = predict(lm(y~., data.frame(y=qfb, X=X[,-Acol])))   # simple linear regression on W to get E_g[Q | W]
  eqfb <- 0 # cancels out
  dc1 <- Haw*(Y - qinit)
  dc2 <- qbinit - eqfb
  dc3 <- eqfb - Y*(estimand != "mean")                # Y doesn't show up in Diaz,vdl 2012 b/c they are estimating mean Y|A+delta
  as.vector(dc1 + dc2 + dc3)*wt
}

.DbAIPW <- function(
  n,
  X,
  Y,
  Acol,
  delta,
  qfun,
  gfun,
  qfit,
  gfits,
  estimand,
  bounded,
  wt,
  isbin=FALSE,
  ...){
  # define shifts
  X0 <- .shift(X,Acol, -X[,Acol])
  X1 <- .shift(X,Acol,  (1-X[,Acol]))
  g0 <- 1-gfun(NULL,Acol,gfits=gfits)

  #
  qinit = qfun(X, Acol,qfit=qfit)
  q1init = qfun(X1, Acol,qfit=qfit)
  q0init = qfun(X0, Acol,qfit=qfit)
  #
  #ga = .enforce_min_dens(ga,eps=1e-8)

  Hawmat <- .Hawb(g0, delta, X, Acol, retcols=3)
  Haw  <- Hawmat[,1]

  eqfb <- 0 # cancels out
  dc1 <- Haw*(Y - qinit)
  dc2 <- qinit - eqfb
  dc3 <- delta*(q1init - q0init) + eqfb - Y*(estimand != "mean")                # Y doesn't show up in Diaz,vdl 2012 b/c they are estimating mean Y|A+delta
  as.vector(dc1 + dc2 + dc3)*wt
}



################################
#: estimating equations
################################

.MakeiAipwEst <- function(dphi){
  #summary(fit <- lm(dphi~1))
  est <- mean(dphi)
  D <- dphi-est
  avar = mean(D^2) # asymptotic variance
  se = sqrt(avar)/sqrt(length(D))
  c(est=est, se = se, z=est/se)
}

.EstEqAIPW <- function(n,
                       X,
                       Y,
                       whichcols=seq_len(ncol(X)),
                       delta,
                       qfun,
                       gfun,
                       qfit,
                       gfits,
                       estimand,
                       bounded=FALSE, # future use
                       wt=rep(1,n),
                       isbin=FALSE,
                       ...){
  if(length(whichcols>1)) {
    isbin_vec <- apply(X[,whichcols, drop=FALSE], 2, function(x) length(unique(x))==2)
  } else isbin_vec = length(unique(X[,whichcols]))==2

  resmat <- matrix(NA, nrow=length(isbin_vec), ncol=3)
  for(Acol in seq_len(length(isbin_vec))){
    if(isbin_vec[Acol]){
      #dphi <- .DbAIPW(n=n,X=X,Y=Y,Acol=Acol,delta=delta,qfun=qfun,gfun=gfun,qfit=qfit,gfits=gfits,estimand=estimand, bounded=bounded, wt=wt,isbin=isbin)
      aipwfun = .DbAIPW
    } else{
      #dphi <- .DcAIPW(n=n,X=X,Y=Y,Acol=Acol,delta=delta,qfun=qfun,gfun=gfun,qfit=qfit,gfits=gfits,estimand=estimand, bounded=bounded, wt=wt,isbin=isbin)
      aipwfun = .DcAIPW
    }
    dphi = aipwfun(n=n,X=X,Y=Y,Acol=Acol,delta=delta,qfun=qfun,gfun=gfun,qfit=qfit,gfits=gfits,estimand=estimand, bounded=bounded, wt=wt,isbin=isbin)
    tm <- .MakeiAipwEst(dphi)
    resmat[Acol,] <- tm
  }
  colnames(resmat) <- names(tm)
  rownames(resmat) <- names(X[,whichcols,drop=FALSE])
  resmat <- data.frame(resmat)
  resmat$p <- stats::pnorm(-abs(resmat$z))*2
  resmat
}


################################
# expert wrappers
################################


.trained_aipw <- function(
  obj,
  X,
  Y,
#  whichcols = seq_len(ncol(X)),
  delta,
  qfun,
  gfun,
  estimand,
  bounded,
  updatetype # place holder here
){
  fittable <- .EstEqAIPW(n=obj$n,X=X,Y=Y, whichcols=obj$whichcols,delta,qfun=.qfunction,gfun=.gfunction,qfit=obj$sl.qfit,gfits=obj$sl.gfits, estimand, bounded=FALSE,wt=obj$weights,isbin=obj$isbin)
  res <- list(
    res = fittable,
    qfit = obj$sl.qfit,
    gfits = obj$sl.gfits,
    binomial = obj$isbin,
    type = "AIPW",
    weights=obj$weights
  )
  class(res) <- c("vibr.fit", class(res))
  res
}

#' @export
.varimp_aipw <- function(
  X,
  Y,
  V=NULL,
  whichcols=seq_len(ncol(X)),
  delta=0.1,
  Y_learners=NULL,
  Xdensity_learners=NULL,
  Xbinary_learners=NULL,
  verbose=TRUE,
  estimand,
  bounded=FALSE,
  isbin=FALSE,
  ...){
  obj = .prelims(X=X, Y=Y, V=V, whichcols=whichcols, delta, Y_learners, Xbinary_learners, Xdensity_learners, verbose=verbose,isbin=isbin, ...)
  #res = .trained_aipw(obj,X,Y,delta,qfun,gfun,estimand,bounded,updatetype)
  res = .trained_aipw(obj=obj,X=X,Y=Y,delta=delta,qfun=.qfunction,gfun=.gfunction,estimand=estimand,bounded=bounded,updatetype=NULL)
  res
}

#' @importFrom future future value
#' @export
.varimp_aipw_boot <- function(X,
                              Y,
                              V=NULL,
                              whichcols=seq_len(ncol(X)),
                              delta=0.1,
                              Y_learners=NULL,
                              Xdensity_learners=NULL,
                              Xbinary_learners=NULL,
                              verbose=TRUE,
                              estimand="diff",
                              isbin=NULL,
                              bounded=FALSE,
                              B=100,
                              showProgress=TRUE,
                              ...){
  if(is.null(isbin)) isbin <- as.logical((length(unique(Y))==2))
  est <- .varimp_aipw(X=X,Y=Y,V=V,whichcols=whichcols,delta,Y_learners,Xdensity_learners,Xbinary_learners,verbose,estimand,bounded,isbin=isbin,...)
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
      obj = .prelims(X=Xi, Y=Yi, V=Vi, whichcols=whichcols, delta, Y_learners, Xbinary_learners, Xdensity_learners, verbose=verbose,isbin=isbin, ...)
      fittable <- .EstEqAIPW(n=obj$n,X=Xi,Y=Yi, whichcols=obj$whichcols,delta=delta,qfun=.qfunction,gfun=.gfunction,qfit=obj$sl.qfit,gfits=obj$sl.gfits, estimand=estimand,bounded=bounded,wt=obj$weights,isbin=obj$isbin)
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
    type = "AIPW"
  )
  class(res) <- c("vibr.bootfit", class(res))
  res
}

