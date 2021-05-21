# todo:
#  allow "whichcols" variables

################################################################################
#
#    g-computation for joint effects
#
################################################################################
.gcjoint <- function(
  n,
  X,
  Y,
  expnms,
  delta,
  qfun,
  gfun=NULL,
  qfit,
  gfits=NULL,
  estimand,
  wt,
  ...
){
  #cat(paste0("column ", names(X)[Acol], ": binary\n"))
  if(length(ncol(X)>1)) {
    isbin_vec <- apply(X, 2, function(x) length(unique(x))==2)
  } else isbin_vec = length(unique(X[,1]))==2
  X1 <- X0 <- X1b <- X0b <- Xb <- X
  nmx = names(X)
  for(Acol in seq_len(length(isbin_vec))){
    if(nmx[Acol] %in% expnms){
      if(isbin_vec[Acol]){
        Xb[,Acol] <- Xb[,Acol]+delta
        X0b[,Acol] <- X0b[,Acol]+delta
        X1b[,Acol] <- X1b[,Acol]+delta
      } else{
        X0 <- .shift(X1, Acol, -X[,Acol])
        X1 <- .shift(X0, Acol,  (1-X[,Acol]))
        X0b <- .shift(X1b, Acol, -X[,Acol])
        X1b <- .shift(X0b, Acol,  (1-X[,Acol]))
      }
    }
  }
  # X1: Acs = A+delta, Ab=1
  # X0: Acs = A+delta, Ab=0
  #qfun(Xb,Acol,qfit=qfit,...)
  # E(Y | Abs,Acs) = E(Y|Ab,Ac) +
  #                  dE(Y|Ab,Ac)/dAc +
  #                  dE(Y|Ab,Ac)/dAb +
  #                  dE(Y|Ab,Acs)/dAb
  #
  # E(Y | Abs,Acs) = E(Y|Ab,Ac)  +
  #                  (E(Y|Ab,Acs)-E(Y|Ab,Ac)) +
  #                  (E(Y|Abs,Ac)-E(Y|Ab,Ac))
  #                  (E(Y|Abs,Acs)-E(Y|Ab,Acs)) - (E(Y|Abs,Ac)-E(Y|Ab,Ac))
  #
  #p1 <- qfun(qfit=qfit,...)           # E(Y|Ab,ac)
  p2 <- qfun(Xb,Acol,qfit=qfit,...)   # E(Y|Ab,ac) + (E(Y|Ab,Acs)-E(Y|Ab,Ac))
  p3 <- 0#delta*(qfun(X1,Acol,qfit=qfit,...) - qfun(X0,Acol,qfit=qfit,...)) # (E(Y|Abs,Ac)-E(Y|Ab,Ac))
  p4 <- delta*(qfun(X1b,Acol,qfit=qfit,...) - qfun(X0b,Acol,qfit=qfit,...)) # (E(Y|Abs,Acs)-E(Y|Ab,Acs))
  (p2 + p3 + p4 - Y*(estimand != "mean"))*wt
}


.EstimatorGcompJoint <- function(
  n,
  X,
  Y,
  expnms,
  delta,
  qfun,
  gfun,
  qfit,
  gfits,
  estimand,
  bounded=FALSE,
  wt=rep(1,n),
  ...
){
  #isbin_vec <- apply(X, 2, function(x) length(unique(x))==2)
  resmat <- matrix(NA, nrow=1, ncol=3)
  phi <- .gcjoint(n=n,X=X,Y=Y,expnms,delta=delta,qfun=qfun,gfun=NULL,qfit=qfit,gfits=NULL,estimand=estimand, wt=wt, ...)
  tm <- .EstGcomp(phi)
  resmat[1,] <- tm
  colnames(resmat) <- names(tm)
  rownames(resmat) <- "Joint"
  resmat <- data.frame(resmat)
  resmat$p <- pnorm(-abs(resmat$z))*2
  resmat
}


.trained_gcomp_joint <- function(
  obj,
  X,
  Y,
  expnms,
  delta,
  qfun,
  gfun,
  estimand,
  bounded,
  updatetype
){
  fittable <- .EstimatorGcompJoint(obj$n,X,Y,expnms,delta,qfun=.qfunction,gfun=NULL,qfit=obj$sl.qfit,gfits=NULL, estimand, bounded,wt=obj$weights)
  res <- list(
    res = fittable,
    qfit = obj$sl.qfit,
    gfits = NULL,
    binomial = obj$isbin,
    type = "GCOMP",
    weights=obj$weights
  )
  class(res) <- c("vibr.fit", class(res))
  res
}


.varimp_gcomp_joint <- function(
  X,
  Y,
  V=NULL,
  expnms = NULL,
  delta=0.1,
  Y_learners=NULL,
  Xdensity_learners=NULL,
  Xbinary_learners=NULL,
  verbose=TRUE,
  estimand,
  bounded=FALSE,
  ...
){
  scale_continuous = TRUE
  if(scale_continuous){
    if(verbose) cat("Scaling all continuous variables by 2*sd\n")
    # divide continuous by 2*sd
    X = .scale_continuous(X)
  }
  obj = .prelims(
    X=X,
    Y=Y,
    V=V,
    delta=delta,
    Y_learners=Y_learners,
    Xbinary_learners=NULL,
    Xdensity_learners=NULL,
    verbose=verbose,
    ...
  )
  res = .trained_gcomp_joint(
    obj=obj,
    X=X,
    Y=Y,
    expnms=expnms,
    delta=delta,
    qfun=.qfunction,
    gfun=.gfunction,
    estimand=estimand,
    bounded=bounded,
    updatetype=updatetype
  )
  res <- .attach_misc(res, scale_continuous=scale_continuous, delta=delta, B=NULL)
  class(res) <- c("vibr.fit", class(res))
  res
}



#'
.varimp_gcomp_joint_boot <- function(
  X,
  Y,
  V=NULL,
  expnms = NULL,
  delta=0.1,
  Y_learners=NULL,
  Xdensity_learners=NULL,
  Xbinary_learners=NULL,
  verbose=TRUE,
  estimand="diff",
  bounded=FALSE,
  B=100,
  showProgress=TRUE,
  ...
){
  scale_continuous = TRUE
  if(scale_continuous){
    if(verbose) cat("Scaling all continuous variables by 2*sd\n")
    # divide continuous by 2*sd
    X = .scale_continuous(X)
  }
  est <- .varimp_gcomp_joint(
    X=X,
    Y=Y,
    V=V,
    expnms=expnms,
    delta=delta,
    Y_learners=Y_learners,
    Xdensity_learners=NULL,
    Xbinary_learners=NULL,
    verbose=verbose,
    estimand=estimand,
    bounded=bounded,
    ...
  )
  rn <- rownames(est$res)
  n = length(Y)
  isbin <- as.character((length(unique(Y))==2))
  ee <- new.env()
  for(b in 1:B){
    ridx <- .bootsample(n)
    ee[[paste0("iter",b)]] <- future::future({
      if(showProgress) cat(".")
      Xi = X[ridx,,drop=FALSE]
      Yi = Y[ridx]
      Vi = V[ridx,,drop=FALSE]
      obj = .prelims(
        X=Xi,
        Y=Yi,
        V=Vi,
        delta=delta,
        Y_learners=Y_learners,
        Xbinary_learners=NULL,
        Xdensity_learners=NULL,
        verbose=verbose,
        ...
      )
      fittable <- .EstimatorGcompJoint(
        obj$n,
        X=Xi,
        Y=Yi,
        expnms=expnms,
        delta=delta,
        qfun=.qfunction,
        gfun=NULL,
        qfit=obj$sl.qfit,
        gfits=NULL,
        estimand=estimand,
        bounded=bounded,
        wt=obj$weights
      )
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
    type = "GCOMP"
  )
  res$est <- .attach_misc(res$est, scale_continuous=scale_continuous, delta=delta, B=NULL)
  res <- .attach_misc(res, scale_continuous=scale_continuous, delta=delta, B=B)

  class(res) <- c("vibr.bootfit", class(res))
  res

}
