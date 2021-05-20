################################
#: estimator
################################

# continuous
.gcc <- function(n,
                 X,
                 Y,
                 Acol,
                 delta,
                 qfun,
                 gfun=NULL,
                 qfit,
                 gfits=NULL,
                 estimand,
                 wt,
                 isbin=FALSE,
                 ...){
  #cat(paste0("column ", names(X)[Acol], ": continuous\n"))
  Xb <- X
  Xb[,Acol] <- X[,Acol]+delta
  p2 <- qfun(Xb,Acol,qfit=qfit,...)
  (p2 - Y*(estimand != "mean"))*wt
}

# binary
.gcb <- function(n,
                 X,
                 Y,
                 Acol,
                 delta,
                 qfun,
                 gfun=NULL,
                 qfit,
                 gfits=NULL,
                 estimand,
                 wt,
                 isbin=FALSE,
                 ...){
  #cat(paste0("column ", names(X)[Acol], ": binary\n"))
  X0 <- .shift(X,Acol, -X[,Acol])
  X1 <- .shift(X,Acol,  (1-X[,Acol]))
  p2 <- qfun(qfit=qfit,...)
  p3 <- delta*(qfun(X1,Acol,qfit=qfit,...) - qfun(X0,Acol,qfit=qfit,...))
  (p2 + p3 - Y*(estimand != "mean"))*wt
}



.EstGcomp <- function(phi){
  #summary(fit <- lm(dphi~1))
  est = mean(phi)
  D <- phi-est
  se = sqrt(mean(D^2)/length(D))
  c(est=est, se = se, z=est/se)
}

.EstimatorGcomp <- function(n,
                            X,
                            Y,
                            delta,
                            qfun,
                            gfun,
                            qfit,
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
      phi <- .gcb(n=n,X=X,Y=Y,Acol=Acol,delta=delta,qfun=qfun,gfun=NULL,qfit=qfit,gfits=NULL,estimand=estimand, wt=wt, isbin=isbin, ...)
    } else{
      phi <- .gcc(n=n,X=X,Y=Y,Acol=Acol,delta=delta,qfun=qfun,gfun=NULL,qfit=qfit,gfits=NULL,estimand=estimand, wt=wt, isbin=isbin,...)
    }
    tm <- .EstGcomp(phi)
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
.trained_gcomp <- function(obj,
                         X,
                         Y,
                         delta,
                         qfun,
                         gfun,
                         estimand,
                         bounded,
                         updatetype){
  fittable <- .EstimatorGcomp(obj$n,X,Y,delta,qfun=.qfunction,gfun=NULL,qfit=obj$sl.qfit,gfits=NULL, estimand, bounded,wt=obj$weights, isbin=obj$isbin)
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


#' @export
.varimp_gcomp <- function(X,
                          Y,
                          V=NULL,
                          delta=0.1,
                          Y_learners=NULL,
                          Xdensity_learners=NULL,
                          Xbinary_learners=NULL,
                          verbose=TRUE,
                          estimand,
                          bounded=FALSE,
                          isbin=NULL,
                          ...){
  obj = .prelims(X, Y, V, delta, Y_learners, Xbinary_learners=NULL, Xdensity_learners=NULL, verbose=verbose, isbin=isbin, ...)
  res = .trained_gcomp(obj,X,Y,delta,qfun,gfun,estimand,bounded,updatetype)
  res
}


#' @importFrom future future value
#' @export
.varimp_gcomp_boot <- function(X,
                               Y,
                               V=NULL,
                               delta=0.1,
                               Y_learners=NULL,
                               Xdensity_learners=NULL,
                               Xbinary_learners=NULL,
                               verbose=TRUE,
                               estimand="diff",
                               bounded=FALSE,
                               isbin=NULL,
                               B=100,
                               showProgress=TRUE,
                               ...){
  if(is.null(isbin)) isbin <- as.logical((length(unique(Y))==2))
  est <- .varimp_gcomp(X,Y,V,delta,Y_learners,Xdensity_learners,Xbinary_learners,verbose,estimand,bounded,isbin=isbin,...)
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
      obj = .prelims(X=Xi, Y=Yi, V=Vi, delta, Y_learners, Xbinary_learners, Xdensity_learners, verbose=verbose,isbin=isbin, ...)
      fittable <- .EstimatorGcomp(obj$n,Xi,Yi,delta,qfun=.qfunction,gfun=.gfunction,qfit=obj$sl.qfit,gfits=obj$sl.gfits, estimand,bounded,wt=obj$weights, isbin=obj$isbin)
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
  class(res) <- c("vibr.bootfit", class(res))
  res
}
