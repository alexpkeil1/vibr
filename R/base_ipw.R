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
  dc1 <- Haw*(Y - 0)
  dc2 <- 0
  dc3 <- - Y*(estimand != "mean")                # Y doesn't show up in Diaz,vdl 2012 b/c they are estimating mean Y|A+delta
  as.vector(dc1 + dc2 + dc3)
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
                   ...){
  # define shifts
  Xa <- .shift(X,Acol, -delta)
  Xb <- .shift(X,Acol,  delta)
  #
  gn <- gfun(X,Acol,gfits=gfits)
  gb <- gfun(Xb,Acol,gfits=gfits)
  #
  Haw = .Hawb(gn, delta, X, Acol)

  dc1 <- Haw*(Y - 0)
  dc2 <- 0
  dc3 <- - Y*(estimand != "mean")                # Y doesn't show up in Diaz,vdl 2012 b/c they are estimating mean Y|A+delta
  as.vector(dc1 + dc2 + dc3)
}


################################
#: estimating equations
################################


.EstEqIPW <- function(n,X,Y,delta,qfun=NULL,gfun,qfit=NULL,gfits,estimand, bounded=FALSE,...){
  isbin_vec <- apply(X, 2, function(x) length(unique(x))==2)
  resmat <- matrix(NA, nrow=length(isbin_vec), ncol=3)
  for(Acol in seq_len(length(isbin_vec))){
    if(isbin_vec[Acol]){
      dphi <- .DbIPW(n,X,Y,Acol,delta,qfun=NULL,gfun,qfit=NULL,gfits,estimand, ...)
    } else{
      dphi <- .DcIPW( n,X,Y,Acol,delta,qfun=NULL,gfun,qfit=NULL,gfits,estimand, ...)
    }
    tm <- .MakeiAipwEst(dphi,n)
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


#' @export
.varimp_ipw <- function(X,Y, delta=0.1, Y_learners=NULL, Xdensity_learners=NULL, Xbinary_learners=NULL, verbose=TRUE,estimand, bounded=FALSE, ...){
  tasklist = .create_tasks(X,Y,delta)
  n = length(Y)
  if(verbose) cat(paste0("Default delta = ", delta, "\n")) # TODO: better interpretation

  isbin <- as.character((length(unique(Y))==2))
  #sl.qfit <- .train_Y(X,Y, Y_learners, verbose=TRUE, isbin)
  sl.gfits <- .train_allX(X, tasklist$slX, Xbinary_learners, Xdensity_learners, verbose=verbose)
  #.gfunction(X=NULL,Acol,gfits=sl.gfits)
  #.qfunction(X=NULL,Acol,qfit=sl.qfit)
  fittable <- .EstEqIPW(n,X,Y,delta,qfun=NULL,gfun=.gfunction,qfit=NULL,gfits=sl.gfits,estimand, bounded=FALSE, ...)
  res <- list(
    res = fittable,
    qfit = NULL,
    gfits = sl.gfits,
    binomial = isbin,
    type = "IPW"
  )
  class(res) <- c("vibr.fit", class(res))
  res
}


#' @export
.varimp_ipw_boot <- function(X,
                              Y,
                              delta=0.1,
                              Y_learners=NULL,
                              Xdensity_learners=NULL,
                              Xbinary_learners=NULL,
                              verbose=TRUE,
                              estimand="diff",
                              bounded=FALSE,
                              B=100,
                             showProgress=TRUE,
                             ...){
  est <- .varimp_ipw(X,Y,delta,Y_learners,Xdensity_learners,Xbinary_learners,verbose,estimand,bounded,...)
  rn <- rownames(est$res)
  bootests <- matrix(NA, nrow=B, ncol = length(rn))
  n = length(Y)
  isbin <- as.character((length(unique(Y))==2))
  for(b in 1:B){
    if(showProgress) cat(".") # TODO: better interpretation
    ridx <- sample(seq_len(n), n, replace=TRUE)
    Xi = X[ridx,,drop=FALSE]
    Yi = Y[ridx]
    tasklist = .create_tasks(Xi,Yi,delta)
    yb = .bound_zero_one(Yi)
    Ybound = yb[[1]]
    sl.qfit <- .train_Y(Xi,Yi, Y_learners, verbose=FALSE, isbin)
    sl.gfits <- .train_allX(Xi, tasklist$slX, Xbinary_learners, Xdensity_learners, verbose=verbose)
    fittable <- .EstEqIPW(n,Xi,Yi,delta,qfun=.qfunction,gfun=.gfunction,qfit=sl.qfit,gfits=sl.gfits, estimand,bounded, ...)
    bootests[b,] <- fittable$est
  }
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
