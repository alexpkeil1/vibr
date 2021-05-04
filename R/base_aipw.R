# estimators obtained by using updated TMLE algorithm (Diaz and vdl 2018, and setting all Q=0)
################################
## efficient influence functions
################################

.DcAIPW <- function(n,
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
  as.vector(dc1 + dc2 + dc3)
}

.DbAIPW <- function(n,
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
  Xbb <- .shift(X,Acol,2*delta)
  #
  gn <- gfun(X,Acol,gfits=gfits)
  gb <- gfun(Xb,Acol,gfits=gfits)
  #
  qinit = qfun(X, Acol,qfit=qfit)
  qbinit = qfun(Xb, Acol,qfit=qfit)
  #
  Haw = .Hawb(gn, delta, X, Acol)

  eqfb <- 0 # cancels out
  dc1 <- Haw*(Y - qinit)
  dc2 <- qbinit - eqfb
  dc3 <- eqfb - Y*(estimand != "mean")                # Y doesn't show up in Diaz,vdl 2012 b/c they are estimating mean Y|A+delta
  as.vector(dc1 + dc2 + dc3)
}



################################
#: estimating equations
################################

.MakeiAipwEst <- function(dphi, n){
  #summary(fit <- lm(dphi~1))
  est <- mean(dphi)
  D <- dphi-est
  avar = mean(D^2) # asymptotic variance
  se = sqrt(avar)/sqrt(n)
  c(est=est, se = se, z=est/se)
}

.EstEqAIPW <- function(n,
                       X,
                       Y,
                       delta,
                       qfun,
                       gfun,
                       qfit,
                       gfits,
                       estimand,
                       bounded=FALSE, # future use
                       ...){
  isbin_vec <- apply(X, 2, function(x) length(unique(x))==2)
  resmat <- matrix(NA, nrow=length(isbin_vec), ncol=3)
  for(Acol in seq_len(length(isbin_vec))){
    if(isbin_vec[Acol]){
      dphi <- .DbAIPW(n,X,Y,Acol,delta,qfun,gfun,qfit,gfits,estimand,...)
    } else{
      dphi <- .DcAIPW( n,X,Y,Acol,delta,qfun,gfun,qfit,gfits,estimand,...)
    }
    tm <- .MakeiAipwEst(dphi, n)
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
.varimp_aipw <- function(X,
                         Y,
                         delta=0.1,
                         Y_learners=NULL,
                         Xdensity_learners=NULL,
                         Xbinary_learners=NULL,
                         verbose=TRUE,
                         estimand,
                         bounded=FALSE,
                         ...){
  tasklist = .create_tasks(X,Y,delta)
  n = length(Y)
  if(verbose) cat(paste0("Default delta = ", delta, "\n"))

  isbin <- as.character((length(unique(Y))==2))
  sl.qfit <- .train_Y(X, Y, Y_learners, verbose=verbose, isbin)
  sl.gfits <- .train_allX(X, tasklist$slX, Xbinary_learners, Xdensity_learners, verbose=verbose)
  fittable <- .EstEqAIPW(n,X,Y,delta,qfun=.qfunction,gfun=.gfunction,qfit=sl.qfit,gfits=sl.gfits, estimand, bounded=FALSE, ...)
  res <- list(
    res = fittable,
    qfit = sl.qfit,
    gfits = sl.gfits,
    binomial = isbin,
    type = "AIPW"
  )
  class(res) <- c("vibr.fit", class(res))
  res
}

#' @export
.varimp_aipw_boot <- function(X,
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
  est <- .varimp_aipw(X,Y,delta,Y_learners,Xdensity_learners,Xbinary_learners,verbose,estimand,bounded,...)
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
    fittable <- .EstEqAIPW(n,Xi,Yi,delta,qfun=.qfunction,gfun=.gfunction,qfit=sl.qfit,gfits=sl.gfits, estimand,bounded, ...)
    bootests[b,] <- fittable$est
  }
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

