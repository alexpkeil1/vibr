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
                 ...){
  #cat(paste0("column ", names(X)[Acol], ": continuous\n"))
  Xb <- X
  Xb[,Acol] <- X[,Acol]+delta
  p2 <- qfun(Xb,Acol,qfit=qfit,...)
  sum(p2 - Y*(estimand != "mean"))/n
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
                 ...){
  #cat(paste0("column ", names(X)[Acol], ": binary\n"))
  X1 <- X0 <- X
  X1[,Acol] <- 1
  X0[,Acol] <- 0
  p2 <- qfun(qfit=qfit,...)
  p3 <- delta*(qfun(X1,Acol,qfit=qfit,...) - qfun(X0,Acol,qfit=qfit,...))
  sum(p2 + p3 - Y*(estimand != "mean"))/n
}


.EstGcomp <- function(est){
  #summary(fit <- lm(dphi~1))
  D <- est
  se = sqrt(mean(D^2)/length(D))
  c(est=est, se = NA, z=NA)
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
                            ...){
  isbin_vec <- apply(X, 2, function(x) length(unique(x))==2)
  resmat <- matrix(NA, nrow=length(isbin_vec), ncol=3)
  for(Acol in seq_len(length(isbin_vec))){
    if(isbin_vec[Acol]){
      est <- .gcb(n,X,Y,Acol,delta,qfun,gfun=NULL,qfit,gfits=NULL,estimand, ...)
    } else{
      est <- .gcc( n,X,Y,Acol,delta,qfun,gfun=NULL,qfit,gfits=NULL,estimand, ...)
    }
    tm <- .EstGcomp(est)
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

# \dontrun{
# XYlist = .dgm(n=100,p=4,ncat=3)
# data(metals, package="qgcomp")
# XYlist = list(X=metals[,1:23], Y=metals$y)
# Y_learners = .default_continuous_learners()
# Xbinary_learners = .default_binary_learners()
# Xdensity_learners = .default_density_learners(n_bins=c(5, 20))
# vi2 <- .varimp_gcomp(X=XYlist$X,Y=XYlist$Y, delta=0.1, Y_learners = Y_learners[1:4],Xdensity_learners=Xdensity_learners[1:2], Xbinary_learners=Xbinary_learners[1:2] )
# vi2
# }
# coef(summary(lm(XYlist$Y~as.matrix(XYlist$X))))[-1,]
# y = rnorm(N, 0, 0.5) + (0.35 + ph*0.1)*calcium*sodium*zinc
# + (0.35 + ph*0.1)*iron*selenium*selenium
# + 1*(calcium>mean(calcium))
# - (0.1 + ph*0.1)*lead*cadmium*arsenic
# - (0.1 + ph*0.1)*chromium*mercury,

#' Variable importance using g-computation
#' @description Not usually called by users
#'
#' @param X data frame of predictors
#' @param Y outcome
#' @param delta change in each column of X corresponding to
#' @param Y_learners list of sl3 learners used to predict the outcome, conditional on all predictors in X
#' @param Xdensity_learners list of sl3 learners used to estimate the density of continuous predictors, conditional on all other predictors in X
#' @param Xbinary_learners list of sl3 learners used to estimate the probability mass of continuous predictors, conditional on all other predictors in X
#' @param verbose (logical) print extra information
#' @param ... passed to sl3::base_predict (rare)
#'
#' @return vi object
#' @export
#'
.varimp_gcomp <- function(X,
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
  if(verbose) cat(paste0("Default delta = ", delta, "\n")) # TODO: better interpretation

  isbin <- as.character((length(unique(Y))==2))
  sl.qfit <- .train_Y(X,Y, Y_learners, verbose=verbose, isbin)
  fittable <- .EstimatorGcomp(n,X,Y,delta,qfun=.qfunction,gfun=NULL,qfit=sl.qfit,gfits=NULL, estimand, bounded,...)
  res <- list(
    res = fittable,
    qfit = sl.qfit,
    gfits = NULL,
    binomial = isbin,
    type = "GCOMP"
  )
  class(res) <- c("vibr.fit", class(res))
  res
}


.varimp_gcomp_boot <- function(X,
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
  est <- .varimp_gcomp(X,Y,delta,Y_learners,Xdensity_learners,Xbinary_learners,verbose,estimand,bounded,...)
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
    fittable <- .EstimatorGcomp(n,Xi,Yi,delta,qfun=.qfunction,gfun=.gfunction,qfit=sl.qfit,gfits=sl.gfits, estimand,bounded, ...)
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
