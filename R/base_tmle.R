# TMLE using approach from Targeted Learning in Data Science (S 14.3), Diaz and van der Laan (2018)
# note this is a much easier method than that given in Diaz and van der Laan (2012)
################################
## efficient influence functions
################################
.enforce_min_dens <- function(x,eps=1e-8){
  ifelse(x<eps, eps, x)
}

# note d(a|w) = A+delta, and d^-1(a|w) = h(a|w) = A-delta

# clever covariate, general or continuous
.Haw <- function(gn, ga, gb){
  # I(A<u(w))g(a-delta|w)/g(a|w) + I(a>=u(w)-delta)
  Haw <- ifelse(gn>0, ga/gn, 0) + as.numeric(gb == 0)
  Haw
}
# clever covariate, binary (based on diaz and vdl 2012)
.Hawb_old <- function(gn, ga, gb, shift, X, Acol){
  # if intervention would push exposure out of the support of A | W, then don't intervene
  # (delta*(I(A=1)-I(A=0)) + gn)/gn + gn/gn
  #
  Haw0 <- ifelse(gn>0, (-shift/gn + 1), 0)+ as.numeric(gb == 0)
  Haw1 <- ifelse(gn>0, ( shift/gn + 1), 0)+ as.numeric(gb == 0)
  Haw <- ifelse(gn>0, (shift*(2*X[,Acol] - 1)/gn + 1), 0)+ as.numeric(gb == 0)
  Haw
}

# clever covariate, binary (based on diaz and vdl 2018, translated to shift in propensity score)
.Hawb <- function(gn, shift, X, Acol){
  # if intervention would push exposure out of the support of A | W, then don't intervene
  Haw1 <- ifelse(gn>0, -shift/gn + 1, 0) + as.numeric(gn + shift > 1)
  Haw0 <- ifelse((1-gn)>0, shift/(1-gn) + 1, 0) + as.numeric((1-gn) + shift > 1)
  Haw <- X[,Acol]*Haw1 + (1-X[,Acol])*Haw0
  Haw
}


.OneStepTmleCont <- function(Y,Qinit,Qdawinit,Haw,Hdaw,isbin=FALSE, weighted=FALSE){
  # Q_update(A,W) = Q(A,W) + eps*H(A,W)
  # 1. estimand epsilon
  if(weighted){
    epsk <- as.numeric(lm(Y~ offset(Qinit), weights = Haw)$coefficients[1])
  } else{
    epsk <- as.numeric(lm(Y~ -1 + offset(Qinit) + Haw)$coefficients[1])
  }
  # 2. update Qk
  Qk1 <- Qinit + epsk*Haw
  Qdawk1 <- Qdawinit + epsk*Hdaw
  cbind(Qk1, Qdawk1)
}

.OneStepTmleBin <- function(Y,Qinit,Qdawinit,Haw,Hdaw,isbin=FALSE){
  # Q_update(A,W) = expit(logit(Q(A,W)) + eps*H(A,W))
  # 1. estimand epsilon
  epsk <- as.numeric(glm(Y~ -1 + offset(.logit(Qinit)) + Haw, family=binomial(link="logit"))$coefficients[1])
  # 2. update Qk
  Qk1 <- .expit(.logit(Qinit) + epsk*Haw)
  Qdawk1 <- .expit(.logit(Qdawinit) + epsk*Hdaw)
  cbind(Qk1, Qdawk1)
}

.bound_zero_one <- function(Y){
  miny = min(Y)
  maxy = max(Y)
  rngy = maxy-miny
  Ybound <- (Y - miny)/rngy
  list(Ybound, miny, rngy)
}


.DcTMLE <- function(n,
                X,
                Y,
                Acol,
                delta,
                qfun,
                gfun,
                qfit,
                gfits,
                est = "mean",
                bounded,
                ...
                ){
  # define shifts
  Xa <- .shift(X,Acol, -delta)
  Xb <- .shift(X,Acol,  delta)
  Xbb <- .shift(X,Acol,2*delta)
  #
  gn <- gfun(X,Acol,gfits=gfits)
  ga <- gfun(Xa,Acol,gfits=gfits)
  gb <- gfun(Xb,Acol,gfits=gfits)
  gbb <- gfun(Xbb,Acol,gfits=gfits)
  #
  qinit = qfun(X, Acol,qfit=qfit)
  qbinit = qfun(Xb, Acol,qfit=qfit)
  #
  #ga = .enforce_min_dens(ga,eps=1e-8)
  Haw = .Haw(gn, ga, gb)  # evaluated at A_i
  Hdaw = .Haw(gb, gn, gbb) # evaluated at d(A_i,W_i)

  QuMat <- .OneStepTmleCont(Y,qinit,qbinit,Haw,Hdaw,isbin=FALSE)
  Qupdate <- QuMat[,1,drop=TRUE]
  Qdawupdate <- QuMat[,2,drop=TRUE]
  #.OneStepTmleBin(Y,qinit,Haw,Hdaw,isbin=FALSE)
  #eqfb = predict(lm(y~., data.frame(y=qfb, X=X[,-Acol])))   # simple linear regression on W to get E_g[Q | W]
  eqfb <- 0 # cancels out
  dc1 <- Haw*(Y - Qupdate)
  dc2 <- Qdawupdate - eqfb
  dc3 <- eqfb - Y*(estimand != "mean")                # Y doesn't show up in Diaz,vdl 2012 b/c they are estimating mean Y|A+delta
  as.vector(dc1 + dc2 + dc3)
}


.DbTMLE <- function(n,
                               X,
                               Y,
                               Acol,
                               delta,
                               qfun,
                               gfun,
                               qfit,
                               gfits,
                               est = "mean",
                               bounded,
                               ...
){
  # define shifts
  Xa <- .shift(X,Acol, -delta)
  Xb <- .shift(X,Acol,  delta)
  Xbb <- .shift(X,Acol,2*delta)
  #
  gn <- gfun(X,Acol,gfits=gfits)
  #ga <- gfun(Xa,Acol,gfits=gfits)
  gb <- gfun(Xb,Acol,gfits=gfits)
  #gbb <- gfun(Xbb,Acol,gfits=gfits)
  #
  qinit = qfun(X, Acol,qfit=qfit)
  qbinit = qfun(Xb, Acol,qfit=qfit)
  #
  #ga = .enforce_min_dens(ga,eps=1e-8)

  Haw = .Hawb(gn, delta, X, Acol)
  Hdaw = .Hawb(gb, 2*delta, X, Acol)

  QuMat <- .OneStepTmleCont(Y,qinit,qbinit,Haw,Hdaw,isbin=FALSE)
  Qupdate <- QuMat[,1,drop=TRUE]
  Qdawupdate <- QuMat[,2,drop=TRUE]
  #.OneStepTmleBin(Y,qinit,Haw,Hdaw,isbin=FALSE)
  #eqfb = predict(lm(y~., data.frame(y=qfb, X=X[,-Acol])))   # simple linear regression on W to get E_g[Q | W]
  eqfb <- 0 # cancels out
  dc1 <- Haw*(Y - Qupdate)
  dc2 <- Qdawupdate - eqfb
  dc3 <- eqfb - Y*(estimand != "mean")                # Y doesn't show up in Diaz,vdl 2012 b/c they are estimating mean Y|A+delta
  as.vector(dc1 + dc2 + dc3)
}

################################
#: estimating equations
################################

.MakeTmleEst <- function(dphi){
  #summary(fit <- lm(dphi~1))
  est <- mean(dphi)
  D <- dphi-est
  se = sqrt(mean(D^2)/length(D))
  c(est=est, se = se, z=est/se)
}

.EstEqTMLE <- function(n,
                       X,
                       Y,
                       delta,
                       qfun,
                       gfun,
                       qfit,
                       gfits,
                       est,
                       bounded,
                       ...){
  #return(NULL)# remove when done
  isbin_vec <- apply(X, 2, function(x) length(unique(x))==2)
  resmat <- matrix(NA, nrow=length(isbin_vec), ncol=3)
  for(Acol in seq_len(length(isbin_vec))){
    if(isbin_vec[Acol]){
      dphi <- .DbTMLE(n,X,Y,Acol,delta,qfun,gfun,qfit,gfits,est,bounded)
      #dphi <- .DbTMLE(n,X,Y,Acol,delta,qfun,gfun,qfit,gfits,...)
    } else{
      dphi <- .DcTMLE( n,X,Y,Acol,delta,qfun,gfun,qfit,gfits,est,bounded)
      #dphi <- .DcTMLE( n,X,Y,Acol,delta,qfun,gfun,qfit,gfits,...)
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
#' Variable importance using targeted minimum loss estimation
#' @description Not usually called by users
#'
#' @param X data frame of predictors
#' @param Y outcome
#' @param delta change in each column of X corresponding to
#' @param Y_learners list of sl3 learners used to predict the outcome, conditional on all predictors in X
#' @param Xdensity_learners list of sl3 learners used to estimate the density of continuous predictors, conditional on all other predictors in X
#' @param Xbinary_learners list of sl3 learners used to estimate the probability mass of continuous predictors, conditional on all other predictors in X
#' @param verbose (logical) print extra information
#' @param estimand (character) "diff" (default, estimate mean difference comparing Y under intervention with observed Y), "mean" (estimate mean Y under intervention)
#' @param bounded (logical) does nothing yet
#' @param ... passed to sl3::base_predict (rare)
#'
#' @return vi object
#' @export
#'
#' @examples
#' \dontrun{
#' XYlist = .dgm(n=100,p=4,ncat=3)
#' data(metals, package="qgcomp")
#' XYlist = list(X=metals[,1:23], Y=metals$y)
#' Y_learners = .default_continuous_learners()
#' Xbinary_learners = .default_binary_learners()
#' Xdensity_learners = .default_density_learners(n_bins=c(5, 20))
#' vi <- .varimp_aipw(X=XYlist$X,Y=XYlist$Y, delta=0.1, Y_learners = Y_learners[1:4],
#' Xdensity_learners=Xdensity_learners[1:2], Xbinary_learners=Xbinary_learners[1:2] )
#' vi
#' }
.varimp_tmle <- function(X,
                         Y,
                         delta=0.1,
                         Y_learners=NULL,
                         Xdensity_learners=NULL,
                         Xbinary_learners=NULL,
                         verbose=TRUE,
                         estimand="diff",
                         bounded=FALSE,
                         ...){
  tasklist = .create_tasks(X,Y,delta)
  n = length(Y)
  if(verbose) cat(paste0("delta = ", delta, "\n")) # TODO: better interpretation

  isbin <- as.character((length(unique(Y))==2))
  yb = .bound_zero_one(Y)
  Ybound = yb[[1]]
  sl.qfit <- .train_Y(X,Y, Y_learners, verbose=verbose, isbin)
  sl.gfits <- .train_allX(X, tasklist$slX, Xbinary_learners, Xdensity_learners, verbose=verbose)
  fittable <- .EstEqTMLE(n,X,Y,delta,qfun=.qfunction,gfun=.gfunction,qfit=sl.qfit,gfits=sl.gfits, estimand,bounded, ...)
  res <- list(
    res = fittable,
    qfit = sl.qfit,
    gfits = sl.gfits,
    binomial = isbin,
    type = "TMLE"
  )
  class(res) <- c("vibr.fit", class(res))
  res
}

.varimp_tmle_boot <- function(X,
                         Y,
                         delta=0.1,
                         Y_learners=NULL,
                         Xdensity_learners=NULL,
                         Xbinary_learners=NULL,
                         verbose=TRUE,
                         estimand="diff",
                         bounded=FALSE,
                         B=100,
                         ...){
  est <- .varimp_tmle(X,Y,delta,Y_learners,Xdensity_learners,Xbinary_learners,verbose,estimand,bounded,...)
  rn <- rownames(est$res)
  bootests <- matrix(NA, nrow=B, ncol = length(rn))
  n = length(Y)
  isbin <- as.character((length(unique(Y))==2))
  for(b in 1:B){
    if(verbose) cat(".") # TODO: better interpretation
    ridx <- sample(seq_len(n), n, replace=TRUE)
    Xi = X[ridx,,drop=FALSE]
    Yi = Y[ridx]
    tasklist = .create_tasks(Xi,Yi,delta)
    yb = .bound_zero_one(Yi)
    Ybound = yb[[1]]
    sl.qfit <- .train_Y(Xi,Yi, Y_learners, verbose=FALSE, isbin)
    sl.gfits <- .train_allX(Xi, tasklist$slX, Xbinary_learners, Xdensity_learners, verbose=FALSE)
    fittable <- .EstEqTMLE(n,Xi,Yi,delta,qfun=.qfunction,gfun=.gfunction,qfit=sl.qfit,gfits=sl.gfits, estimand,bounded, ...)
    bootests[b,] <- fittable$est
  }
  colnames(bootests) <- rn
  if(verbose) cat("\n")
  res <- list(
    est = est,
    boots = bootests,
    binomial = isbin,
    type = "TMLE"
  )
  class(res) <- c("vibr.bootfit", class(res))
  res
}
