# TMLE using approach from Targeted Learning in Data Science (S 14.3), Diaz and van der Laan (2018)
# note this is a much easier method than that given in Diaz and van der Laan (2012)
################################
## efficient influence functions
################################
.enforce_min_dens <- function(x,eps=1e-8){
  ifelse(x<eps, eps, x)
}

# note d(a|w) = A+delta, and d^-1(a|w) = h(a|w) = A-delta

.OneStepTmleCont <- function(Y,
                             Qinit,
                             Qdawinit,
                             Haw,
                             Hdaw,
                             isbin=FALSE,
                             weighted=FALSE,
                             wt=wt,
                             .link=.identity,
                             .ilink=.invidentity
                             ){
  # Q_update(A,W) = Q(A,W) + eps*H(A,W)
  # Q_update(A,W) = expit(logit(Q(A,W)) + eps*H(A,W))
  # 1. estimate epsilon
  if(.link(.3)==.ilink(.3)) fam <- gaussian()
  if(.link(.3)!=.ilink(.3)) fam <- binomial()
  if(weighted){
    #
    epsk <- as.numeric(glm(Y~ offset(.link(Qinit)), weights = pmax(0,Haw*wt), family=fam)$coefficients[1])
    #epsk <- as.numeric(glm(Y~ offset(.link(Qinit)), weights = pmax(0,Haw), family=fam)$coefficients[1])
    # 2. update Qk
    Qk1 <- .ilink(.link(Qinit) + epsk)        # Q(A       | W) + eps*H(A      |W)
    Qdawk1 <- .ilink(.link(Qdawinit) + epsk)  # Q(A+delta | W) + eps*H(A+delta|W)
  } else{
    epsk <- as.numeric(glm(Y~ -1 + offset(.link(Qinit)) + Haw, weights=wt, family=fam)$coefficients[1])
    #epsk <- as.numeric(glm(Y~ -1 + offset(.link(Qinit)) + Haw, family=fam)$coefficients[1])
    # 2. update Qk
    Qk1 <- .ilink(.link(Qinit) + epsk*Haw)        # Q(A       | W) + eps*H(A      |W)
    Qdawk1 <- .ilink(.link(Qdawinit) + epsk*Hdaw) # Q(A+delta | W) + eps*H(A+delta|W)
  }
  cbind(Qk1, Qdawk1)
}

.OneStepTmleBin <- function(Y,
                            Qinit,
                            Q1init,
                            Q0init,
                            Haw,
                            H1,
                            H0,
                            isbin=FALSE,
                            weighted=FALSE,
                            wt=wt,
                            .link=.identity,
                            .ilink=.invidentity
                            ){
  # Q_update(A,W) = Q(A,W) + eps*H(A,W)
  # Q_update(A,W) = expit(logit(Q(A,W)) + eps*H(A,W))
  # 1. estimand epsilon
  if(.link(.3)==.ilink(.3)) fam <- gaussian()
  if(.link(.3)!=.ilink(.3)) fam <- binomial()
  if(weighted){
    epsk <- as.numeric(glm(Y~ offset(.link(Qinit)), weights = pmax(0,Haw*wt), family=fam)$coefficients[1])
    #epsk <- as.numeric(glm(Y~ offset(.link(Qinit)), weights = pmax(0,Haw), family=fam)$coefficients[1])
    # 2. update Qk
    Qk1 <- .ilink(.link(Qinit) + epsk)
    Q1k1 <- .ilink(.link(Q1init) + epsk)
    Q0k1 <- .ilink(.link(Q0init) + epsk)
  } else{
    epsk <- as.numeric(glm(Y~ -1 + offset(.link(Qinit)) + Haw, weights=wt, family=fam)$coefficients[1])
    #epsk <- as.numeric(glm(Y~ -1 + offset(.link(Qinit)) + Haw, family=fam)$coefficients[1])
    # 2. update Qk
    Qk1 <- .ilink(.link(Qinit) + epsk*Haw)
    Q1k1 <- .ilink(.link(Q1init) + epsk*H1)
    Q0k1 <- .ilink(.link(Q0init) + epsk*H0)
  }
  cbind(Qk1, Q1k1, Q0k1)
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
                estimand = "mean",
                bounded,
                wt,
                updatetype="weighted",
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
  qinit = qfun(X, Acol, qfit=qfit)
  qbinit = qfun(Xb, Acol, qfit=qfit)
  #
  #ga = .enforce_min_dens(ga,eps=1e-8)
  Haw = .Haw(gn, ga, gb)  # evaluated at A_i (ga/gn) + I(gb=0)
  Hdaw = .Haw(gb, gn, gbb) # evaluated at d(A_i,W_i)
  #
  QuMat <- .OneStepTmleCont(Y=Y,Qinit=qinit,Qdawinit=qbinit,Haw=Haw,Hdaw=Hdaw,isbin=FALSE, weighted=(updatetype=="weighted"), wt=wt, .link=.identity, .ilink=.invidentity)
  Qupdate <- QuMat[,1,drop=TRUE]     # predict at observed A
  Qdawupdate <- QuMat[,2,drop=TRUE]  # predict at observed A + delta
  #.OneStepTmleBin(Y,qinit,Haw,Hdaw,isbin=FALSE)
  #eqfb = predict(lm(y~., data.frame(y=qfb, X=X[,-Acol])))   # simple linear regression on W to get E_g[Q | W]
  # TODO: FIGURE OUT WEIGHTS!
  eqfb <- 0 # cancels out
  psi <- mean(wt*(Qdawupdate - Y*(estimand != "mean")))
  dc1 <- Haw*(Y - Qupdate)
  dc2 <- Qdawupdate - eqfb
  dc3 <- eqfb - Y*(estimand != "mean") - psi               # Y doesn't show up in Diaz,vdl 2012 b/c they are estimating mean Y|A+delta
  list(eif = as.vector(dc1 + dc2 + dc3)*wt, psi=psi)
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
                    estimand = "mean",
                    bounded,
                    wt,
                    updatetype="weighted",
                    ...
){
  # define shifts
  X0 <- .shift(X,Acol, -X[,Acol])
  X1 <- .shift(X,Acol,  (1-X[,Acol]))
  g0 <- 1-gfun(X,Acol,gfits=gfits)

  #
  qinit = qfun(X, Acol,qfit=qfit)
  q1init = qfun(X1, Acol,qfit=qfit)
  q0init = qfun(X0, Acol,qfit=qfit)
  #
  #ga = .enforce_min_dens(ga,eps=1e-8)

  Hawmat = .Hawb(g0, delta, X, Acol, retcols=3)
  Haw  <- Hawmat[,1]
  H1 <- Hawmat[,2]
  H0 <- Hawmat[,3]

  QuMat <- .OneStepTmleBin(Y=Y,Qinit=qinit, Q1init=q1init, Q0init=q0init,Haw=Haw,H1=H1,H0=H0,isbin=FALSE, weighted=(updatetype=="weighted"), wt=wt,.link=.identity, .ilink=.invidentity)
  Qupdate <- QuMat[,1,drop=TRUE]
  Q1update <- QuMat[,2,drop=TRUE]
  Q0update <- QuMat[,3,drop=TRUE]
  #.OneStepTmleBin(Y,qinit,Haw,Hdaw,isbin=FALSE)
  #eqfb = predict(lm(y~., data.frame(y=qfb, X=X[,-Acol])))   # simple linear regression on W to get E_g[Q | W]
  eqfb <- 0 # cancels out
  psi <- mean(wt*(Qupdate + delta*(Q1update - Q0update) - Y*(estimand != "mean")))
  dc1 <- Haw*(Y - Qupdate)
  dc2 <- Qupdate - eqfb
  dc3 <- delta*(Q1update - Q0update) + eqfb - Y*(estimand != "mean") - psi                # Y doesn't show up in Diaz,vdl 2012 b/c they are estimating mean Y|A+delta
  list(eif = as.vector(dc1 + dc2 + dc3)*wt, psi=psi)
}

################################
#: estimating equations
################################

.MakeTmleEst <- function(dphi){
  est <- dphi$psi
  eif <- dphi$eif
  #summary(fit <- lm(dphi~1))
  #est <- mean(eif)
  D <- eif
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
                       estimand,
                       bounded,
                       wt=rep(1,n),
                       updatetype="weighted",
                       ...){
  isbin_vec <- apply(X, 2, function(x) length(unique(x))==2)
  resmat <- matrix(NA, nrow=length(isbin_vec), ncol=3)
  for(Acol in seq_len(length(isbin_vec))){

    if(isbin_vec[Acol]){
      dphi <- .DbTMLE(n,X,Y,Acol,delta,qfun,gfun,qfit,gfits,estimand,bounded,wt,updatetype)
    } else{
      dphi <- .DcTMLE( n,X,Y,Acol,delta,qfun,gfun,qfit,gfits,estimand,bounded,wt,updatetype)
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
.trained_tmle <- function(obj,
                          X,
                          Y,
                          delta,
                          qfun,
                          gfun,
                          estimand,
                          bounded,
                          updatetype){
  fittable <- .EstEqTMLE(n=obj$n,X=X,Y=Y,delta=delta,qfun=.qfunction,gfun=.gfunction,qfit=obj$sl.qfit,gfits=obj$sl.gfits, estimand=estimand,bounded=bounded,wt=obj$weights,updatetype=updatetype)
  res <- list(
    res = fittable,
    qfit = obj$sl.qfit,
    gfits = obj$sl.gfits,
    binomial = obj$isbin,
    type = "TMLE",
    weights=obj$weights
    )
  class(res) <- c("vibr.fit", class(res))
  res
}

#' @export
.varimp_tmle <- function(X,
                         Y,
                         V=NULL,
                         delta=0.1,
                         Y_learners=NULL,
                         Xdensity_learners=NULL,
                         Xbinary_learners=NULL,
                         verbose=TRUE,
                         estimand,
                         bounded=FALSE,
                         updatetype="weighted",
                         ...){
  obj = .prelims(X=X, Y=Y, V=V, delta=delta, Y_learners, Xbinary_learners, Xdensity_learners, verbose=verbose, ...)
  res = .trained_tmle(obj,X,Y,delta,qfun,gfun,estimand,bounded,updatetype)
  res
}

#' @importFrom future future value
#' @export
.varimp_tmle_boot <- function(X,
                              Y,
                              V=NULL,
                              delta=0.1,
                              Y_learners=NULL,
                              Xdensity_learners=NULL,
                              Xbinary_learners=NULL,
                              verbose=TRUE,
                              estimand="diff",
                              bounded=FALSE,
                              updatetype="weighted",
                              B=100,
                              showProgress=TRUE,
                              ...){
  est <- .varimp_tmle(X,Y,V,delta,Y_learners,Xdensity_learners,Xbinary_learners,verbose,estimand,bounded,updatetype,...)
  rn <- rownames(est$res)
  #bootests <- matrix(NA, nrow=B, ncol = length(rn))
  n = length(Y)
  isbin <- as.character((length(unique(Y))==2))
  ee <- new.env()
  for(b in 1:B){
    #ridx <- sample(seq_len(n), n, replace=TRUE)
    ridx <- .bootsample(n)
    ee[[paste0("iter",b)]] <- future::future( {
      if(showProgress) cat(".")
      Xi = X[ridx,,drop=FALSE]
      Yi = Y[ridx]
      Vi = V[ridx,,drop=FALSE]
      obj = .prelims(X=Xi, Y=Yi, V=Vi, delta, Y_learners, Xbinary_learners, Xdensity_learners, verbose=verbose, ...)
      fittable <- .EstEqTMLE(n=obj$n,X=Xi,Y=Yi,delta=delta,qfun=.qfunction,gfun=.gfunction,qfit=obj$sl.qfit,gfits=obj$sl.gfits, estimand=estimand,bounded=bounded,wt=obj$weights,updatetype=updatetype)
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
    type = "TMLE"
  )
  class(res) <- c("vibr.bootfit", class(res))
  res
}
