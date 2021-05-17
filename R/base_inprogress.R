

# .Haw <- function(gn, ga, gb){
#   # I(A<u(w))g(a-delta|w)/g(a|w) + I(a>=u(w)-delta)
#   Haw <- ifelse(gn>0, ga/gn, 0) + as.numeric(gb == 0)
#   Haw
# }

.OneStepTmleContBounded <- function(Y,
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
    epsk <- as.numeric(glm(Y~ offset(.link(Qinit)), weights = Haw*wt, family=fam)$coefficients[1])
  } else{
    epsk <- as.numeric(glm(Y~ -1 + offset(.link(Qinit)) + Haw, weights=wt, family=fam)$coefficients[1])
  }
  # 2. update Qk
  Qk1 <- .ilink(.link(Qinit) + epsk*Haw)
  Qdawk1 <- .ilink(.link(Qdawinit) + epsk*Hdaw)
  cbind(Qk1, Qdawk1)
}

.plotfun <- function(){
  data(metals, package="qgcomp")
  XYlist = list(X=metals[,1:23], Y=metals$y)
  set.seed(123123)
  (vimp <- varimp(data.frame(XYlist$X),XYlist$Y, delta=.1, Y_learners=.default_continuous_learners_big(),
                  Xdensity_learners=.default_density_learners_big(), Xbinary_learners=.default_binary_learners_big(),
                  verbose=FALSE, estimator="TMLE", estimand="diff", updatetype = "unweighted"))
  (vimp2 <- varimp(data.frame(XYlist$X),XYlist$Y, delta=.1, Y_learners=.default_continuous_learners_big(),
                  Xdensity_learners=.default_density_learners_big(), Xbinary_learners=.default_binary_learners_big(),
                  verbose=FALSE, estimator="AIPW", estimand="diff", updatetype = "unweighted"))

  gf <-  vimp$gfits[[2]]
  target <- gf$training_task$Y
  gf$training_task$column_names
  gf$training_task$nodes$outcome
  pred <- gf$predict()[[1]]
  plot(target,pred, cex=.3, pch=19)
  points(target+0.1,pred, cex=.3, pch=19, col="red")

  gf2 <-  vimp$gfits[[2]]
  target2 <- gf2$training_task$Y
  pred2 <- gf2$predict()[[1]]
  # should be a note that p is uniform
  hist(pred2)
  points(target+0.1,pred, cex=.3, pch=19, col="red")

}



.gcjoint <- function(n,
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
                     ...){
  #cat(paste0("column ", names(X)[Acol], ": binary\n"))
  isbin_vec <- apply(X, 2, function(x) length(unique(x))==2)
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
  # E(Y | Abs,Acs) = E(Y|Ab,Ac)  +
  #                  (E(Y|Ab,Acs)-E(Y|Ab,Ac)) +
  #                  (E(Y|Abs,Ac)-E(Y|Ab,Ac))
  #                  (E(Y|Abs,Acs)-E(Y|Ab,Acs))
  #   = b_0 + b_1*Acs + b_2*Abs + b_3*Abs*Acs
  #
  #p1 <- qfun(qfit=qfit,...)           # E(Y|Ab,ac)
  p2 <- qfun(Xb,Acol,qfit=qfit,...)   # E(Y|Ab,ac) + (E(Y|Ab,Acs)-E(Y|Ab,Ac))
  p3 <- delta*(qfun(X1,Acol,qfit=qfit,...) - qfun(X0,Acol,qfit=qfit,...)) # (E(Y|Abs,Ac)-E(Y|Ab,Ac))
  p4 <- delta*(qfun(X1b,Acol,qfit=qfit,...) - qfun(X0b,Acol,qfit=qfit,...)) # (E(Y|Abs,Acs)-E(Y|Ab,Acs))
  (p2 + p3 + p4 - Y*(estimand != "mean"))*wt
}


.EstimatorGcompJoint <- function(n,
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
                                 ...){
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

.trained_gcomp_joint <- function(obj,
                           X,
                           Y,
                           expnms,
                           delta,
                           qfun,
                           gfun,
                           estimand,
                           bounded,
                           updatetype){
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


.varimp_gcomp_joint <- function(X,
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
                          ...){
  obj = .prelims(X, Y, V, delta, Y_learners, Xbinary_learners=NULL, Xdensity_learners=NULL, verbose=verbose, ...)
  res = .trained_gcomp_joint(obj,X,Y,expnms,delta,qfun,gfun,estimand,bounded,updatetype)
  res
}


.varimp_gcomp_joint_boot <- function(X,
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
                               ...){
  est <- .varimp_gcomp_joint(X,Y,V,expnms,delta,Y_learners,Xdensity_learners,Xbinary_learners,verbose,estimand,bounded,...)
  rn <- rownames(est$res)
  n = length(Y)
  isbin <- as.character((length(unique(Y))==2))
  ee <- new.env()
  for(b in 1:B){
    ridx <- .bootsample(n)
    ee[[paste0("iter",b)]] <- future::future( {
      if(showProgress) cat(".")
      Xi = X[ridx,,drop=FALSE]
      Yi = Y[ridx]
      Vi = V[ridx,,drop=FALSE]
      obj = .prelims(X=Xi, Y=Yi, V=Vi, delta, Y_learners, Xbinary_learners, Xdensity_learners, verbose=verbose, ...)
      fittable <- .EstimatorGcompJoint(obj$n,Xi,Yi,expnms,delta,qfun=.qfunction,gfun=NULL,qfit=obj$sl.qfit,gfits=NULL, estimand, bounded,wt=obj$weights)
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
