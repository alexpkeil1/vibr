.make_xfit_folds <- function(fold, set1, set2, set3){
  fold <- list(fold = fold,
               set1 = set1,
               set2 = set2,
               set3 = set3
               )
  fold
}


.xfitfolds_from_foldvec <- function (r, folds, ordermat)
{
  nfolds <- length(unique(folds))
  remfolds = (nfolds-1)/2
  set1 <- which(folds %in% ordermat[1:remfolds,r])
  set2 <- which(folds %in% ordermat[(remfolds+1):(2*remfolds),r])
  set3 <- which(folds == ordermat[nfolds,r])
  .make_xfit_folds(r, set1, set2, set3)
}

# cross fitting (nothing here yet)
.xfitsplit <- function(r=1,n, V=5){
  folds <- rep(seq_len(V), length = n)
  folds <- sample(folds)
  combinations <- combn(V,V-1)
  combinations <- rbind(combinations, apply(combinations, 2, function(x) setdiff(1:V,x)))
  lapply(1:V, .xfitfolds_from_foldvec, folds=folds, ordermat=combinations)
}

.checkeven <- function(val){
  !as.logical(val %% 2)
}
# efficient influence curve + estimate from generalized cross fit
# 3 splits: train IPW, Train GCOMP, fit both

.varimp_tmle_xfit <- function(X,
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
                         xfitfolds=3,
                         foldrepeats=10,
                         ...){
  ee = new.env()
  # todo: ensure that outcome is always typed correctly (isbin should be set locally)
  n = length(Y)
  if(.checkeven(xfitfolds) || xfitfolds < 3) stop("xfitfolds must be an odd number >2")
  allpartitions <- lapply(seq_len(foldrepeats), .xfitsplit,n=n,V=xfitfolds)
  order = list()
  idx = 1
  for(r in seq_len(foldrepeats)){
    partitions <- allpartitions[[r]]
    for(fold in partitions){ # xfitfolds combinations
      order[[idx]] = paste0("p",r, "_", fold$fold)
      idx = idx + 1
      ee[[paste0("p",r, "_", fold$fold)]] <- future::future( {
        X1 = X[fold$set1,,drop=FALSE]
        Y1 = Y[fold$set1]
        V1 = V[fold$set1,,drop=FALSE]
        X2 = X[fold$set2,,drop=FALSE]
        Y2 = Y[fold$set2]
        V2 = V[fold$set2,,drop=FALSE]
        X3 = X[fold$set3,,drop=FALSE]
        Y3 = Y[fold$set3]
        V3 = V[fold$set3,,drop=FALSE]

        obj_G = .prelims(X=X1, Y=Y1, V=V1, delta=delta, Y_learners=NULL, Xbinary_learners, Xdensity_learners, verbose=verbose, ...)
        obj_Y = .prelims(X=X2, Y=Y2, V=V2, delta=delta, Y_learners, Xbinary_learners=NULL, Xdensity_learners=NULL, verbose=verbose, ...)
        obj <- .prelims(X=X3, Y=Y3, V=V3, delta=delta, Y_learners=NULL, Xbinary_learners=NULL, Xdensity_learners=NULL, verbose=verbose, ...)
        obj$sl.qfit = obj_Y$sl.qfit
        obj$sl.gfits = obj_G$sl.gfits
        fittable <- .EstEqTMLE(n=obj$n,X=X3,Y=Y3,delta=delta,qfun=.qfunction,gfun=.gfunction,qfit=obj$sl.qfit,gfits=obj$sl.gfits, estimand=estimand,bounded=bounded,wt=obj$weights,updatetype=updatetype)
        ft <- fittable[,1:2]
        ft[,2] <- obj$n*ft[,2]^2 # asymptotic variance of sqrt(n)(psi_0 - psi_n)
        names(ft)[2] <- "sumd2"
        ft
      }, seed=TRUE, lazy=TRUE)
    }
  }
  xfitres <- as.list(future::value(ee))
  ord <- do.call(c, order)
  xfitres <- xfitres[ord]
  varnm <- rownames(xfitres[[1]])
  allests <- do.call(rbind, lapply(xfitres, function(x) x$est))
  allvars <- do.call(rbind, lapply(xfitres, function(x) x$sumd2))
  partitions <- as.numeric(gsub("p([0-9]+)_([0-9]+)", "\\1", names(xfitres)))
  # take means of every three rows for estimates and variances
  ests <- apply(allests, 2, function (x) tapply(x, partitions, mean))
  vars <- apply(allvars, 2, function (x) tapply(x, partitions, mean)) # each partition has sum IF^2 rather than 1/n*sum IF^2
  ##
  #
  colnames(ests) <- varnm
  #est <- apply(ests, 2, mean) # median also used
  est <- apply(ests, 2, median)

  resid <- sweep(ests, 2, est, check.margin = FALSE)
  vars <- vars/n + resid^2


  #var <-apply(vars, 2, mean)
  V <- apply(vars, 2, median)
  resmat <- data.frame(est=est, se=sqrt(V), z = est/sqrt(V))
  resmat$p <- pnorm(-abs(resmat$z))*2
  res <- list(
    res = resmat,
    qfit = NULL,
    gfits = NULL,
    binomial = NULL,
    type = "TMLEX",
    weights=NULL
  )
  class(res) <- c("vibr.fit", class(res))
  res

}


.varimp_tmle_xfit_boot <- function(X,
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
                              foldrepeats=10,
                              xfitfolds=5,
                              B=100,
                              showProgress=TRUE,
                              ...){
  est <- .varimp_tmle_xfit(X=X,Y=Y,V=V,delta,Y_learners=Y_learners,Xdensity_learners=Xdensity_learners,Xbinary_learners=Xbinary_learners,verbose=verbose,estimand=estimand,bounded=bounded,updatetype=updatetype,
                           foldrepeats=foldrepeats,
                           xfitfolds=xfitfolds, ...)
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
      bootfit <- .varimp_tmle_xfit(X=Xi,Y=Yi,V=Vi,delta=delta,Y_learners=Y_learners,Xdensity_learners=Xdensity_learners,Xbinary_learners,verbose=verbose,estimand=estimand,bounded=bounded,updatetype=updatetype,foldrepeats=foldrepeats,
                                   xfitfolds=xfitfolds,...)
      fittable <- bootfit$res
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
    type = "TMLEX"
  )
  class(res) <- c("vibr.bootfit", class(res))
  res
}
