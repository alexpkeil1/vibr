# efficient influence curve + estimate from generalized cross fit
# 3 splits: train IPW, Train GCOMP, fit both

.varimp_aipw_xfit <- function(X,
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
                         xfitfolds=3,
                         foldrepeats=10,
                         ...){
  ee = new.env()
  # todo: ensure that outcome is always typed correctly (isbin should be set locally)
  n = length(Y)
  if(.checkeven(xfitfolds)) stop("xfitfolds must be an odd number >2")
  if(xfitfolds==1) message("xfitfolds = 1 implies averaging over multiple standard AIPW fits (determined by fold repeats), rather than cross fitting")
  allpartitions <- lapply(seq_len(foldrepeats), .xfitsplit,n=n,V=xfitfolds)
  order = list()
  idx = 1
  sd = runif(foldrepeats*xfitfolds, -.Machine$integer.max, .Machine$integer.max)
  for(r in seq_len(foldrepeats)){
    partitions <- allpartitions[[r]]
    for(fold in partitions){ # xfitfolds combinations
      order[[idx]] = paste0("p",r, "_", fold$fold)
      idx = idx + 1
      ee[[paste0("p",r, "_", fold$fold)]] <- future::future( {
        set.seed(sd[idx-1])
        X1 = X[fold$set1,,drop=FALSE]
        Y1 = Y[fold$set1]
        V1 = V[fold$set1,,drop=FALSE]
        X2 = X[fold$set2,,drop=FALSE]
        Y2 = Y[fold$set2]
        V2 = V[fold$set2,,drop=FALSE]
        X3 = X[fold$set3,,drop=FALSE]
        Y3 = Y[fold$set3]
        V3 = V[fold$set3,,drop=FALSE]

        obj_G = .prelims(X=X1, Y=Y1, V=V1, whichcols=whichcols, delta=delta, Y_learners=NULL, Xbinary_learners, Xdensity_learners, verbose=verbose, ...)
        obj_Y = .prelims(X=X2, Y=Y2, V=V2, whichcols=whichcols, delta=delta, Y_learners, Xbinary_learners=NULL, Xdensity_learners=NULL, verbose=verbose, ...)
        obj <- .prelims(X=X3, Y=Y3, V=V3, whichcols=whichcols, delta=delta, Y_learners=NULL, Xbinary_learners=NULL, Xdensity_learners=NULL, verbose=verbose,isbin=isbin, ...)
        obj$sl.qfit = obj_Y$sl.qfit
        obj$sl.gfits = obj_G$sl.gfits
        fittable <- try(
          .EstEqAIPW(n=obj$n,X=X3,Y=Y3,whichcols=obj$whichcols,delta=delta,qfun=.qfunction,gfun=.gfunction,qfit=obj$sl.qfit,gfits=obj$sl.gfits, estimand=estimand,bounded=FALSE,wt=obj$weights,isbin=obj$isbin)
        )
        if(inherits(fittable,"try-error")){
          warning(fittable[1])
          fill=rep(NA,length(obj$whichcols))
          ft <- data.frame(est=fill,se=fill, row.names = names(X)[obj$whichcols])
        } else{
          ft <- fittable[,1:2]
          ft[,2] <- obj$n*ft[,2]^2 # asymptotic variance of sqrt(n)(psi_0 - psi_n)
        }
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
  # check for missing
  if(any(is.na(ests))){
    if(foldrepeats==1){
      stop("Error in weight calculation resulted in missing weights for some variables (try using 'W' for potentially problematic variables)")
    } else{
      whichmiss = which(is.na(ests[,1]))
      nmiss = length(whichmiss)
      if(nmiss==foldrepeats)
        stop("Error in weight calculation for all partitionings (foldrepeats)")
      warning(paste0("vibr: Error in weight calculation for ", nmiss, " of ", foldrepeats, " partitionings (foldrepeats) - these are excluded from calculation"))
      ests = ests[-whichmiss,,drop=FALSE]
      vars = vars[-whichmiss,,drop=FALSE]
    }
  }
  #
  #colnames(ests) <- varnm
  #est <- apply(ests, 2, mean) # median also used
  est <- .safeapply(ests, 2, median)
  names(est) <- varnm

  resid <- .safesweepminus(ests, 2, est, check.margin = FALSE)
  vars <- vars/n + resid^2


  #var <-apply(vars, 2, mean)
  V <- .safeapply(vars, 2, median)
  resmat <- data.frame(est=est, se=sqrt(V), z = est/sqrt(V))
  if(any(dim(resmat)==0)) stop("No valid results were generated as a result of errors during estimation")
  resmat$p <- stats::pnorm(-abs(resmat$z))*2
  res <- list(
    res = resmat,
    qfit = NULL,
    gfits = NULL,
    binomial = NULL,
    type = "AIPWX",
    weights=NULL,
    ests = ests,
    vars = vars
  )
  class(res) <- c("vibr_fit", class(res))
  res

}


.varimp_aipw_xfit_boot <- function(X,
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
                              foldrepeats=10,
                              xfitfolds=5,
                              B=100,
                              showProgress=TRUE,
                              ...){
  est <- .varimp_aipw_xfit(X=X,Y=Y,V=V,whichcols=whichcols,delta,Y_learners=Y_learners,Xdensity_learners=Xdensity_learners,Xbinary_learners=Xbinary_learners,verbose=verbose,estimand=estimand,isbin=isbin,bounded=bounded,
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
      bootfit <- .varimp_aipw_xfit(X=Xi,Y=Yi,V=Vi,whichcols=whichcols,delta=delta,Y_learners=Y_learners,Xdensity_learners=Xdensity_learners,Xbinary_learners,verbose=verbose,estimand=estimand,isbin=isbin,bounded=bounded,foldrepeats=foldrepeats,
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
    type = "AIPWX"
  )
  class(res) <- c("vibr_bootfit", class(res))
  res
}
