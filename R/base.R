# user facing functions

#' Variable importance in cross-sectional data
#'
#' @param X data frame of predictors
#' @param Y outcome
#' @param V (default NULL) a data frame of other variables that contain weights or offsets, if used ("weights" and "offset" are both passed to other functions via extra arguments represented by ...: see example below)
#' @param delta change in each column of X corresponding to
#' @param Y_learners list of sl3 learners used to predict the outcome, conditional on all predictors in X
#' @param Xdensity_learners list of sl3 learners used to estimate the density of continuous predictors, conditional on all other predictors in X
#' @param Xbinary_learners list of sl3 learners used to estimate the probability mass of continuous predictors, conditional on all other predictors in X
#' @param verbose (logical) print extra information
#' @param estimator (character) "AIPW" (default), "TMLE", "GCOMP", "IPW"
#' @param bounded (logical) not used
#' @param updatetype (character) (used for estimator = "TMLE" only) weighted" or any other valid character. If "weighted" then uses weighting by clever covariate in update step of TMLE, otherwise fits a generalized linear model with no intercept and clever covariate as a sole predictor
#' @param estimand (character) "diff" (default, estimate mean difference comparing Y under intervention with observed Y), "mean" (estimate mean Y under intervention)
#' @param B (NULL or integer) Numer of bootstrap iterations (NULL = asymptotic variance only)
#' @param showProgress show progress of bootstrapping (only relevant if B is not NULL)
#' @param ... passed to sl3::make_sl3_Task (e.g. weights)
#'
#' @return vi object
#' @export
#' @examples
#' \dontrun{
#' data(metals, package="qgcomp")
#' XYlist = list(X=metals[,1:23], Y=metals$y)
#' Y_learners = .default_continuous_learners()
#' Xbinary_learners = list(Lrnr_stepwise$new(name="SW"))
#' Xdensity_learners = .default_density_learners(n_bins=c(10))
#' vi <- varimp(X=XYlist$X,Y=XYlist$Y, delta=0.1, Y_learners = Y_learners,
#'        Xdensity_learners=Xdensity_learners[1:2], Xbinary_learners=Xbinary_learners,
#'        estimator="TMLE")
#' vi
#' V = data.frame(wt=runif(nrow(metals)))
#' viw <- varimp(X=XYlist$X,Y=XYlist$Y, V=V, delta=0.1, Y_learners = Y_learners,
#'        Xdensity_learners=Xdensity_learners[1:2], Xbinary_learners=Xbinary_learners,
#'        estimator="TMLE", weights="wt")
#' viw
#' }
varimp <- function(X,
                   Y,
                   V = NULL,
                   delta=0.1,
                   Y_learners=NULL,
                   Xdensity_learners=NULL,
                   Xbinary_learners=NULL,
                   verbose=TRUE,
                   estimator="AIPW",
                   bounded=FALSE,
                   updatetype="weighted",
                   estimand="diff",
                   B=NULL,
                   showProgress=TRUE,
                   scale_continuous = TRUE,
                   ...){
  if(scale_continuous){
    if(verbose) cat("Scaling all continuous variables by 2*sd\n")
    # divide continuous by 2*sd
    X = .scale_continuous(X)
  }
  if(is.null(Y_learners)) Y_learners = .default_continuous_learners()
  if(is.null(Xbinary_learners)) Xbinary_learners = .default_binary_learners()
  if(is.null(Xdensity_learners)) Xdensity_learners = .default_density_learners(n_bins=c(5, 20))
  if(is.null(B)){
    res = switch(estimator,
                 AIPW=.varimp_aipw(  X=X,Y=Y,V=V, delta=delta, Y_learners=Y_learners, Xdensity_learners=Xdensity_learners, Xbinary_learners=Xbinary_learners, verbose=verbose, estimand=estimand, bounded=bounded,...),
                 GCOMP=.varimp_gcomp(X=X,Y=Y,V=V, delta=delta, Y_learners=Y_learners, Xdensity_learners=Xdensity_learners, Xbinary_learners=Xbinary_learners, verbose=verbose, estimand=estimand, bounded=bounded,...),
                 IPW=.varimp_ipw(    X=X,Y=Y,V=V, delta=delta, Y_learners=Y_learners, Xdensity_learners=Xdensity_learners, Xbinary_learners=Xbinary_learners, verbose=verbose, estimand=estimand, bounded=bounded,...),
                 TMLE=.varimp_tmle(  X=X,Y=Y,V=V, delta=delta, Y_learners=Y_learners, Xdensity_learners=Xdensity_learners, Xbinary_learners=Xbinary_learners, verbose=verbose, estimand=estimand, bounded=bounded, updatetype=updatetype,...)
    )
  } else{
    res = switch(estimator,
                 AIPW=.varimp_aipw_boot(  X=X,Y=Y,V=V, delta=delta, Y_learners=Y_learners, Xdensity_learners=Xdensity_learners, Xbinary_learners=Xbinary_learners, verbose=verbose, estimand=estimand, bounded=bounded, B=B, showProgress=showProgress,...),
                 GCOMP=.varimp_gcomp_boot(X=X,Y=Y,V=V, delta=delta, Y_learners=Y_learners, Xdensity_learners=Xdensity_learners, Xbinary_learners=Xbinary_learners, verbose=verbose, estimand=estimand, bounded=bounded, B=B, showProgress=showProgress,...),
                 IPW=.varimp_ipw_boot(    X=X,Y=Y,V=V, delta=delta, Y_learners=Y_learners, Xdensity_learners=Xdensity_learners, Xbinary_learners=Xbinary_learners, verbose=verbose, estimand=estimand, bounded=bounded, B=B, showProgress=showProgress,...),
                 TMLE=.varimp_tmle_boot(  X=X,Y=Y,V=V, delta=delta, Y_learners=Y_learners, Xdensity_learners=Xdensity_learners, Xbinary_learners=Xbinary_learners, verbose=verbose, estimand=estimand, bounded=bounded, updatetype=updatetype, B=B, showProgress=showProgress,...)
    )
    res$est <- .attach_misc(res$est, scale_continuous=scale_continuous, delta=delta)
  }
  res <- .attach_misc(res, scale_continuous=scale_continuous, delta=delta)
  res
}


#' Variable importance in cross-sectional data
#'
#' @description Refitting a variable importance model with outcome regression and propensity scores trained from another model
#' @param vibr.fit a vibr.fit object
#' @param X data frame of predictors
#' @param Y outcome
#' @param delta change in each column of X corresponding to
#' @param verbose (logical) print extra information
#' @param estimator (character) "AIPW" (default), "TMLE", "GCOMP", "IPW"
#' @param bounded (logical) not used
#' @param updatetype (character) "weighted" or any other valid character. If "weighted" then uses weighting by clever covariate in update step of TMLE, otherwise fits a generalized linear model with no intercept and clever covariate as a fixed effect.
#' @param estimand (character) "diff" (default, estimate mean difference comparing Y under intervention with observed Y), "mean" (estimate mean Y under intervention)
#'
#' @return vibr.fit object
#' @export
#' @examples
#' \dontrun{
#' data(metals, package="qgcomp")
#' XYlist = list(X=metals[,c(1:22, 23)], Y=metals$y)
#' Y_learners = .default_continuous_learners_big()
#' Xbinary_learners = .default_binary_learners_big()
#' Xdensity_learners = .default_density_learners_big()[c(1:4,6:7)]
#' vi <- varimp(X=XYlist$X,Y=XYlist$Y, delta=0.1, Y_learners = Y_learners,
#'        Xdensity_learners=Xdensity_learners, Xbinary_learners=Xbinary_learners,
#'        estimator="TMLE", updatetype="unweighted",estimand="diff")
#' vi
#' vi1 <- varimp_refit(vi, X=XYlist$X,Y=XYlist$Y, delta=0.1,
#'                     estimator="TMLE", updatetype="weighted", estimand="diff")
#' vi1
#' vi2 <- varimp_refit(vi, X=XYlist$X,Y=XYlist$Y, delta=0.1,
#'                     estimator="AIPW")
#' vi2
#' vi3 <- varimp_refit(vi, X=XYlist$X,Y=XYlist$Y, delta=0.1,
#'                     estimator="GCOMP", estimand="mean")
#' vi3
#' vi4 <- varimp_refit(vi, X=XYlist$X,Y=XYlist$Y, delta=0.1,
#'                     estimator="IPW")
#' vi4
#'
#' hist(metals$y)
#' hist(metals$calcium)
#' hist(metals$total_hardness)
#' caidx <- which(names(XYlist$X)=="calcium")
#' plot(metals$calcium, vi1$gfits[[4]]$predict()[[1]], pch=19, cex=0.2)
#' plot(metals$calcium, metals$y, pch=19, cex=0.2)
#' plot(metals$total_hardness, vi1$gfits[[21]]$predict()[[1]], pch=19, cex=0.2)
#' plot(metals$total_hardness, metals$y, pch=19, cex=0.2)
#' }
varimp_refit <- function(vibr.fit,
                         X,
                         Y,
                         delta=0.1,
                         verbose=TRUE,
                         estimator="AIPW",
                         bounded=FALSE,
                         updatetype="weighted",
                         estimand="diff"
                         ){
  obj <- vibr:::.ChExstractFit(vibr.fit)
  if(obj$scaled){
    if(verbose) cat("Scaling all continuous variables by 2*sd\n")
    # divide continuous by 2*sd
    X = .scale_continuous(X)
  }
  #list(yfit=yfit, xfit=xfit, sl.qfit=obj$qfit, sl.gfits=obj$gfit,
  #     tasklist=tasklist, binomial=obj$binomial,oldtype=obj$type)
  if(!is.na(match(estimator, c("TMLE", "AIPW")))){
    stopifnot(obj$yfit & obj$xfit)
    stopifnot(length(obj$sl.gfits) == ncol(X))
  }
  if(!is.na(match(estimator, c("GCOMP")))) stopifnot(obj$yfit)
  if(!is.na(match(estimator, c("IPW")))) stopifnot(obj$xfit)

  res = switch(estimator,
               IPW =     .trained_ipw(obj=obj,X=X,Y=Y,delta=delta,qfun=.qfunction,gfun=.gfunction,estimand=estimand,bounded=bounded,updatetype=updatetype),
               GCOMP = .trained_gcomp(obj=obj,X=X,Y=Y,delta=delta,qfun=.qfunction,gfun=.gfunction,estimand=estimand,bounded=bounded,updatetype=updatetype),
               AIPW =   .trained_aipw(obj=obj,X=X,Y=Y,delta=delta,qfun=.qfunction,gfun=.gfunction,estimand=estimand,bounded=bounded,updatetype=updatetype),
               TMLE =  .trained_tmle(obj=obj,X=X,Y=Y,delta=delta,qfun=.qfunction,gfun=.gfunction,estimand=estimand,bounded=bounded,updatetype=updatetype)
  )
  res <- .attach_misc(res, scale_continuous=obj$scaled, delta=delta)
  res
}




#' @importFrom stats printCoefmat
#' @export
print.vibr.fit <- function(x, ...){
  xr = x$res
  xr$rnk = rank(-abs(xr$est))
  cat(paste0("Variable importance estimates (", x$type, "): \n"))
  names(xr) <- c("Estimate", "Std. Error (asymptotic)", "z value", "Pr(>|z|)", "Rank")
  printCoefmat(xr, P.values=TRUE, has.Pvalue=TRUE, signif.stars=FALSE, cs.ind=c(1,2), tst.ind=3)
  if(x$scaled){
    cat("\nScaling:\n")
    cat(paste0("  Continuous predictors are scaled by 2*Std. Dev\n"))
  }
  invisible(x)
}




#' @importFrom stats printCoefmat
#' @export
print.vibr.bootfit <- function(x, ...){
  asest = x$est
  print(asest)
  cat("\n")
  #
  est = asest$res$est
  rnk = rank(-abs(est))
  sds <- apply(x$boots,2,sd)
  zz <- est/sds
  xr2 <- as.data.frame(cbind(est, sds, zz, pnorm(-abs(zz))*2, rnk))
  names(xr2) <- c("Estimate", "Std. Error (bootstrap)", "z value", "Pr(>|z|)", "Rank")
  printCoefmat(xr2, P.values=TRUE, has.Pvalue=TRUE, signif.stars=FALSE, cs.ind=c(1,2), tst.ind=3)
  invisible(x)
}

