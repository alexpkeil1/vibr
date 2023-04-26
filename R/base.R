# user facing functions

#' Stochastic intervention based variable importance in cross-sectional data
#'
#' @description Estimate variable importance based on a stochastic intervention to increase each predictor by a small amount (continuous) or increase the probability of a predictor by the same small amount (binary). The underlying methodology is based on papers by Ivan Diaz and Mark van der Laan. This function supports doubly robust estimation (Targeted maximum likelihood and augmented inverse probability weighting) as well as modern techniques to reduce bias and speed up convergence (in sample size) of these methods (double cross-fitting and repeated double cross-fitting). The underlying statistical models draw on \code{sl3}, an R package with a unified framework for machine learning and ensemble machine learning based around the super learner algorithm (stacking).
#' @param X data frame of variables for which variable importance will be estimated
#' @param W data frame of covariates (e.g. potential confounders) for which variable importance will not be estimated
#' @param Y outcome
#' @param V (default NULL) a data frame of other variables that contain weights or offsets, if used ("weights" and "offset" are both passed to other functions via extra arguments represented by ...: see example below)
#' @param delta change in each column of X corresponding to
#' @param Y_learners list of sl3 learners used to predict the outcome, conditional on all predictors in X
#' @param Xdensity_learners list of sl3 learners used to estimate the density of continuous predictors, conditional on all other predictors in X
#' @param Xbinary_learners list of sl3 learners used to estimate the probability mass of continuous predictors, conditional on all other predictors in X
#' @param verbose (logical) print extra information
#' @param estimator (character) "AIPW" (default), "TMLE", "GCOMP", "IPW", "TMLEX" (cross fit TMLE), "AIPWX" (cross fit AIPW), "GCOMPX" (cross fit GCOMP), "IPWX" (cross fit IPW)
#' @param bounded (logical) not yet implemented
#' @param updatetype (character) (used for estimator = "TMLE" only) "weighted" or "predictor." If "weighted" then uses weighting by clever covariate in update step of TMLE, otherwise fits a generalized linear model with no intercept and clever covariate as a sole predictor
#' @param estimand (character) "diff" (default, estimate mean difference comparing Y under intervention with observed Y), "mean" (estimate mean Y under intervention)
#' @param family (character or glm families binomial or gaussian, default = gaussian()) Outcome type: can be gaussian(), binomial(), "gaussian" or "binomial"; will be guessed if left NULL
#' @param xfitfolds (odd integer, default=3) (used for estimator = "TMLEX" only) number of cross-fit folds (must be odd number - last fold is used for validation while the rest of the data are split in two for fitting treatment or outcome models)
#' @param foldrepeats (integer, default=10) (used for estimator = "TMLEX" only) number of times to repeat cross-fitting (higher numbers = more stable)
#' @param B (NULL or integer) Number of bootstrap iterations (NULL = asymptotic variance only)
#' @param showProgress show progress of bootstrapping (only relevant if B is not NULL)
#' @param scale_continuous (logical, default: TRUE) scale all continuous variables in X to have a standard deviation of 0.5
#' @param ... passed to \code{sl3::make_sl3_Task} (e.g. weights)
#'
#' @return vibr_fit object
#' @export
#' @importFrom stats pnorm gaussian binomial lm median quantile rnorm sd
#' @importFrom utils combn
#' @examples
#' \dontrun{
#' data(metals, package="qgcomp")
#' XYlist = list(X=metals[,1:23], Y=metals$y)
#' Y_learners = .default_continuous_learners()
#' Xbinary_learners = list(Lrnr_stepwise$new(name="SW"))
#' Xdensity_learners = .default_density_learners(n_bins=c(10))
#' set.seed(1231)
#' vi <- varimp(X=XYlist$X,Y=XYlist$Y, delta=0.1, Y_learners = Y_learners,
#'        Xdensity_learners=Xdensity_learners[1:2], Xbinary_learners=Xbinary_learners,
#'        estimator="TMLE")
#' vi
#' set.seed(1231)
#' V = data.frame(wt=runif(nrow(metals)))
#' viw <- varimp(X=XYlist$X,Y=XYlist$Y, V=V, delta=0.1, Y_learners = Y_learners,
#'        Xdensity_learners=Xdensity_learners[1:2], Xbinary_learners=Xbinary_learners,
#'        estimator="TMLE", weights="wt")
#' viw
#' }
varimp <- function(X,
                   W = NULL,
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
                   family = gaussian(),
                   xfitfolds=3,
                   foldrepeats=10,
                   B=NULL,
                   showProgress=TRUE,
                   scale_continuous = TRUE,
                   ...){
  if(is.null(family)) {
    isbin <- as.logical((length(unique(Y))==2))
  } else{
    fam <- glm(c(1,0) ~ 1, family=family)$family
    if(fam$family=="binomial"){
      isbin=TRUE
    } else if(fam$family=="gaussian"){
      isbin=FALSE
    } else{
      stop("only gaussian (continuous) and binomial (binary) allowed for 'family'")
    }
  }
  if(scale_continuous){
    if(verbose) cat("Scaling all continuous variables in X by 2*sd\n")
    # divide continuous by 2*sd
    X = .scale_continuous(X)
  }
  whichcols = seq_len(ncol(X))
  if(!is.null(W)){
    X = data.frame(X,W)
  }
  if(is.null(Y_learners)) Y_learners = .default_continuous_learners()
  if(is.null(Xbinary_learners)) Xbinary_learners = .default_binary_learners()
  if(is.null(Xdensity_learners)) Xdensity_learners = .default_density_learners(n_bins=c(5, 20))
  if(is.null(B)){
    res = switch(estimator,
                 # standard estimators
                 IPW=.varimp_ipw(          X=X,Y=Y,V=V, whichcols=whichcols, delta=delta, Y_learners=Y_learners, Xdensity_learners=Xdensity_learners, Xbinary_learners=Xbinary_learners, verbose=verbose, estimand=estimand, bounded=bounded, isbin=isbin,...),
                 GCOMP=.varimp_gcomp(      X=X,Y=Y,V=V, whichcols=whichcols, delta=delta, Y_learners=Y_learners, Xdensity_learners=Xdensity_learners, Xbinary_learners=Xbinary_learners, verbose=verbose, estimand=estimand, bounded=bounded, isbin=isbin,...),
                 AIPW=.varimp_aipw(        X=X,Y=Y,V=V, whichcols=whichcols, delta=delta, Y_learners=Y_learners, Xdensity_learners=Xdensity_learners, Xbinary_learners=Xbinary_learners, verbose=verbose, estimand=estimand, bounded=bounded, isbin=isbin,...),
                 TMLE=.varimp_tmle(        X=X,Y=Y,V=V, whichcols=whichcols, delta=delta, Y_learners=Y_learners, Xdensity_learners=Xdensity_learners, Xbinary_learners=Xbinary_learners, verbose=verbose, estimand=estimand, bounded=bounded, isbin=isbin, updatetype=updatetype,...),
                 # cross fit estimators
                 IPWX=.varimp_ipw_xfit(    X=X,Y=Y,V=V, whichcols=whichcols, delta=delta, Y_learners=Y_learners, Xdensity_learners=Xdensity_learners, Xbinary_learners=Xbinary_learners, verbose=verbose, estimand=estimand, bounded=bounded, isbin=isbin, xfitfolds=xfitfolds,foldrepeats=foldrepeats,...),
                 GCOMPX=.varimp_gcomp_xfit(X=X,Y=Y,V=V, whichcols=whichcols, delta=delta, Y_learners=Y_learners, Xdensity_learners=Xdensity_learners, Xbinary_learners=Xbinary_learners, verbose=verbose, estimand=estimand, bounded=bounded, isbin=isbin, xfitfolds=xfitfolds,foldrepeats=foldrepeats,...),
                 AIPWX=.varimp_aipw_xfit(  X=X,Y=Y,V=V, whichcols=whichcols, delta=delta, Y_learners=Y_learners, Xdensity_learners=Xdensity_learners, Xbinary_learners=Xbinary_learners, verbose=verbose, estimand=estimand, bounded=bounded, isbin=isbin, xfitfolds=xfitfolds,foldrepeats=foldrepeats,...),
                 TMLEX=.varimp_tmle_xfit(  X=X,Y=Y,V=V, whichcols=whichcols, delta=delta, Y_learners=Y_learners, Xdensity_learners=Xdensity_learners, Xbinary_learners=Xbinary_learners, verbose=verbose, estimand=estimand, bounded=bounded, isbin=isbin, updatetype=updatetype, xfitfolds=xfitfolds,foldrepeats=foldrepeats,...)
    )
  } else{
    res = switch(estimator,
                 # standard estimators
                 IPW=.varimp_ipw_boot(          X=X,Y=Y,V=V, whichcols=whichcols, delta=delta, Y_learners=Y_learners, Xdensity_learners=Xdensity_learners, Xbinary_learners=Xbinary_learners, verbose=verbose, estimand=estimand, bounded=bounded, isbin=isbin, B=B, showProgress=showProgress,...),
                 GCOMP=.varimp_gcomp_boot(      X=X,Y=Y,V=V, whichcols=whichcols, delta=delta, Y_learners=Y_learners, Xdensity_learners=Xdensity_learners, Xbinary_learners=Xbinary_learners, verbose=verbose, estimand=estimand, bounded=bounded, isbin=isbin, B=B, showProgress=showProgress,...),
                 AIPW=.varimp_aipw_boot(        X=X,Y=Y,V=V, whichcols=whichcols, delta=delta, Y_learners=Y_learners, Xdensity_learners=Xdensity_learners, Xbinary_learners=Xbinary_learners, verbose=verbose, estimand=estimand, bounded=bounded, isbin=isbin, B=B, showProgress=showProgress,...),
                 TMLE=.varimp_tmle_boot(        X=X,Y=Y,V=V, whichcols=whichcols, delta=delta, Y_learners=Y_learners, Xdensity_learners=Xdensity_learners, Xbinary_learners=Xbinary_learners, verbose=verbose, estimand=estimand, bounded=bounded, isbin=isbin, updatetype=updatetype, B=B, showProgress=showProgress,...),
                 # cross fit estimators
                 IPWX=.varimp_ipw_xfit_boot(    X=X,Y=Y,V=V, whichcols=whichcols, delta=delta, Y_learners=Y_learners, Xdensity_learners=Xdensity_learners, Xbinary_learners=Xbinary_learners, verbose=verbose, estimand=estimand, bounded=bounded, isbin=isbin, xfitfolds=xfitfolds,foldrepeats=foldrepeats, B=B, showProgress=showProgress,...),
                 GCOMPX=.varimp_gcomp_xfit_boot(X=X,Y=Y,V=V, whichcols=whichcols, delta=delta, Y_learners=Y_learners, Xdensity_learners=Xdensity_learners, Xbinary_learners=Xbinary_learners, verbose=verbose, estimand=estimand, bounded=bounded, isbin=isbin, xfitfolds=xfitfolds,foldrepeats=foldrepeats, B=B, showProgress=showProgress,...),
                 AIPWX=.varimp_aipw_xfit_boot(  X=X,Y=Y,V=V, whichcols=whichcols, delta=delta, Y_learners=Y_learners, Xdensity_learners=Xdensity_learners, Xbinary_learners=Xbinary_learners, verbose=verbose, estimand=estimand, bounded=bounded, isbin=isbin, xfitfolds=xfitfolds,foldrepeats=foldrepeats, B=B, showProgress=showProgress,...),
                 TMLEX=.varimp_tmle_xfit_boot(  X=X,Y=Y,V=V, whichcols=whichcols, delta=delta, Y_learners=Y_learners, Xdensity_learners=Xdensity_learners, Xbinary_learners=Xbinary_learners, verbose=verbose, estimand=estimand, bounded=bounded, isbin=isbin, updatetype=updatetype, xfitfolds=xfitfolds,foldrepeats=foldrepeats, B=B, showProgress=showProgress,...)
    )
    res$est <- .attach_misc(res$est, scale_continuous=scale_continuous, delta=delta, B=NULL, whichcols=whichcols)
  }
  res <- .attach_misc(res, scale_continuous=scale_continuous, delta=delta, B=B, whichcols=whichcols)
  res
}


#' Efficiently refitting certain vibr models
#'
#' @description Refitting a variable importance model with outcome regression and propensity scores trained from another model. Works for any estimator which does not utilize sample splitting (i.e. cross-fit estimators)
#' @param vibr_fit a vibr_fit object
#' @param X data frame of variables for which variable importance will be estimated
#' @param W data frame of covariates (e.g. potential confounders) for which variable importance will not be estimated
#' @param Y outcome
#' @param delta change in each column of X corresponding to
#' @param verbose (logical) print extra information
#' @param estimator (character) "AIPW" (default), "TMLE", "GCOMP", "IPW" (note this function does not permit re-fitting cross-fit estimators)
#' @param bounded (logical) not used
#' @param updatetype (character) "weighted" or any other valid character. If "weighted" then uses weighting by clever covariate in update step of TMLE, otherwise fits a generalized linear model with no intercept and clever covariate as a fixed effect.
#' @param estimand (character) "diff" (default, estimate mean difference comparing Y under intervention with observed Y), "mean" (estimate mean Y under intervention)
#'
#' @return vibr_fit object
#' @export
#' @importFrom stats pnorm gaussian binomial
#' @examples
#'
#' library(future)
#' currplan = plan()
#' plan(multisession) # fit models in parallel
#' data(metals, package="qgcomp")
#' XYlist = list(X=metals[,c(1:10, 15:23)], Y=metals$y)
#' Y_learners = .default_continuous_learners_big()
#' Xbinary_learners = .default_binary_learners_big()
#' Xdensity_learners = .default_density_learners_big()[c(1:4,6:7)]
#' vi <- varimp(X=XYlist$X,Y=XYlist$Y, delta=0.1, Y_learners = Y_learners,
#'        Xdensity_learners=Xdensity_learners, Xbinary_learners=Xbinary_learners,
#'        estimator="TMLE", updatetype="unweighted",estimand="diff")
#' vi
#' plan(currplan) # go back to standard evaluation
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
#' # find the fit corresponding to calcium
#' caidx <- which(names(XYlist$X)=="calcium")
#' thidx <- which(names(XYlist$X)=="total_hardness")
#' # can confirm
#' # vi1$gfits[[caidx]]$training_task$nodes$outcome
#' calpredict = vi1$gfits[[caidx]]$predict()[[1]]
#' thpredict = vi1$gfits[[thidx]]$predict()[[1]]
#' # plot predicted density (not predicted value!) against original value,
#' # compare with kernel density
#' plot(metals$calcium, calpredict/max(calpredict), pch=19, cex=0.2,
#'   ylab="scaled conditional density")
#' lines(density(metals$calcium))
#' plot(metals$total_hardness, thpredict/max(thpredict), pch=19, cex=0.2,
#'   ylab="scaled conditional density")
#' lines(density(metals$total_hardness))
#' # note these are effectively measuring much of the same quantity
#' plot(metals$calcium, metals$total_hardness)
#' plot(calpredict, thpredict)
#'
varimp_refit <- function(vibr_fit,
                         X,
                         W = NULL,
                         Y,
                         delta=0.1,
                         verbose=TRUE,
                         estimator="AIPW",
                         bounded=FALSE,
                         updatetype="weighted",
                         estimand="diff"
                         ){
  obj <- .ChExstractFit(vibr_fit)
  if(obj$scaled){
    if(verbose) cat("Scaling all continuous variables in X by 2*sd\n")
    # divide continuous by 2*sd
    X = .scale_continuous(X)
  }
  #list(yfit=yfit, xfit=xfit, sl.qfit=obj$qfit, sl.gfits=obj$gfit,
  #     tasklist=tasklist, binomial=obj$binomial,oldtype=obj$type)
  if(!is.na(match(estimator, c("TMLE", "AIPW")))){
    stopifnot(obj$yfit & obj$xfit)
    stopifnot(length(obj$sl.gfits) == ncol(X))
  }
  whichcols = obj$whichcols
  if(!is.null(W)){
    X = data.frame(X,W)
  }
  if(!is.na(match(estimator, c("GCOMP")))) stopifnot(obj$yfit)
  if(!is.na(match(estimator, c("IPW")))) stopifnot(obj$xfit)

  res = switch(estimator,
               IPW =     .trained_ipw(obj=obj,X=X,Y=Y,delta=delta,qfun=.qfunction,gfun=.gfunction,estimand=estimand,bounded=bounded,updatetype=updatetype),
               GCOMP = .trained_gcomp(obj=obj,X=X,Y=Y,delta=delta,qfun=.qfunction,gfun=.gfunction,estimand=estimand,bounded=bounded,updatetype=updatetype),
               AIPW =   .trained_aipw(obj=obj,X=X,Y=Y,delta=delta,qfun=.qfunction,gfun=.gfunction,estimand=estimand,bounded=bounded,updatetype=updatetype),
               TMLE =  .trained_tmle(obj=obj,X=X,Y=Y,delta=delta,qfun=.qfunction,gfun=.gfunction,estimand=estimand,bounded=bounded , updatetype=updatetype)
  )
  res <- .attach_misc(res, scale_continuous=obj$scaled, delta=delta, whichcols=obj$whichcols)
  res
}




#' @importFrom stats printCoefmat
#' @export
print.vibr_fit <- function(x, ...){
  xr = x$res
  xr$rnk = x$rank
  xr$rnkz = x$rankz
  cat(paste0("Variable importance estimates (", x$type, "): \n"))
  names(xr) <- c("Estimate", "Std. Error (asymptotic)", "z value", "Pr(>|z|)", "Rank(|estimate|)", "Rank(|z|)")
  if(substr(x$type, 1, 5)=="GCOMP"){
    cat("Note: no valid asymptotic std. error estimates are available for this estimator.\n")
    xr = xr[c(1,5)]
    printCoefmat(xr, P.values=FALSE, has.Pvalue=FALSE, signif.stars=FALSE, cs.ind=c(1))
  } else{
    printCoefmat(xr, P.values=TRUE, has.Pvalue=TRUE, signif.stars=FALSE, cs.ind=c(1,2), tst.ind=3)
  }

  if(x$scaled){
    cat("\nScaling:\n")
    cat(paste0("  Continuous predictors are scaled by 2*Std. Dev\n"))
  }
  invisible(x)
}




#' @importFrom stats printCoefmat
#' @export
print.vibr_bootfit <- function(x, ...){
  #asest = x$est
  print(x$est)
  cat("\n")
  #
  est = x$est$res$est
  rnk = x$est$rank
  rnkz = x$rankz
  sds <- x$sds# apply(x$boots,2,sd)
  zz <- est/sds
  xr2 <- as.data.frame(cbind(est, sds, zz, pnorm(-abs(zz))*2, rnk, rnkz))
  names(xr2) <- c("Estimate", "Std. Error (bootstrap)", "z value", "Pr(>|z|)", "Rank(|estimate|)", "Rank(|z|)")
  printCoefmat(xr2, P.values=TRUE, has.Pvalue=TRUE, signif.stars=FALSE, cs.ind=c(1,2), tst.ind=3)
  invisible(x)
}

