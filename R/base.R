# user facing functions

#' Variable importance in cross-sectional data
#'
#' @param X data frame of predictors
#' @param Y outcome
#' @param delta change in each column of X corresponding to
#' @param Y_learners list of sl3 learners used to predict the outcome, conditional on all predictors in X
#' @param Xdensity_learners list of sl3 learners used to estimate the density of continuous predictors, conditional on all other predictors in X
#' @param Xbinary_learners list of sl3 learners used to estimate the probability mass of continuous predictors, conditional on all other predictors in X
#' @param verbose (logical) print extra information
#' @param estimator (character) "AIPW" (default), "GCOMP", "IPW"
#' @param estimand (character) "diff" (default, estimate mean difference comparing Y under intervention with observed Y), "mean" (estimate mean Y under intervention)
#' @param B (NULL or integer) Numer of bootstrap iterations (NULL = asymptotic variance only)
#' @param ... passed to sl3::make_sl3_Task (e.g. weights)
#'
#' @return vi object
#' @export
#'
varimp <- function(X,Y, delta=0.1, Y_learners=NULL, Xdensity_learners=NULL, Xbinary_learners=NULL, verbose=TRUE, estimator="AIPW", bounded=FALSE, estimand="diff", B=NULL, ...){
  if(is.null(Y_learners)) Y_learners = .default_continuous_learners()
  if(is.null(Xbinary_learners)) Xbinary_learners = .default_binary_learners()
  if(is.null(Xdensity_learners)) Xdensity_learners = .default_density_learners(n_bins=c(5, 20))
  if(is.null(B)){
    res = switch(estimator,
                 AIPW=.varimp_aipw(X,Y, delta=delta, Y_learners=Y_learners, Xdensity_learners=Xdensity_learners, Xbinary_learners=Xbinary_learners, verbose=verbose, estimand=estimand, bounded=bounded,...),
                 GCOMP=.varimp_gcomp(X,Y, delta=delta, Y_learners=Y_learners, Xdensity_learners=Xdensity_learners, Xbinary_learners=Xbinary_learners, verbose=verbose, estimand=estimand, bounded=bounded,...),
                 IPW=.varimp_ipw(X,Y, delta=delta, Y_learners=Y_learners, Xdensity_learners=Xdensity_learners, Xbinary_learners=Xbinary_learners, verbose=verbose, estimand=estimand, bounded=bounded,...),
                 TMLE=.varimp_tmle(X,Y, delta=delta, Y_learners=Y_learners, Xdensity_learners=Xdensity_learners, Xbinary_learners=Xbinary_learners, verbose=verbose, estimand=estimand, bounded=bounded,...)
    )
  } else{
    res = switch(estimator,
                 AIPW=.varimp_aipw_boot(X,Y, delta=delta, Y_learners=Y_learners, Xdensity_learners=Xdensity_learners, Xbinary_learners=Xbinary_learners, verbose=verbose, estimand=estimand, bounded=bounded, B=B,...),
                 GCOMP=.varimp_gcomp_boot(X,Y, delta=delta, Y_learners=Y_learners, Xdensity_learners=Xdensity_learners, Xbinary_learners=Xbinary_learners, verbose=verbose, estimand=estimand, bounded=bounded, B=B,...),
                 IPW=.varimp_ipw_boot(X,Y, delta=delta, Y_learners=Y_learners, Xdensity_learners=Xdensity_learners, Xbinary_learners=Xbinary_learners, verbose=verbose, estimand=estimand, bounded=bounded, B=B,...),
                 TMLE=.varimp_tmle_boot(X,Y, delta=delta, Y_learners=Y_learners, Xdensity_learners=Xdensity_learners, Xbinary_learners=Xbinary_learners, verbose=verbose, estimand=estimand, bounded=bounded, B=B,...)
    )
  }
  res
}


#' @importFrom stats printCoefmat
#' @export
print.vibr.fit <- function(x, ...){
  xr = x$res
  cat(paste0("Variable importance estimates (", x$type, "): \n"))
  names(xr) <- c("Estimate", "Std. Error (asymptotic)", "z value", "Pr(>|z|)")
  printCoefmat(xr, P.values=TRUE, has.Pvalue=TRUE, signif.stars=FALSE, cs.ind=c(1,2))
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
  sds <- apply(x$boots,2,sd)
  zz <- est/sds
  xr2 <- as.data.frame(cbind(est, sds, zz, pnorm(-abs(zz))*2))
  names(xr2) <- c("Estimate", "Std. Error (bootstrap)", "z value", "Pr(>|z|)")
  printCoefmat(xr2, P.values=TRUE, has.Pvalue=TRUE, signif.stars=FALSE, cs.ind=c(1,2))
  invisible(x)
}

