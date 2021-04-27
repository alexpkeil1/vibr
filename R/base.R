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
#' @param ... passed to sl3::make_sl3_Task (e.g. weights)
#'
#' @return vi object
#' @export
#'
varimp <- function(X,Y, delta=0.1, Y_learners=NULL, Xdensity_learners=NULL, Xbinary_learners=NULL, verbose=TRUE, estimator="AIPW", ...){
  if(is.null(Y_learners)) Y_learners = .default_continuous_learners()
  if(is.null(Xbinary_learners)) Xbinary_learners = .default_binary_learners()
  if(is.null(Xdensity_learners)) Xdensity_learners = .default_density_learners(n_bins=c(5, 20))

  res = switch(estimator,
    AIPW=.varimp_aipw(X,Y, delta=delta, Y_learners=Y_learners, Xdensity_learners=Xdensity_learners, Xbinary_learners=Xbinary_learners, verbose=verbose,...),
    GCOMP=.varimp_gcomp(X,Y, delta=delta, Y_learners=Y_learners, Xdensity_learners=Xdensity_learners, Xbinary_learners=Xbinary_learners, verbose=verbose,...),
    IPW=.varimp_ipw(X,Y, delta=delta, Y_learners=Y_learners, Xdensity_learners=Xdensity_learners, Xbinary_learners=Xbinary_learners, verbose=verbose,...)
  )
  res
}


#' @importFrom stats printCoefmat
#' @export
print.vibr.fit <- function(x, ...){
  xr = x$res
  names(xr) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
  cat(paste0("Variable importance estimates (", x$type, "): \n"))
  printCoefmat(xr, P.values=TRUE, has.Pvalue=TRUE, signif.stars=FALSE, cs.ind=c(1,2))
  invisible(x)
}

