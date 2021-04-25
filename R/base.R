# user facing functions


#' @export
varimp <- function(X,Y, delta=0.1, Y_learners=NULL, Xdensity_learners=NULL, Xbinary_learners=NULL, verbose=TRUE, estimator="AIPW"){

  res = switch(estimator,
    AIPW=.varimp_aipw(X,Y, delta=delta, Y_learners=Y_learners, Xdensity_learners=Xdensity_learners, Xbinary_learners=Xbinary_learners, verbose=verbose,...),
    GCOMP=.varimp_gcomp(X,Y, delta=delta, Y_learners=Y_learners, Xdensity_learners=Xdensity_learners, Xbinary_learners=Xbinary_learners, verbose=verbose,...),
    IPW=.varimp_ipw(X,Y, delta=delta, Y_learners=Y_learners, Xdensity_learners=Xdensity_learners, Xbinary_learners=Xbinary_learners, verbose=verbose,...)
  )
}
