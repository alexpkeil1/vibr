
#' Density based graphical identifiability diagnostics and visualization
#'
#' Plot a sorted estimated density function of pre- and post-stochastic-intervention values of a predictor.
#' The closer the two functions are for pre- and post-intervention will typically
#' mean that the stochastic intervention is better identified for a given column
#' of the predictor matrix. This plot will generally give a useful way to visualize
#' the implied stochastic intervention and may diagnose some issues with the intervention.
#' This can also help diagnose issues with modeling inverse probability weights, when the
#' estimated density is close to zero over a large portion of the data.
#' @param vibr_fit a fit from varimp
#' @param Acol (integer) which column of predictors in call to varimp to diagnose
#' @param delta (numeric, default=0.01) change in each column of predictors in call to varimp corresponding to stochastic intervention
#'
#' @return ggplot2 plot object
#' @export
#' @import ggplot2
#'
#' @examples
#' \dontrun{
#' data(metals, package="qgcomp")
#' # subset of predictors
#' XYlist = list(X=metals[,1:4], Y=metals$y)
#' Y_learners = .default_continuous_learners()
#' Xbinary_learners = list(Lrnr_stepwise$new(name="SW"))
#' Xdensity_learners = .default_density_learners(n_bins=c(10))
#' set.seed(1231)
#' vi_ipw <- varimp(X=XYlist$X,Y=XYlist$Y, delta=0.1, Y_learners = Y_learners,
#'        Xdensity_learners=Xdensity_learners[1:2], Xbinary_learners=Xbinary_learners,
#'        estimator="TMLE")
#' plotshift_dens(vi_ipw, Acol=1, delta=0.01)
#' }
plotshift_dens <- function(vibr_fit, Acol=1, delta){
  if(is.null(delta)) delta <- vibr_fit$delta
  if(!is.null(vibr_fit$qfit)){
    task <- vibr_fit$qfit$training_task
    varnms <- task$nodes$covariates
  }else{
    task <- vibr_fit$gfits[[1]]$training_task
    varnms <- c(task$nodes$outcome, task$nodes$covariates)
  }
  dat <- task$data
  X <- data.frame(dat)[,varnms,drop=FALSE]
  ft <- vibr_fit$gfits[[Acol[1]]]
  xnm = names(X)
  Xc <- .shift(X,Acol,shift = -delta)
  X$set="obs"
  Xc$set="int"
  X2 <- as.data.frame(rbind(X,Xc))
  tsk <- sl3_Task$new(data = X,
                      covariates = xnm[-Acol[1]],
                      outcome = xnm[Acol[1]])
  tskc <- sl3_Task$new(data = Xc,
                       covariates = xnm[-Acol[1]],
                       outcome = xnm[Acol[1]])
  X$dens = ft$predict(tsk)[[1]]
  X$densshift = ft$predict(tskc)[[1]]
  xo=order(X$dens)
  xo2=order(X$densshift)
  X1 = X[xo,]
  X2 = X[xo2,]
  X1$ord = order(X1$dens)
  X2$ord = order(X2$densshift)

  p1 <-
  ggplot() + theme_classic() + scale_color_grey(name="") +
    geom_step(data=X1,aes(x=ord, y=dens, color="Observed"))+
    geom_step(data=X2,aes(x=ord, y=densshift, color="Shifted"))+
    scale_y_continuous(name=paste0("density(",xnm[Acol[1]],")"), expand = expansion(0))+
    scale_x_continuous(name="Sorted index", expand = expansion(.01))
  print(p1)
  invisible(p1)
}

#' Weight-based graphical identifiability diagnostics and visualization (continuous predictors only)
#'
#' This produces a plot of estimated density (horizontal axis) by a ratio of estimated density under a stochastic intervention (intervention exposure = observed exposure -delta) over the estimated density of the observed data. This gives an indication of influential observations (high weights) and allows one to cross-check against estimated density in the observed data.
#'
#' Note that the "weights" used in vibr use a more difficult calculation than the weight given here (the actual weights used by \code{varimp} depend on whether the stochastic intervention pushes exposure outside of the support of the data), so this function is useful for identifying potentially influential points that may turn out to be a non-issue in analysis.
#' @param vibr_fit a vibr_fit object from varimp
#' @param Acol (integer) which column of predictors in call to varimp to diagnose (can only be continuous column of predictors in call to varimp)
#' @param delta (numeric, default=0.01) change in each column of predictors in call to varimp corresponding to stochastic intervention
#'
#' @export
#' @import ggplot2
#'
#' @examples
#' \dontrun{
#' data(metals, package="qgcomp")
#' # subset of predictors
#' XYlist = list(X=metals[,1:4], Y=metals$y)
#' Y_learners = .default_continuous_learners()
#' Xbinary_learners = list(Lrnr_stepwise$new(name="SW"))
#' Xdensity_learners = .default_density_learners(n_bins=c(10))
#' set.seed(1231)
#' vi_ipw <- varimp(X=XYlist$X,Y=XYlist$Y, delta=0.1, Y_learners = Y_learners,
#'        Xdensity_learners=Xdensity_learners[1:2], Xbinary_learners=Xbinary_learners,
#'        estimator="TMLE")
#' plotshift_wt(vi_ipw, Acol=3, delta=0.01)
#' # Note many weights are close to 1.0, meaning that the intervention effect
#' # (in terms of weighted estimators) is based on the relatively smaller number
#' # of observations with weights > 1.0. This could be reflected by a large
#' # overlap of the observed and "shifted" density functions, which is reflected
#' # by the call to plotshift_dens
#' plotshift_dens(vi_ipw, Acol=3, delta=0.01)
#' # There is not necessarily an ideal appearance for plotshift_wt plots, but
#' # a smooth and well estimated density will tend to lead to y-axes where the
#' # values concentrate away from 1.0 but do not get too large. Thus,
#' # this plot can diagnose issues in identifiability (really large weights),
#' # or density estimation (not very many non-zero weights even with large
#' # value of delta, or really large weights)
#' }
plotshift_wt <- function(vibr_fit, Acol=1, delta=0.01){
  if(is.null(delta)) delta <- vibr_fit$delta
  if(!is.null(vibr_fit$qfit)){
    task <- vibr_fit$qfit$training_task
    varnms <- task$nodes$covariates
  }else{
    task <- vibr_fit$gfits[[1]]$training_task
    varnms <- c(task$nodes$outcome, task$nodes$covariates)
  }
  dat <- task$data
  ft <- vibr_fit$gfits[[Acol[1]]]
  X <- data.frame(dat)[,varnms,drop=FALSE]
  xnm = names(X)
  Xc <- .shift(X,Acol,shift = -delta)
  X$set="obs"
  Xc$set="int"
  tsk <- sl3_Task$new(data = X,
                      covariates = xnm[-Acol[1]],
                      outcome = xnm[Acol[1]])
  tskc <- sl3_Task$new(data = Xc,
                      covariates = xnm[-Acol[1]],
                      outcome = xnm[Acol[1]])
  X$dens = ft$predict(tsk)[[1]]
  X$densshift = ft$predict(tskc)[[1]]
  xo=order(X$dens)
  X = X[xo,]
  X$wt = X$densshift/X$dens
  p1 <-
  ggplot(data = X) + theme_classic() + scale_color_grey() +
    geom_point(aes(x=dens, y=wt), pch=19, size=1, alpha=0.5)+
    scale_x_continuous(name=paste0("density(observed ",xnm[Acol[1]],")"))+
    scale_y_continuous(name="density(shifted)/density(observed)")
  print(p1)
  invisible(p1)
}

#' Visualizing stochastic interventions via bivariate plots
#'
#' This helps identify some potential identifiability issues that may not be picked up via weight diagnostics. A bivariate plot (where the X-axis variable is plotted with both the observed and "intervened" values). If the intervention implies predictors that are pushed outside the support of the data, this will manifest as shifted observations that are far from the observed observations. This is a low-dimensional way to roughly check identifiability of interventions, though users should be aware that identifiability is much more problematic with larger dimensions (number of predictors).
#' @param vibr_fit a vibr_fit object from varimp
#' @param Acol (integer) which column of predictors in call to varimp to diagnose
#' @param Bcol (integer) second column of predictors in call to varimp to diagnose
#' @param delta (numeric, default=0.01) change in each column of predictors in call to varimp corresponding to stochastic intervention
#'
#' @export
#' @import ggplot2
#'
#' @examples
#' \dontrun{
#' data(metals, package="qgcomp")
#' # subset of predictors
#' XYlist = list(X=metals[,1:7], Y=metals$y)
#' Y_learners = .default_continuous_learners()
#' Xbinary_learners = list(Lrnr_stepwise$new(name="SW"))
#' Xdensity_learners = .default_density_learners(n_bins=c(10))
#' set.seed(1231)
#' vi_ipw <- varimp(X=XYlist$X,Y=XYlist$Y, delta=0.1, Y_learners = Y_learners,
#'        Xdensity_learners=Xdensity_learners, Xbinary_learners=Xbinary_learners,
#'        estimator="TMLE")
#' plotshift_scatter(vi_ipw, Acol=1, Bcol=2, delta=0.1)
#' # delta = 0.1, seems reasonable in two dimensions, but may lead to
#' extrapolation in more dimensions. Good to check with other "Bcol" variables
#' plotshift_scatter(vi_ipw, Acol=1, Bcol=5, delta=0.1)
#' # contrast with an intervention level that leads to obvious extrapolation
#' plotshift_scatter(vi_ipw, Acol=1, Bcol=2, delta=1)
#' plotshift_scatter(vi_ipw, Acol=1, Bcol=5, delta=1)
#' }
plotshift_scatter <- function(vibr_fit, Acol, Bcol, delta=NULL, joint=FALSE){
  if(is.null(delta)) delta <- vibr_fit$delta
  if(!is.null(vibr_fit$qfit)){
    task <- vibr_fit$qfit$training_task
    varnms <- task$nodes$covariates
  }else{
    task <- vibr_fit$gfits[[1]]$training_task
    varnms <- c(task$nodes$outcome, task$nodes$covariates)
  }
  dat <- task$data
  X <- data.frame(dat)[,varnms,drop=FALSE]
  #
  requireNamespace("ggplot2")
  Xint <- data.frame(
    x1a = X[,Acol]+delta,
    x2a = X[,Bcol]+ifelse(joint, delta, 0)
  )
  X$col <- "Observed"
  Xint$col <- "Shifted"

  p1 <- ggplot() + theme_classic() +
    geom_point(data = Xint, aes_string(x="x1a",y="x2a", color="col"), size=1) +
    geom_point(data = X,    aes_string(x=varnms[Acol],y=varnms[Bcol], color="col"), size=1) +
    scale_x_continuous(name=varnms[Acol]) +
    scale_y_continuous(name=varnms[Bcol]) +
    scale_color_grey(name = "", start=0, end=.7)+
    theme(legend.position = c(.99,0.01), legend.justification = c(1,0))
  print(p1)
  invisible(p1)
}



#' Density based numerical identifiability diagnostics
#'
#' Give a numerical (but cruder) version of the diagnostics in plotshift_dens, where one can track the change in estimated exposure mass/density following a stochastic intervention on exposure.
#' @param vibr_fit a fit from varimp
#' @param X predictors from a varimp fit
#' @param Acol (integer) which column of predictors in call to varimp to diagnose
#' @param delta (numeric, default=0.01) change in each column of predictors in call to varimp corresponding to stochastic intervention
#' @param quantiles (numeric vector, default=c(0, 0.1, 0.9, 1)) cutpoints in the closed interval [0,1] that correspond to quantiles of the estimated density of observed values of a predictor. The length of this vector determines the size of the table. Using values close to 0 or 1 allows one to track whether "intervened" predictors are pushed toward the extreme of the estimated predictor density, which could indicate lack of support for the scale of the implied intervention (e.g. delta is too big).
#'
#' @export
#' @examples
#' \dontrun{
#' data(metals, package="qgcomp")
#' # subset of predictors
#' XYlist = list(X=metals[,1:4], Y=metals$y)
#' Y_learners = .default_continuous_learners()
#' Xbinary_learners = list(Lrnr_stepwise$new(name="SW"))
#' Xdensity_learners = .default_density_learners(n_bins=c(10))
#' set.seed(1231)
#' vi_ipw <- varimp(X=XYlist$X,Y=XYlist$Y, delta=0.1, Y_learners = Y_learners,
#'        Xdensity_learners=Xdensity_learners[1:2], Xbinary_learners=Xbinary_learners,
#'        estimator="IPW")
#' dx_dens(vi_ipw, Acol=1, delta=0.01)
#' # this shows that most observations keep a similar density range
#' dx_dens(vi_ipw, Acol=1, delta=0.01, quantiles=c(0, 0.1, 0.3, 0.7, 0.9, 1))
#' }
dx_dens <- function(vibr_fit, Acol=1, delta=0.01, quantiles=c(0, 0.1, 0.9, 1)){
  #X = dat[,c(mixturela, "ridageyr")]
  if(!is.null(vibr_fit$qfit)){
    task <- vibr_fit$qfit$training_task
    varnms <- task$nodes$covariates
  }else{
    task <- vibr_fit$gfits[[1]]$training_task
    varnms <- c(task$nodes$outcome, task$nodes$covariates)
  }
  ft <- vibr_fit$gfits[[Acol[1]]]
  if(is.null(delta)) delta = vibr_fit$delta
  xnm = names(X)
  Xc <- .shift(X,Acol,shift = -delta)
  X$set="obs"
  Xc$set="int"
  tsk <- sl3_Task$new(data = X,
                      covariates = xnm[-Acol[1]],
                      outcome = xnm[Acol[1]])
  tskc <- sl3_Task$new(data = Xc,
                       covariates = xnm[-Acol[1]],
                       outcome = xnm[Acol[1]])
  X$dens = ft$predict(tsk)[[1]]
  X$densshift = ft$predict(tskc)[[1]]
  xo=order(X$dens)
  X = X[xo,]
  X$wt = X$densshift/X$dens
  q <- c(-Inf,0,quantile(X$dens, quantiles), Inf)
  denscut <- cut(X$dens, breaks = q)
  densshiftcut <- cut(X$densshift, breaks = q)
  cat(paste0("# obs by quantile-based categories of estimated density: ", xnm[Acol[1]], "\n"))
  print(tab <- table("Density: observed"=as.character(denscut),
              "Density: shifted"=as.character(densshiftcut)
              ))
  invisible(tab)
}


.wtB <- function(X,Acol,gfun,gfits, delta){
  X0 <- .shift(X,Acol, -X[,Acol])
  #OT = sl3::variable_type(type="binary", levels=c(0,max(X[,Acol])))
  g0 <- gfun(X0,Acol,gfits=gfits)
  #function(g0, shift, X, Acol, retcols=1)
  Haw = .Hawb(g0, delta, X, Acol, retcols=1)
  abs(Haw) # shift/p(a|w) = avg proportion of trt probability in the shift
}

.wtC <- function(X,Acol,gfun,gfits, delta){
  # define shifts
  Xa <- .shift(X,Acol, -delta)
  Xb <- .shift(X,Acol,  delta)
  #
  gn <- gfun(X,Acol,gfits=gfits)
  ga <- gfun(Xa,Acol,gfits=gfits)
  gb <- gfun(Xb,Acol,gfits=gfits)
  #
  Haw = .Haw(gn, ga, gb)
  Haw # g(a-delta | w)/g(a | w) =
}

.summarizewts <- function(wts, threshold){
  qs <- as.numeric(quantile(wts, c(0.0, 0.99, 1.0)))
  c(mean=mean(wts), min=qs[1], p99=qs[2], max=qs[3], pctabovethreshold=mean(wts>threshold))
}

.EstWts <- function(X,
                    delta,
                    gfun,
                    gfits,
                    threshold,
                    ...
){
  #return(NULL)# remove when done
  isbin_vec <- apply(X, 2, function(x) length(unique(x))==2)
  resmat <- matrix(NA, nrow=length(isbin_vec), ncol=5)
  for(Acol in seq_len(length(isbin_vec))){
    if(isbin_vec[Acol]){
      wts <- .wtB(X,Acol,gfun,gfits, delta)
    } else{
      wts <- .wtC(X,Acol,gfun,gfits, delta)
    }
    tm <- .summarizewts(wts, threshold)
    resmat[Acol,] <- tm
  }
  colnames(resmat) <- names(tm)
  rownames(resmat) <- names(X)
  resmat <- data.frame(resmat)
  resmat
}


#' Weight-based numerical identifiability diagnostics
#'
#' This function allows one to check, prior to performing inference with \code{vibr::varimp}, whether implied stochastic interventions may be subject to sparsity. Primarily, this approach is based off of estimates of generalized-propensity score weights, where extreme values can suggest highly influential observations due to sparsity.
#'
#' Generally, the identifiability will not be obtained if there are some values of the implied stochastic intervention that have a probability mass/density = 0. This will often not occur in fitted models due to some form of local parametric smoothing, so instead looking for extreme values inverse mass/density based weights can help to suggest where the implied stochastic intervention is extrapolating beyond the observed predictor data.
#'
#' @param X data frame of predictors
#' @param delta change in each column of predictors in call to varimp corresponding to stochastic intervention
#' @param Xdensity_learners list of sl3 learners used to estimate the density of continuous predictors, conditional on all other predictors in X
#' @param Xbinary_learners list of sl3 learners used to estimate the probability mass of continuous predictors, conditional on all other predictors in X
#' @param verbose (logical) print extra information
#' @param ... passed to sl3::make_sl3_Task (e.g. weights)
#'
#' @export
#' @examples
#' \dontrun{
#' data(metals, package="qgcomp")
#' XYlist = list(X=metals[,1:23], Y=metals$y)
#' Xbinary_learners = .default_binary_learners()
#' Xdensity_learners = .default_density_learners(n_bins=c(5, 20))
#' set.seed(12321)
#' # check for intervention = 0.02 standard deviations (scale_continuous=TRUE
#' # will scale continuous predictors to have sd=0.5)
#' ident <- precheck_identification(X=XYlist$X[,1:23], delta=0.01,
#'        Xdensity_learners=Xdensity_learners[c(1,2,3)],
#'        Xbinary_learners=Xbinary_learners, threshold=10,
#'        scale_continuous = TRUE)
#' ident
#' # some extreme weights suggest using a smaller delta. This can be done
#' # by manually scaling variables with extreme weights to have a larger standard deviation
#' # (so that delta would imply a smaller effect size), or one can simply set
#' # delta to a smaller value.
#' }
#'
#'
#'
precheck_identification <- function(
  X,
  delta=0.1,
  Xdensity_learners=NULL,
  Xbinary_learners=NULL,
  verbose=FALSE,
  scale_continuous = TRUE,
  threshold=10,
  ...
){
  # identify the number of observations for which the estimated propensity score is below a threshold
  # in 2x2 grids: identify
  if(scale_continuous){
    if(verbose) cat("Scaling all continuous variables by 2*sd\n")
    # divide continuous by 2*sd
    X = .scale_continuous(X)
  }
  tasklist = .create_tasks(X,NULL,delta, ...)
  if(verbose) cat(paste0("delta = ", delta, "\n")) # TODO: better interpretation
  if(verbose) cat(paste0("Training models\n"))
  sl.gfits <- .train_allX(X, tasklist$slX, Xbinary_learners, Xdensity_learners, verbose=verbose)
  wttable <- .EstWts(X,delta,gfun=.gfunction,gfits=sl.gfits, threshold=threshold,...)
  res <- list(
    gfits = sl.gfits,
    threshold=threshold,
    delta = delta,
    res = wttable
  )
  class(res) <- c("vibr_identify")
  res
}


#' @export
print.vibr_identify <- function(x,...){
  cat("Weight Learners: \n")
  for(idx in seq_len(length(x$gfits))){
    ff = x$gfits[[idx]]
    pr = paste0(row.names(x$res)[idx], ": ", paste0(names(ff$learner_fits$Stack$learner_fits), collapse=", "))
    cat(pr);cat("\n")
  }
  cat(paste0("\nEstimated weights under shift of delta = ", x$delta, " : \n"))
  cat("Binary predictor weights: abs(delta)/Pr(A=1|W=w)\n")
  cat("Continuous predictor weights: g(x-delta|w)/g(a|w)\n")
  print(x$res)
}
