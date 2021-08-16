
#' In progress: visualizing a stochastic intervention
#'
#' @param vibr.fit a fit from varimp
#' @param X predictors from a varimp fit
#' @param Acol (integer) which column of X to diagnose
#' @param delta (numeric, default=0.01) change in each column of X corresponding to stochastic intervention
#'
#' @return ggplot2 plot object
#' @export
#' @import ggplot2
#'
#' @examples
#' set.seed(123)
plotshift_dens <- function(vibr.fit, X, Acol=1, delta){
  if(is.null(delta)) delta <- vibr.fit$delta
  if(1==1){
    task <- vibr.fit$qfit$training_task
    varnms <- task$nodes$covariates
  }
  if(1==2){
    task <- vibr.fit$gfits[[1]]$training_task
    varnms <- c(task$nodes$outcome, task$nodes$covariates)
  }
  dat <- task$data
  X <- data.frame(dat)[,varnms,drop=FALSE]
  ft <- vibr.fit$gfits[[Acol[1]]]
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

#' In progress: diagnosing weight problems
#'
#' @param vibr.fit a fit from varimp
#' @param X predictors from a varimp fit
#' @param Acol (integer) which column of X to diagnose
#' @param delta (numeric, default=0.01) change in each column of X corresponding to stochastic intervention
#'
#' @export
#' @import ggplot2
#'
#' @examples
#' set.seed(123)
plotshift_wt <- function(vibr.fit, X, Acol=1, delta=0.01){
  if(is.null(delta)) delta <- vibr.fit$delta
  if(1==1){
    task <- vibr.fit$qfit$training_task
    varnms <- task$nodes$covariates
  }
  if(1==2){
    task <- vibr.fit$gfits[[1]]$training_task
    varnms <- c(task$nodes$outcome, task$nodes$covariates)
  }
  dat <- task$data
  ft <- vibr.fit$gfits[[Acol[1]]]
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

#' In progress: diagnosing fit problems
#'
#'
#' @param vibr.fit a fit from varimp
#' @param X predictors from a varimp fit
#' @param Acol (integer) which column of X to diagnose
#' @param Bcol (integer) second column of X to diagnose
#' @param delta (numeric, default=0.01) change in each column of X corresponding to stochastic intervention
#'
#' @export
#' @import ggplot2
#'
#' @examples
#' set.seed(123)
plotshift_scatter <- function(vibr.fit, Acol, Bcol, delta=NULL, joint=FALSE){
  if(is.null(delta)) delta <- vibr.fit$delta
  if(1==1){
    task <- vibr.fit$qfit$training_task
    varnms <- task$nodes$covariates
  }
  if(1==2){
    task <- vibr.fit$gfits[[1]]$training_task
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



#' In progress: diagnosing problems in estimating densities
#'
#' @param vibr.fit a fit from varimp
#' @param X predictors from a varimp fit
#' @param Acol (integer) which column of X to diagnose
#' @param delta (numeric, default=0.01) change in each column of X corresponding to stochastic intervention
#'
#' @export
#' @examples
#' set.seed(123)
dx_dens <- function(vibr.fit, X, Acol=1, delta=0.01){
  #X = dat[,c(mixturela, "ridageyr")]
  if(vibr.fit$scaled){
    X = .scale_continuous(X)
  }
  ft <- vibr.fit$gfits[[Acol[1]]]
  if(is.null(delta)) delta = vibr.fit$delta
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
  q <- c(-Inf,0,quantile(X$dens, c(0, 0.1, 0.9, 1)), Inf)
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

.summarizewts <- function(wts){
  qs <- as.numeric(quantile(wts, c(0.0, 1.0)))
  c(mean=mean(wts), min=qs[1], max=qs[2])
}

.EstWts <- function(X,
                    delta,
                    gfun,
                    gfits,
                    ...
){
  #return(NULL)# remove when done
  isbin_vec <- apply(X, 2, function(x) length(unique(x))==2)
  resmat <- matrix(NA, nrow=length(isbin_vec), ncol=3)
  for(Acol in seq_len(length(isbin_vec))){
    if(isbin_vec[Acol]){
      wts <- .wtB(X,Acol,gfun,gfits, delta)
    } else{
      wts <- .wtC(X,Acol,gfun,gfits, delta)
    }
    tm <- .summarizewts(wts)
    resmat[Acol,] <- tm
  }
  colnames(resmat) <- names(tm)
  rownames(resmat) <- names(X)
  resmat <- data.frame(resmat)
  resmat
}


#' In progress: check for identifiability of stochastic interventions
#'
#' @param X data frame of predictors
#' @param delta change in each column of X corresponding to stochastic intervention
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
#' ident <- check_identification(X=XYlist$X, delta=0.1, Y_learners = Y_learners,
#'        Xdensity_learners=Xdensity_learners[1:2], Xbinary_learners=Xbinary_learners,
#'        estimator="TMLE", verbose=TRUE)
#' ident
#' }
#'
#'
#'
check_identification <- function(
  X,
  delta=0.1,
  Xdensity_learners=NULL,
  Xbinary_learners=NULL,
  verbose=FALSE,
  scale_continuous = TRUE,
  ...
){
  # identify the number of observations for which the estimated propensity score is zero
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
  wttable <- .EstWts(X,delta,gfun=.gfunction,gfits=sl.gfits,...)
  res <- list(
    res = wttable
  )
  class(res) <- c("vibr.identify", class(res))
  res
}



