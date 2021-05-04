.wtB <- function(X,Acol,gfun,gfits, delta){
  gn <- gfun(X,Acol,gfits=gfits)
  #
  Haw = .Hawb(gn, delta, X, Acol)
  Haw
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
  Haw
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


#' Variable importance in cross-sectional data
#'
#' @param X data frame of predictors
#' @param delta change in each column of X corresponding to
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
#' ident <- VI_identify(X=XYlist$X, delta=0.1, Y_learners = Y_learners,
#'        Xdensity_learners=Xdensity_learners[1:2], Xbinary_learners=Xbinary_learners,
#'        estimator="TMLE", verbose=TRUE)
#' ident
#' }
VI_identify <- function(X,
                        delta=0.1,
                        Xdensity_learners=NULL,
                        Xbinary_learners=NULL,
                        verbose=FALSE,
                        ...){
  tasklist = .create_tasks(X,NULL,delta)
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
