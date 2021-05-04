################################
## efficient influence functions
################################

.DcAIPW <- function(n,
                X,
                Y,
                Acol,
                delta,
                qfun,
                gfun,
                qfit,
                gfits,
                estimand,
                ...){
  # define shifts
  Xa <- .shift(X,Acol, -delta)
  Xb <- .shift(X,Acol,  delta)
  Hk <- gfun(Xa,Acol,gfits=gfits,...)/gfun(X,Acol,gfits=gfits,...)
  qfb = qfun(Xb, Acol,qfit=qfit,...)
  #eqfb = predict(lm(y~., data.frame(y=qfb, X=X[,-Acol])))   # simple linear regression on W to get E_g[Q | W]
  eqfb <- 0 # cancels out
  dc1 <- Hk*(Y - qfun(qfit=qfit,...))
  dc2 <- qfb - eqfb
  dc3 <- eqfb - Y*(estimand != "mean")                # Y doesn't show up in Diaz,vdl 2012 b/c they are estimating mean Y|A+delta
  as.vector(dc1 + dc2 + dc3)
}


.DbAIPW <- function(n,
                X,
                Y,
                Acol,
                delta,
                qfun,
                gfun,
                qfit,
                gfits,
                estimand,
                ...){
  X1 <- X0 <- X
  X1[,Acol] <- 1
  X0[,Acol] <- 0
  #Hk <- gfun(Xa,Acol,gfits=gfits,...)/gfun(X,Acol,gfits=gfits,...) # unsure why this one is not used
  Hk <- (delta*(2*X[,Acol] - 1)/gfun(X,Acol,gfits=gfits,...) + 1)
  qa <- qfun(qfit=qfit,...)
  #eqa <- predict(lm(y~., data.frame(y=qfb, X=X[,-Acol])))   # simple linear regression on W to get E_g[Q | W]
  eqa <- 0 # cancels out
  q1 <- qfun(X1,Acol,qfit=qfit,...)
  q0 <- qfun(X0,Acol,qfit=qfit,...)
  db1 <- Hk*(Y - qa)
  db2 <- qa - eqa
  db3 <- delta*(q1-q0) + eqa - Y*(estimand != "mean") # Y doesn't show up in Diaz,vdl 2012 b/c they are estimating mean Y|A+delta
  as.vector(db1 + db2 + db3)
}



################################
#: estimating equations
################################

.MakeiAipwEst <- function(dphi, n){
  #summary(fit <- lm(dphi~1))
  est <- mean(dphi)
  D <- dphi-est
  avar = mean(D^2) # asymptotic variance
  se = sqrt(avar)/sqrt(n)
  c(est=est, se = se, z=est/se)
}

.EstEqAIPW <- function(n,
                       X,
                       Y,
                       delta,
                       qfun,
                       gfun,
                       qfit,
                       gfits,
                       estimand,
                       ...){
  isbin_vec <- apply(X, 2, function(x) length(unique(x))==2)
  resmat <- matrix(NA, nrow=length(isbin_vec), ncol=3)
  for(Acol in seq_len(length(isbin_vec))){
    if(isbin_vec[Acol]){
      dphi <- .DbAIPW(n,X,Y,Acol,delta,qfun,gfun,qfit,gfits,estimand,...)
    } else{
      dphi <- .DcAIPW( n,X,Y,Acol,delta,qfun,gfun,qfit,gfits,estimand,...)
    }
    tm <- .MakeiAipwEst(dphi, n)
    resmat[Acol,] <- tm
  }
  colnames(resmat) <- names(tm)
  rownames(resmat) <- names(X)
  resmat <- data.frame(resmat)
  resmat$p <- pnorm(-abs(resmat$z))*2
  resmat
}


################################
# expert wrappers
################################
#' Variable importance using estimating equations
#' @description Not usually called by users
#'
#' @param X data frame of predictors
#' @param Y outcome
#' @param delta change in each column of X corresponding to
#' @param Y_learners list of sl3 learners used to predict the outcome, conditional on all predictors in X
#' @param Xdensity_learners list of sl3 learners used to estimand the density of continuous predictors, conditional on all other predictors in X
#' @param Xbinary_learners list of sl3 learners used to estimand the probability mass of continuous predictors, conditional on all other predictors in X
#' @param verbose (logical) print extra information
#' @param estimand (character) "diff" (default, estimate mean difference comparing Y under intervention with observed Y), "mean" (estimate mean Y under intervention)
#' @param ... passed to sl3::base_predict (rare)
#'
#' @return vi object
#' @export
#'
#' @examples
#' \dontrun{
#' XYlist = .dgm(n=100,p=4,ncat=3)
#' data(metals, package="qgcomp")
#' XYlist = list(X=metals[,1:23], Y=metals$y)
#' Y_learners = .default_continuous_learners()
#' Xbinary_learners = .default_binary_learners()
#' Xdensity_learners = .default_density_learners(n_bins=c(5, 20))
#' vi <- .varimp_aipw(X=XYlist$X,Y=XYlist$Y, delta=0.1, Y_learners = Y_learners[1:4],
#' Xdensity_learners=Xdensity_learners[1:2], Xbinary_learners=Xbinary_learners[1:2] )
#' vi
#' }
.varimp_aipw <- function(X,
                         Y,
                         delta=0.1,
                         Y_learners=NULL,
                         Xdensity_learners=NULL,
                         Xbinary_learners=NULL,
                         verbose=TRUE,
                         estimand,
                         ...){
  tasklist = .create_tasks(X,Y,delta)
  n = length(Y)
  if(verbose) cat(paste0("Default delta = ", delta, "\n"))

  isbin <- as.character((length(unique(Y))==2))
  sl.qfit <- .train_Y(X, Y, Y_learners, verbose=verbose, isbin)
  sl.gfits <- .train_allX(X, tasklist$slX, Xbinary_learners, Xdensity_learners, verbose=verbose)
  fittable <- .EstEqAIPW(n,X,Y,delta,qfun=.qfunction,gfun=.gfunction,qfit=sl.qfit,gfits=sl.gfits, estimand, ...)
  res <- list(
    res = fittable,
    qfit = sl.qfit,
    gfits = sl.gfits,
    binomial = isbin,
    type = "AIPW"
  )
  class(res) <- c("vibr.fit", class(res))
  res
}
