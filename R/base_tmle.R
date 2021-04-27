################################
## efficient influence functions
################################

.hawc1 <- function(gn, ga){
  ga/gn
}


.hawb1 <- function(X,Acol,delta, gn){
  # above this will need to be move dout
  delta*(2*X[,Acol] - 1)/gn + 1
}

.haw1 <- function(X,Acol,delta,gn,ga,isbin=FALSE){
  if(isbin) return(.hawc1(gn, ga))
  if(!isbin) return(.hawb1(X,Acol,delta, gn))
}


.hawc2 <- function(X,Acol,qfit, ...){
  qa <- qfun(qfit=qfit,...)
  eqa <- predict(lm(y~., data.frame(y=qfb, X=X[,-Acol])))   # simple linear regression on W to get E_g[Q | W]
  db2 <- qa - eqa
  db2
}


.hawb2 <- function(X,Acol,qfit, ...){
  qfb = qfun(Xb, Acol,qfit=qfit,...)
  eqfb = predict(lm(y~., data.frame(y=qfb, X=X[,-Acol])))   # simple linear regression on W to get E_g[Q | W]
  dc2 <- qfb - eqfb
  dc2
}



.haw2 <- function(X,Acol,qfit, ...,isbin=FALSE){
  if(isbin) return(.hawc2(X,Acol,qfit, ...))
  if(!isbin) return(.hawb2(X,Acol,qfit, ...))
}


#.haw1(X,Acol,delta,gn,ga,isbin=FALSE)
#.haw2(X,Acol,qfit, ...,isbin=FALSE)
.StepTmle <- function(Y,Qk=qk,Hk=hk,gnk=gnk,gAk=gAk,isbin=FALSE){
  # 1. estimate epsilon
  epsk <- as.numeric(lm(Y~ -1 + offset(Qk) + hk)$coefficients[1])
  # 2. update Qk
  Qk1 <- Qk + epsk*hk
  Xs <- Xsd <- X
  # update gk evaluated at A
  gnk <- exp(epk*hk)*gnk
  # update gk evaluated at A - delta
  gAk <- exp(epk*hk)*gAk
  #update hk
  qfunction(Xsd, Acol)*gfunction(Xs,Acol)
}


.DcTMLE <- function(n,
                X,
                Y,
                Acol,
                delta,
                qfun,
                gfun,
                qfit,
                gfits,
                ...){
  Xb <- Xa <- X
  Xa[,Acol] <- X[,Acol]-delta
  Xb[,Acol] <- X[,Acol]+delta
  qfb = qfun(Xb, Acol,qfit=qfit,...)
  eqfb = predict(lm(y~., data.frame(y=qfb, X=X[,-Acol])))   # simple linear regression on W to get E_g[Q | W]
  dc1 <- gfun(Xa,Acol,gfits=gfits,...)/gfun(X,Acol,gfits=gfits,...)*(Y - qfun(qfit=qfit,...))
  dc2 <- qfun(qfit=qfit,...) - eqfb
  dc3 <- eqfb - Y                # Y doesn't show up in Diaz,vdl 2012
  as.vector(dc1 + dc2 + dc3)
}

.DbTMLE <- function(n,
                X,
                Y,
                Acol,
                delta,
                qfun,
                gfun,
                qfit,
                gfits,
                ...){
  X1 <- X0 <- X
  X1[,Acol] <- 1
  X0[,Acol] <- 0
  qa <- qfun(qfit=qfit,...)
  eqa <- predict(lm(y~., data.frame(y=qfb, X=X[,-Acol])))   # simple linear regression on W to get E_g[Q | W]
  q1 <- qfun(X1,Acol,qfit=qfit,...)
  q0 <- qfun(X0,Acol,qfit=qfit,...)
  db1 <- (delta*(2*X[,Acol] - 1)/gfun(X,Acol,gfits=gfits,...) + 1)*(Y - qa)
  db2 <- qa - eqa
  db3 <- delta*(q1-q0) + eqa - Y # Y doesn't show up in Diaz,vdl 2012
  as.vector(db1 + db2 + db3)
}


################################
#: estimating equations
################################

.MakeTmleEst <- function(dphi){
  #summary(fit <- lm(dphi~1))
  est <- mean(dphi)
  D <- dphi-est
  se = sqrt(mean(D^2)/length(D))
  c(est=est, se = se, z=est/se)
}

.EstEqTMLE <- function(n,
                       X,
                       Y,
                       delta,
                       qfun,
                       gfun,
                       qfit,
                       gfits,
                       ...){
  return(NULL)# remove when done
  isbin_vec <- apply(X, 2, function(x) length(unique(x))==2)
  resmat <- matrix(NA, nrow=length(isbin_vec), ncol=3)
  for(Acol in seq_len(length(isbin_vec))){
    if(isbin_vec[Acol]){
      dphi <- .DbTMLE(n,X,Y,Acol,delta,qfun,gfun,qfit,gfits,...)
    } else{
      dphi <- .DcTMLE( n,X,Y,Acol,delta,qfun,gfun,qfit,gfits,...)
    }
    tm <- .MakeTmleEst(dphi)
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
#' Variable importance using targeted minimum loss estimation
#' @description Not usually called by users
#'
#' @param X data frame of predictors
#' @param Y outcome
#' @param delta change in each column of X corresponding to
#' @param Y_learners list of sl3 learners used to predict the outcome, conditional on all predictors in X
#' @param Xdensity_learners list of sl3 learners used to estimate the density of continuous predictors, conditional on all other predictors in X
#' @param Xbinary_learners list of sl3 learners used to estimate the probability mass of continuous predictors, conditional on all other predictors in X
#' @param verbose (logical) print extra information
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
.varimp_tmle <- function(X,
                         Y,
                         delta=0.1,
                         Y_learners=NULL,
                         Xdensity_learners=NULL,
                         Xbinary_learners=NULL,
                         verbose=TRUE,
                         ...){
  tasklist = .create_tasks(X,Y,delta)
  n = length(Y)
  if(verbose) cat(paste0("Default delta = ", delta, "\n")) # TODO: better interpretation

  isbin <- as.character((length(unique(Y))==2))
  sl.qfit <- .train_Y(X,Y, Y_learners, verbose=verbose, isbin)
  sl.gfits <- .train_allX(X, tasklist$slX, Xbinary_learners, Xdensity_learners, verbose=verbose)
  #.gfunction(X=NULL,Acol,gfits=sl.gfits)
  #.qfunction(X=NULL,Acol,qfit=sl.qfit)
  fittable <- .EstEqTMLE(n,X,Y,delta,qfun=.qfunction,gfun=.gfunction,qfit=sl.qfit,gfits=sl.gfits, ...)
  res <- list(
    res = fittable,
    qfit = sl.qfit,
    gfits = sl.gfits,
    binomial = isbin,
    type = "TMLE"
  )
  class(res) <- c("vibr.fit", class(res))
  res
}
