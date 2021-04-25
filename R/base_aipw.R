################################
## efficient influence functions
################################

.Dc <- function(n,X,Y,Acol,delta,qfun,gfun,qfit,gfits, ...){
  Xb <- Xa <- X
  Xa[,Acol] <- X[,Acol]-delta
  Xb[,Acol] <- X[,Acol]+delta
  p2 <- gfun(Xa,Acol,gfits=gfits,...)/gfun(X,Acol,gfits=gfits,...)
  p3 <- Y - qfun(qfit=qfit,...)
  as.vector(p2*p3) + qfun(Xb, Acol,qfit=qfit,...) - Y  # Y doesn't show up in Diaz,vdl 2012
}

.Db <- function(n,X,Y,Acol,delta,qfun,gfun,qfit,gfits, ...){
  X1 <- X0 <- X
  X1[,Acol] <- 1
  X0[,Acol] <- 0
  qa <- qfun(qfit=qfit,...)
  q1 <- qfun(X1,Acol,qfit=qfit,...)
  q0 <- qfun(X0,Acol,qfit=qfit,...)
  p2 <- delta*(2*X[,Acol] - 1)/gfun(X,Acol,gfits=gfits,...) + 1
  p3 <- Y - qa
  p2*p3 + qa + delta*(q1-q0) - Y  # Y doesn't show up in Diaz,vdl 2012
}


################################
#: estimating equations
################################

.EstEq <- function(dphi){
  #summary(fit <- lm(dphi~1))
  est <- mean(dphi)
  D <- dphi-est
  se = sqrt(mean(D^2)/length(D))
  c(est=est, se = se, z=est/se)
}

.EstEqAIPW <- function(n,X,Y,delta,qfun,gfun,qfit,gfits,...){
  isbin_vec <- apply(X, 2, function(x) length(unique(x))==2)
  resmat <- matrix(NA, nrow=length(isbin_vec), ncol=3)
  for(Acol in seq_len(length(isbin_vec))){
    if(isbin_vec[Acol]){
      #dphi <- .Db(n,X,Y,Acol,delta,qfun,gfun, ...)
      dphi <- .Db(n,X,Y,Acol,delta,qfun,gfun,qfit,gfits)
    } else{
      #dphi <- .Dc( n,X,Y,Acol,delta,qfun,gfun, ...)
      dphi <- .Dc( n,X,Y,Acol,delta,qfun,gfun,qfit,gfits)
    }
    tm <- .EstEq(dphi)
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
# TODO: add "Y learners"
#' Title
#'
#' @param X
#' @param Y
#' @param delta
#' @param Y_learners
#' @param Xdensity_learners
#' @param Xbinary_learners
#' @param verbose
#'
#' @return
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
#' vi <- .varimp_aipw(X=XYlist$X,Y=XYlist$Y, delta=0.1, Y_learners = Y_learners[1:4],Xdensity_learners=Xdensity_learners[1:2], Xbinary_learners=Xbinary_learners[1:2] )
#' vi
#' }
# coef(summary(lm(XYlist$Y~as.matrix(XYlist$X))))[-1,]
# y = rnorm(N, 0, 0.5) + (0.35 + ph*0.1)*calcium*sodium*zinc
# + (0.35 + ph*0.1)*iron*selenium*selenium
# + 1*(calcium>mean(calcium))
# - (0.1 + ph*0.1)*lead*cadmium*arsenic
# - (0.1 + ph*0.1)*chromium*mercury,

.varimp_aipw <- function(X,Y, delta=0.1, Y_learners=NULL, Xdensity_learners=NULL, Xbinary_learners=NULL, verbose=TRUE, ...){
  tasklist = .create_tasks(X,Y,delta)
  n = length(Y)
  if(verbose) cat(paste0("Default delta = ", delta, "\n")) # TODO: better interpretation

  isbin <- as.character((length(unique(Y))==2))
  sl.qfit <- .train_Y(X,Y, Y_learners, verbose=TRUE, isbin)
  sl.gfits <- .train_allX(X, tasklist$slX, Xbinary_learners, Xdensity_learners, verbose=TRUE)
  #.gfunction(X=NULL,Acol,gfits=sl.gfits)
  #.qfunction(X=NULL,Acol,qfit=sl.qfit)
  .EstEqAIPW(n,X,Y,delta,qfun=.qfunction,gfun=.gfunction,qfit=sl.qfit,gfits=sl.gfits)
}
