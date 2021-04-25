################################
## efficient influence functions - based on 2014 paper (same representation in binary and continuous exposures)
################################

.Dcw <- function(n,X,Y,Acol,delta,qfun,gfun,qfit=NULL,gfits, ...){
  Xb <- Xa <- X
  Xa[,Acol] <- X[,Acol]-delta
  #Xb[,Acol] <- X[,Acol]+delta
  p2 <- gfun(Xa,Acol,gfits=gfits,...)/gfun(X,Acol,gfits=gfits,...)
  p3 <- Y
  as.vector(p2*p3)
}

.Dbw <- .Dcw



################################
#: estimating equations
################################


.EstEqIPW <- function(n,X,Y,delta,qfun=NULL,gfun,qfit=NULL,gfits,...){
  isbin_vec <- apply(X, 2, function(x) length(unique(x))==2)
  resmat <- matrix(NA, nrow=length(isbin_vec), ncol=3)
  for(Acol in seq_len(length(isbin_vec))){
    if(isbin_vec[Acol]){
      dphi <- .Dbw(n,X,Y,Acol,delta,qfun=NULL,gfun,qfit=NULL,gfits)
    } else{
      dphi <- .Dcw( n,X,Y,Acol,delta,qfun=NULL,gfun,qfit=NULL,gfits)
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


# \dontrun{
# XYlist = .dgm(n=100,p=4,ncat=3)
# data(metals, package="qgcomp")
# XYlist = list(X=metals[,1:23], Y=metals$y)
# Y_learners = .default_continuous_learners()
# Xbinary_learners = .default_binary_learners()
# Xdensity_learners = .default_density_learners(n_bins=c(5, 20))
# vi3 <- .varimp_ipw(X=XYlist$X,Y=XYlist$Y, delta=0.1, Y_learners = Y_learners[1:4],Xdensity_learners=Xdensity_learners[1:2], Xbinary_learners=Xbinary_learners[1:2] )
# vi3
# }
# coef(summary(lm(XYlist$Y~as.matrix(XYlist$X))))[-1,]
# y = rnorm(N, 0, 0.5) + (0.35 + ph*0.1)*calcium*sodium*zinc
# + (0.35 + ph*0.1)*iron*selenium*selenium
# + 1*(calcium>mean(calcium))
# - (0.1 + ph*0.1)*lead*cadmium*arsenic
# - (0.1 + ph*0.1)*chromium*mercury,

.varimp_ipw <- function(X,Y, delta=0.1, Y_learners=NULL, Xdensity_learners=NULL, Xbinary_learners=NULL, verbose=TRUE, ...){
  tasklist = .create_tasks(X,Y,delta)
  n = length(Y)
  if(verbose) cat(paste0("Default delta = ", delta, "\n")) # TODO: better interpretation

  isbin <- as.character((length(unique(Y))==2))
  #sl.qfit <- .train_Y(X,Y, Y_learners, verbose=TRUE, isbin)
  sl.gfits <- .train_allX(X, tasklist$slX, Xbinary_learners, Xdensity_learners, verbose=TRUE)
  #.gfunction(X=NULL,Acol,gfits=sl.gfits)
  #.qfunction(X=NULL,Acol,qfit=sl.qfit)
  .EstEqIPW(n,X,Y,delta,qfun=NULL,gfun=.gfunction,qfit=NULL,gfits=sl.gfits)
}
