#' Title
#'
#' @param n sample size
#' @param p number of predictors (with random correlation matrix)
#' @param ncat number of dichotomous predictors (must be <= p)
#'
#' @return
#' @export
#' @importFrom mvtnorm rmvnorm
#'
#' @examples
.dgm <- function(n,p,ncat){
  if(p<ncat) stop("p must be >= ncat")
  # continuous
  mu = rep(0, p)
  Sigma = diag(rep(1,p))
  for(j in 1:p) for(k in 1:(j-1)) Sigma[j,k] <- Sigma[k,j] <- runif(1)
  Zc = mvtnorm::rmvnorm(n, mu)
  # categorize ncat of the continuous predictors
  Zc[,1:ncat] <- 1.0*(Zc[,1:ncat]>0)
  Ey = Zc %*% runif(p,0,3) + (Zc*Zc) %*% runif(p,0,1)
  Y = -2 + Ey + rnorm(n, 0, 1)
  list(Y=Y,X=data.frame(Zc))
}
