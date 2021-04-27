library(mvtnorm)
library(vibr)

options(mc.cores=12)
dgm <- function(n, delta, beta){
  z <- runif(n)
  x <- rbinom(n,1,0.5)
  y <- rnorm(n,0,1) + z*beta[1] + x*beta[2] + x*z*beta[3]
  trdiff <- apply(cbind(beta[1]*delta + beta[3]*delta*z, beta[2]*delta + beta[3]*delta*x), 2, mean)
  list(X = cbind(z,x), y=y, tr = trdiff)
}

density_learners <- function(n_bins=c(10, 20), histtypes=c("equal.mass", "equal.length")){
  density_learners=list()
  idx = 1
  for(nb in n_bins){
    for(histtype in histtypes){
      #nm <- paste0("hist_xgb_", nb, histtype)
      #density_learners[[idx]] <- Lrnr_density_discretize$new(name=nm, categorical_learner = Lrnr_xgboost$new(), n_bins = nb, bin_method=histtype)
      #idx  = idx+1
      nm = paste0("hist_Unadj_", nb, histtype)
      density_learners[[idx]] <- Lrnr_density_discretize$new(name=nm, categorical_learner = Lrnr_mean$new(), n_bins = nb, bin_method=histtype)
      idx  = idx+1
    }
  }
  density_learners
}

continuous_learners <- function(){
  continuous_learners=list(
    Lrnr_glm$new(name="OLS"),
    #Pipeline$new(customize_chain(Lrnr_glm$new(name="OLS"), Lrnr_define_interactions$new(name="INT", list(list(1,2))))),
    Pipeline$new(Lrnr_define_interactions$new(name="INT", list(c("x", "z"), c("z", "z"))), Lrnr_glm$new(name="OLS")),
    Pipeline$new(Lrnr_define_interactions$new(name="INT", list(c("x", "z"), c("z", "z"))), Lrnr_glmnet$new(name="LASSO", alpha=1.0)),
    Lrnr_glmnet$new(name="LASSO", alpha=1.0)
  )
  continuous_learners
}

binary_learners <- function(){
  bin_learners=list(
    Lrnr_glm$new(name="LOGIT")
  )
  bin_learners
}



analyze <- function(i, ...){
  dat = dgm(...)
  #dat = dgm( n=1000, delta = 0.1, beta = c(2,1, .3))
  (vimp <- varimp(data.frame(dat$X),dat$y, delta=0.1, Y_learners=continuous_learners(),
         Xdensity_learners=density_learners(), Xbinary_learners=binary_learners(),
         verbose=FALSE, estimator="AIPW"))
  #vimp$qfit$learner_fits$Stack$learner_fits$`Pipeline(INT->OLS)`$learner_fits$OLS$coefficients
  #vimp$qfit$learner_fits$Stack$learner_fits$OLS$coefficients
  #vimp$qfit$learner_fits$Stack$learner_fits$`Pipeline(INT->LASSO)`$learner_fits$LASSO$coefficients
  #vimp$qfit$learner_fits$Stack$learner_fits$OLS$fit_object$coefficients
  obj <- as.matrix(vimp$res)
  tr = dat$tr
  names(tr) <- colnames(dat$X)
  lmfit <- summary(lm(y~., data.frame(y=dat$y, dat$X/0.1)))
  # vimp$qfit$learner_fits$Stack$learner_fits$OLS$fit_object$R
  # chol(solve(lmfit$cov)) # close to cholesky decomposition of the Hessian, but varying in signs
  c(
    drest = obj[1:2,1],
    drse = obj[1:2,2],
    lmest = lmfit$coefficients[2:3, 1],
    lmse = lmfit$coefficients[2:3, 2],
    tr = tr
    )
}

analyze(1, n=100, delta = 0.1, beta = c(2,1, .3))

future::availableCores()
future::plan("multicore")
res = future.apply::future_lapply(1:1000, analyze, n=100, delta = 0.1, beta = c(2,1), future.seed=TRUE)
res = as.data.frame(do.call(rbind, res))
truthx = matrix(NA, nrow=nrow(res), ncol = ncol(res))
truthx[,1] <- truthx[,5] <- 0.2
truthx[,2] <- truthx[,6] <- 0.1
print(apply(res[,c(1:2, 5:6)], 2, function(x) c(mean=mean(x), sd=sd(x))))
print(apply((res-truthx)[,c(1:2, 5:6)], 2, function(x) c(bias=mean(x), mse=mean(x^2), sd.bias=sd(x))))
print(apply(res[,c(3:4, 7:8)], 2, function(x) c(mean=mean(x))))



#lnr  = Pipeline$new(Lrnr_define_interactions$new(name="INT", list(c(1, 2))), Lrnr_glm$new(name="OLS"))
#XX <- sl3_Task$new(
#  data=data.frame(dat$X, y=dat$y),
#  outcome="y",
#  covariates=names(data.frame(dat$X))
#)
#
#res = lnr$train(XX)
#res$fit_object$learner_fits$OLS$coefficients


