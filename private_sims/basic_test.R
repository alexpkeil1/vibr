library(mvtnorm)
library(vibr)
library(txshift)
library(sl3)

options(mc.cores=12)
dgm <- function(n, delta, beta){
  px <- 0.5
  z <- runif(n, min=0, max=1)
  x <- rbinom(n,1,px)
  zint <- runif(n, min=0+delta, max=1+delta)
  xint <- rbinom(n,1,px+delta)
  yfun <- function(z,x,beta){
    z*beta[1] + x*beta[2] + x*z*beta[3]
  }
  yerr <- rnorm(n,0,1)
  y <- yerr + yfun(z,x,beta)
  yintz <- yerr + yfun(zint,x,beta)
  yintx <- yerr + yfun(z,xint,beta)
  (trdiff <- c(mean(yintz-y), mean(yintx-y)))
  list(X = cbind(z,x), y=y, tr = trdiff)
}

density_learners <- function(n_bins=c(3,8), histtypes=c("equal.mass")){
  #sl3::sl3_list_learners("density") # see all available
  density_learners=list()
  idx = 1
  nm <- paste0("dens_sp_gam_")
  density_learners[[idx]] <- Lrnr_density_semiparametric$new(name=nm, mean_learner = Lrnr_gam$new(), var_learner = NULL)
  idx  = idx+1
  nm <- paste0("dens_sp_mean_")
  density_learners[[idx]] <- Lrnr_density_semiparametric$new(name=nm, mean_learner = Lrnr_mean$new(), var_learner = NULL)
  idx  = idx+1
  for(nb in n_bins){
    for(histtype in histtypes){
      nm = paste0("hist_Unadj_", nb, histtype)
      density_learners[[idx]] <- Lrnr_density_discretize$new(name=nm, categorical_learner = Lrnr_mean$new(), n_bins = nb, bin_method=histtype)
      idx  = idx+1
    }
  }
  density_learners
}

continuous_learners <- function(){
  #sl3::sl3_list_learners("continuous")
  continuous_learners=list(
    Lrnr_mean$new(name="Mean"),
    #Pipeline$new(Lrnr_define_interactions$new(name="INT", list(c("x", "z"), c("z", "z"))), Lrnr_glmnet$new(name="LASSO", alpha=1.0)),
    #Lrnr_glmnet$new(name="LASSO", alpha=1.0),
    #Pipeline$new(Lrnr_define_interactions$new(name="INT", list(c("x", "z"), c("z", "z"))), Lrnr_glm$new(name="OLS"),
    Lrnr_glm$new(name="OLS", family=gaussian())
  )
  continuous_learners
}

binary_learners <- function(){
  #sl3::sl3_list_learners("binomial")
  bin_learners=list(
    Lrnr_glm$new(name="LOGIT", family=binomial()),
    Lrnr_mean$new(name="Mean")
  )
  bin_learners
}



testtxshift <- function(){
  set.seed(12312)
  dat = dgm( n=100, delta = 0.1, beta = c(2,1, .3))
  X = dat$X
  A = X[,"z", drop=TRUE]
  W = X[,"x", drop=FALSE]
  ####
  .dl <- continuous_learners()[1:2]
  .tsk <- sl3_Task$new(data = data.frame(dat$X, y=dat$y),
                       covariates=c("x", "z"),
                       outcome="y")
  .tskx <- sl3_Task$new(data = data.frame(dat$X),
                        covariates=c("z"),
                        outcome="x")
  .tskz <- sl3_Task$new(data = data.frame(dat$X),
                        covariates=c("x"),
                        outcome="z")

  f1 <- function(dl=.dl,tsk=.tsk){
    learner_stack <- Stack$new(dl)
    learner_fit <- learner_stack$train(tsk)
    #cross validation
    cv_stack <- Lrnr_cv$new(learner_stack)
    cv_fit <- cv_stack$train(tsk)
    cv_task <- cv_fit$chain()
    # train super learner
    type = "cont"
    metalearner <- Lrnr_nnls$new(convex=TRUE)
    sl_fit <- metalearner$train(cv_task)
    #sl_pipeline <- make_learner(Pipeline, learner_fit, sl_fit)
    sl_pipeline <- Pipeline$new(learner_fit, sl_fit)
    list(sl_pipeline, metalearner)
  }
  res1 <- f1(density_learners(), .tskz)
  res2 <- f1(continuous_learners(), .tsk)
  sl_learner = res1[[1]]
  #sl_learner$train(.tskz)
  A = X[,"z", drop=TRUE]
  W = as.matrix((X[,"x", drop=FALSE]))

  set.seed(123123)

  sl_learner_density <- Lrnr_sl$new(learners = Stack$new(density_learners()), metalearner = Lrnr_solnp_density_quiet$new())
  sl_learner          <-Lrnr_sl$new(learners = Stack$new(continuous_learners()[1:2]), metalearner = Lrnr_nnls$new(convex=TRUE))
  sl_learner_bin      <-Lrnr_sl$new(learners = Stack$new(binary_learners()), metalearner = Lrnr_nnls$new(convex=TRUE))

  (tmle <- txshift(W = data.frame(dat$X[,"x", drop=FALSE]), A = dat$X[,"z"], Y = (dat$y), delta = .1,
          estimator="tmle",
          g_exp_fit_args = list(fit_type = "sl",
                                sl_learners_density = sl_learner_density
          ),
          Q_fit_args = list(fit_type = "sl",
                              sl_learners = sl_learner
          )
  ))
  (tmle2 <- txshift(W = dat$X[,"z"], A = dat$X[,"x"], Y = (dat$y), delta =.1,
                   estimator="tmle",
                   g_exp_fit_args = list(fit_type = "sl",
                                         sl_learners_density = sl_learner_bin
                   ),
                   Q_fit_args = list(fit_type = "sl",
                                     sl_learners = sl_learner
                   )
  ))
  # delta method to get mean difference comparing Y^A+delta to Y^A
  # if x-theta -> n^-1/2 N(0,sig^2)
  # then g(x)-g(theta) -> n^-1/2 N(0,sig^2 * g'(theta)^2)
  # g(x) \equiv X - Y
  # g(theta) \equiv theta - Y
  # g'(theta) = ?
  dphi <- tmle$eif + tmle$psi - dat$y
  mean(dphi)
  dphi2 <- tmle2$eif + tmle2$psi - dat$y
  mean(dphi2)
  sd(dphi)/sqrt(length(dat$y))
  sd(dphi2)/sqrt(length(dat$y))

  set.seed(123123)
  tmle
  tmle2
  # bounded estimator in progress, will look like tmle package
  (vimp <- varimp(data.frame(dat$X),dat$y, delta=.1, Y_learners=continuous_learners()[1:2],
                  Xdensity_learners=density_learners(), Xbinary_learners=binary_learners(),
                  verbose=FALSE, estimator="TMLE", estimand="diff"))
  (vimp2 <- varimp(data.frame(dat$X),dat$y, delta=.1, Y_learners=continuous_learners()[1:2],
                  Xdensity_learners=density_learners()[[1]], Xbinary_learners=binary_learners(),
                  verbose=TRUE, estimator="TMLE", estimand="diff", B=5))
  (vimp <- varimp(data.frame(dat$X),dat$y, delta=.1, Y_learners=continuous_learners()[1:2],
                  Xdensity_learners=density_learners(), Xbinary_learners=binary_learners(),
                  verbose=FALSE, estimator="AIPW", estimand="diff"))
  (vimp2 <- varimp(data.frame(dat$X),dat$y, delta=.1, Y_learners=continuous_learners()[1:2],
                  Xdensity_learners=density_learners(), Xbinary_learners=binary_learners(),
                  verbose=FALSE, estimator="AIPW", estimand="diff", B=5))
  (vimp <- varimp(data.frame(dat$X),dat$y, delta=.1, Y_learners=continuous_learners()[1:2],
                  Xdensity_learners=density_learners(), Xbinary_learners=binary_learners(),
                  verbose=FALSE, estimator="GCOMP", estimand="diff"))
  (vimp2 <- varimp(data.frame(dat$X),dat$y, delta=.1, Y_learners=continuous_learners()[1:2],
                  Xdensity_learners=density_learners(), Xbinary_learners=binary_learners(),
                  verbose=FALSE, estimator="GCOMP", estimand="diff", B=5))
  (vimp <- varimp(data.frame(dat$X),dat$y, delta=.1, Y_learners=continuous_learners()[1:2],
                  Xdensity_learners=density_learners(), Xbinary_learners=binary_learners(),
                  verbose=FALSE, estimator="IPW", estimand="diff"))
  (vimp2 <- varimp(data.frame(dat$X),dat$y, delta=.1, Y_learners=continuous_learners()[1:2],
                  Xdensity_learners=density_learners(), Xbinary_learners=binary_learners(),
                  verbose=FALSE, estimator="IPW", estimand="diff", B=5))

  print(tmle)
  print(vimp)
  summary(lm(y ~ I(z/.1) + I(x/.1), data = data.frame(y=dat$y, dat$X)))

}



analyze <- function(i, ...){
  dat = dgm(...)
  #set.seed(12312); dat = dgm( n=1000, delta = 0.1, beta = c(2,1, .3))
  (vimp <- varimp(data.frame(dat$X),dat$y, delta=0.1, Y_learners=continuous_learners(),
                  Xdensity_learners=density_learners(), Xbinary_learners=binary_learners(),
                  verbose=FALSE, estimator="AIPW"))
  (vimp2 <- varimp(data.frame(dat$X),dat$y, delta=0.1, Y_learners=continuous_learners(),
                  Xdensity_learners=density_learners(), Xbinary_learners=binary_learners(),
                  verbose=FALSE, estimator="TMLE", B=200))
  #vimp$qfit$learner_fits$Stack$learner_fits$`Pipeline(INT->OLS)`$learner_fits$OLS$coefficients
  #vimp$qfit$learner_fits$Stack$learner_fits$OLS$coefficients
  #vimp$qfit$learner_fits$Stack$learner_fits$`Pipeline(INT->LASSO)`$learner_fits$LASSO$coefficients
  #vimp$qfit$learner_fits$Stack$learner_fits$OLS$fit_object$coefficients
  obj <- as.matrix(vimp$res)
  obj2 <- as.matrix(vimp2$est$res)
  tr = dat$tr
  names(tr) <- colnames(dat$X)
  lmfit <- summary(lm(y~., data.frame(y=dat$y, dat$X/0.1)))
  # vimp$qfit$learner_fits$Stack$learner_fits$OLS$fit_object$R
  # chol(solve(lmfit$cov)) # close to cholesky decomposition of the Hessian, but varying in signs
  c(
    AIPWest = obj[1:2,1],
    AIPWse = obj[1:2,2],
    TMLEest = obj2[1:2,1],
    TMLEse = apply(vimp2$boots,2,sd),
    lmest = lmfit$coefficients[2:3, 1],
    lmse = lmfit$coefficients[2:3, 2],
    tr = tr,
    TMLEseasymp = obj2[1:2,2]
  )
}


#dgm(n=1000000, delta = 0.1, beta = c(2,1, .0))$tr
t(res1 <- analyze(1231321, n=100, delta = 0.1, beta = c(2,1, -.25)))

future::availableCores()
future::plan("multicore")
resL = future.apply::future_lapply(1:100, analyze, n=100, delta = 0.1, beta = c(2,1, -.25), future.seed=TRUE)
res = as.data.frame(do.call(rbind, resL))
rm <- apply(res, 2, mean)

cipow <- function(res, root="AIPW", exp="x", type="cover"){
  nm = paste0(root, type, ".", exp)
  est = res[,paste0(root, "est.", exp)]
  se = res[,paste0(root, "se.", exp)]
  tr = rm[paste0("tr.", exp)]
  res[,nm] <<- switch(type,
                   cover= as.numeric(((est + 1.96*se) > tr) & ((est - 1.96*se) < tr)),
                   power= as.numeric(((est - 1.96*se) > 0) | ((est + 1.96*se) < 0))
 )
}
cipow(res, "TMLE", "x", "cover")
cipow(res, "TMLE", "z", "cover")
cipow(res, "AIPW", "x", "cover")
cipow(res, "AIPW", "z", "cover")
cipow(res, "lm", "x", "cover")
cipow(res, "lm", "z", "cover")
cipow(res, "TMLE", "x", "power")
cipow(res, "TMLE", "z", "power")
cipow(res, "AIPW", "x", "power")
cipow(res, "AIPW", "z", "power")
cipow(res, "lm", "x", "power")
cipow(res, "lm", "z", "power")

truthx = matrix(NA, nrow=nrow(res), ncol = ncol(res))
truthx[,c(1,5,9)] <- 0.2
truthx[,c(2,6,10)] <- 0.1
print(apply(res[,c(1:2, 5:6)], 2, function(x) c(mean=mean(x), sd=sd(x))))
print(apply((res-truthx)[,c(1:2, 5:6, 9:10)], 2, function(x) c(bias=mean(x), mse=mean(x^2), sd.bias=sd(x))))
print(apply(res[,c(3:4, 7:8, 11:12)], 2, function(x) c(mean=mean(x))))
print(apply(res[,c(15:25)], 2, function(x) c(mean=mean(x))))



#lnr  = Pipeline$new(Lrnr_define_interactions$new(name="INT", list(c(1, 2))), Lrnr_glm$new(name="OLS"))
#XX <- sl3_Task$new(
#  data=data.frame(dat$X, y=dat$y),
#  outcome="y",
#  covariates=names(data.frame(dat$X))
#)
#
#res = lnr$train(XX)
#res$fit_object$learner_fits$OLS$coefficients


