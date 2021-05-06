library(mvtnorm)
library(vibr)
library(txshift)
library(sl3)

simspline <- function(k, x, degree=1){
  basis.mat <- matrix(NA,nrow=length(x), ncol=(1+length(k))*degree)
  colidx = 1
  for(knot in 0:length(k)){
    for(deg in 1:(degree)){
      if(knot==0) b <- x
      if(knot >0) b <- pmax(0, x-k[knot])
      basis.mat[,colidx] <- b^deg
      colidx = colidx + 1
    }
  }
  basis.mat
}
#simspline(c(-1, 0, 1), rnorm(100), degree=1)


options(mc.cores=12)
dgm <- function(n, delta, beta, degree=1, zk = c(-1, 0, 1)){
  fullbeta <- rep(0, 20)
  fullbeta[1:length(beta)] <- beta
  px <- 0.5
  z <- rnorm(n)
  x <- rbinom(n,1,px)
  zint <- rnorm(n, delta)
  xint <- rbinom(n,1,px+delta)
  yfun <- function(z,x,beta){
    zbasis <- simspline(zk, z, degree)
    pz <- ncol(zbasis)
    x*beta[1] + x*z*beta[2] + zbasis %*% beta[3:(3+pz-1)]
  }
  yerr <- rnorm(n,0,1)
  y <- yerr + yfun(z,x,fullbeta)
  yintz <- yerr + yfun(zint,x,fullbeta)
  yintx <- yerr + yfun(z,xint,fullbeta)
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
    Lrnr_glm$new(name="OLS", family=gaussian()),
    Pipeline$new(Lrnr_define_interactions$new(name="INT", list(c("x", "z"), c("z", "z"))), Lrnr_glmnet$new(name="LASSO", alpha=1.0)),
    Lrnr_gam$new(name="GAM")
    #Pipeline$new(Lrnr_define_interactions$new(name="INT", list(c("x", "z"), c("z", "z"))), Lrnr_glm$new(name="OLS"))
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
  dat = dgm( n=500, delta = 0.1, beta = c(1,2, -1.0, 1.0, -1.0), degree=1, zk = c(-1.5, 0, 1.5))
  #plot(dat$X[,1], dat$y)
  #plot(dat$X[,2], dat$y)
  ####
  .dl <- continuous_learners()
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
  (res1 <- f1(density_learners(), .tskz))
  (res2 <- f1(continuous_learners(), .tsk))
  sl_learner = res1[[1]]
  #sl_learner$train(.tskz)

  set.seed(123123)

  sl_learner_density <- Lrnr_sl$new(learners = Stack$new(density_learners()), metalearner = Lrnr_solnp_density_quiet$new())
  sl_learner          <-Lrnr_sl$new(learners = Stack$new(continuous_learners()[1:2]), metalearner = Lrnr_nnls$new(convex=TRUE))
  sl_learner_bin      <-Lrnr_sl$new(learners = Stack$new(binary_learners()), metalearner = Lrnr_nnls$new(convex=TRUE))
  #sl_learner_density$train_sublearners(.tskz)
  #sl_learner_density$train(.tskz)
  #sl_learner_bin$train_sublearners(.tskx)
  #sl_learner_bin$train(.tskx)
  #sl_learner$train_sublearners(.tsk)
  #sl_learner$train(.tsk)

  (tmle <- txshift(W = data.frame(dat$X[,"x", drop=FALSE]), A = dat$X[,"z"], Y = as.numeric(dat$y), delta = .1,
          estimator="tmle",
          g_exp_fit_args = list(fit_type = "sl",
                                sl_learners_density = sl_learner_density
          ),
          Q_fit_args = list(fit_type = "sl",
                              sl_learners = sl_learner
          )
  ))
  (tmle2 <- txshift(W = dat$X[,"z"], A = dat$X[,"x"], Y = as.numeric(dat$y), delta =.1,
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
  (vimp <- varimp(data.frame(dat$X),dat$y, delta=.1, Y_learners=continuous_learners()[1:3],
                  Xdensity_learners=density_learners(), Xbinary_learners=binary_learners(),
                  verbose=FALSE, estimator="TMLE", estimand="diff"))
  (vimp2 <- varimp(data.frame(dat$X),dat$y, delta=.1, Y_learners=continuous_learners()[1:3],
                  Xdensity_learners=density_learners()[[1]], Xbinary_learners=binary_learners(),
                  verbose=FALSE, estimator="TMLE", estimand="diff", B=5))
  (vimp <- varimp(data.frame(dat$X),dat$y, delta=.1, Y_learners=continuous_learners()[1:3],
                  Xdensity_learners=density_learners(), Xbinary_learners=binary_learners(),
                  verbose=FALSE, estimator="AIPW", estimand="diff"))
  (vimp2 <- varimp(data.frame(dat$X),dat$y, delta=.1, Y_learners=continuous_learners()[1:3],
                  Xdensity_learners=density_learners(), Xbinary_learners=binary_learners(),
                  verbose=FALSE, estimator="AIPW", estimand="diff", B=5))
  (vimp <- varimp(data.frame(dat$X),dat$y, delta=.1, Y_learners=continuous_learners()[1:3],
                  Xdensity_learners=density_learners(), Xbinary_learners=binary_learners(),
                  verbose=FALSE, estimator="GCOMP", estimand="diff"))
  (vimp2 <- varimp(data.frame(dat$X),dat$y, delta=.1, Y_learners=continuous_learners()[1:3],
                  Xdensity_learners=density_learners(), Xbinary_learners=binary_learners(),
                  verbose=FALSE, estimator="GCOMP", estimand="diff", B=5))
  (vimp <- varimp(data.frame(dat$X),dat$y, delta=.1, Y_learners=continuous_learners()[1:3],
                  Xdensity_learners=density_learners(), Xbinary_learners=binary_learners(),
                  verbose=FALSE, estimator="IPW", estimand="diff"))
  (vimp2 <- varimp(data.frame(dat$X),dat$y, delta=.1, Y_learners=continuous_learners()[1:3],
                  Xdensity_learners=density_learners(), Xbinary_learners=binary_learners(),
                  verbose=FALSE, estimator="IPW", estimand="diff", B=5))

  print(tmle)
  print(vimp)
  dat$tr
  summary(lm(y ~ I((.5*z/sd(z))/.1) + I(x/.1), data = data.frame(y=dat$y, dat$X)))$coefficients

}


stabilitytest <- function(...){
  data(metals, package="qgcomp")
  dat = list(X=metals[,1:10], y=metals$y) # up t 23

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
      Lrnr_glm$new(name="OLS", family=gaussian()),
      #Pipeline$new(Lrnr_define_interactions$new(name="INT", list(c("x", "z"), c("z", "z"))), Lrnr_glmnet$new(name="LASSO", alpha=1.0)),
      Lrnr_gam$new(name="GAM")
      #Pipeline$new(Lrnr_define_interactions$new(name="INT", list(c("x", "z"), c("z", "z"))), Lrnr_glm$new(name="OLS"))
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

  set.seed(NULL)
  #set.seed(1231)
  (vimp <- varimp(data.frame(dat$X),dat$y, delta=.1, Y_learners=continuous_learners()[1:3],
                  Xdensity_learners=density_learners(), Xbinary_learners=binary_learners(),
                  verbose=FALSE, estimator="TMLE", estimand="diff", folds=20))


  #set.seed(1231)
  (vimp2 <- varimp(data.frame(dat$X[,1:3]),dat$y, delta=.1, Y_learners=continuous_learners()[1:3],
                  Xdensity_learners=density_learners()[1:2], Xbinary_learners=binary_learners()[1:2],
                  verbose=FALSE, estimator="TMLE", estimand="diff", folds=5))


  set.seed(1231)
  t2 <- sl3::sl3_Task$new(data = dat$X, folds = 20, covariates = names(dat$X))
  set.seed(1231)
  t1 <- sl3::sl3_Task$new(data = dat$X, folds = 20, covariates = names(dat$X))
  t1$folds[[10]]$validation_set
  t2$folds[[10]]$validation_set

  #outcome
  plot(vimp$qfit$predict(), vimp2$qfit$predict()) # stable under a seed
  plot(vimp$qfit$learner_fits$Stack$predict(), vimp2$qfit$learner_fits$Stack$predict()) # stable under a seed
  # density
  plot(vimp$gfits[[1]]$predict()[[1]], vimp2$gfits[[1]]$predict()[[1]])
  plot(vimp$gfits[[2]]$predict()[[1]], vimp2$gfits[[2]]$predict()[[1]])
  plot(vimp$gfits[[6]]$predict()[[1]], vimp2$gfits[[6]]$predict()[[1]])
  plot(vimp$gfits[[6]]$learner_fits$Stack$predict()[[1]], vimp2$gfits[[6]]$learner_fits$Stack$predict()[[1]])
  plot(vimp$gfits[[6]]$learner_fits$Stack$predict()[[2]], vimp2$gfits[[6]]$learner_fits$Stack$predict()[[2]])
  plot(vimp$gfits[[6]]$learner_fits$Stack$predict()[[3]], vimp2$gfits[[6]]$learner_fits$Stack$predict()[[3]])
  plot(vimp$gfits[[6]]$learner_fits$Stack$predict()[[4]], vimp2$gfits[[6]]$learner_fits$Stack$predict()[[4]])
  cbind(vimp$gfits[[6]]$learner_fits$Stack$predict()[[4]], vimp2$gfits[[6]]$learner_fits$Stack$predict()[[4]])
  plot(vimp$res$est, vimp2$res$est)
  plot(vimp$res$se, vimp2$res$se)
  plot(vimp$res$z, vimp2$res$z)
  plot(order(abs(vimp$res$est)), order(abs(vimp2$res$est)))
  plot(vimp$res$est-vimp2$res$est, vimp$res$se-vimp2$res$se)
}

analyze <- function(i, B=0, ...){
  dat = dgm(...)
  #set.seed(12312); dat = dgm( n=1000, delta = 0.1, beta = c(2,1, .3))
  (vimp <- varimp(data.frame(dat$X),dat$y, delta=0.1, Y_learners=continuous_learners(),
                  Xdensity_learners=density_learners(), Xbinary_learners=binary_learners(),
                  verbose=FALSE, estimator="AIPW"))
  (vimp2 <- varimp(data.frame(dat$X),dat$y, delta=0.1, Y_learners=continuous_learners(),
                  Xdensity_learners=density_learners(), Xbinary_learners=binary_learners(),
                  verbose=FALSE, estimator="TMLE"))
  #vimp$qfit$learner_fits$Stack$learner_fits$`Pipeline(INT->OLS)`$learner_fits$OLS$coefficients
  #vimp$qfit$learner_fits$Stack$learner_fits$OLS$coefficients
  #vimp$qfit$learner_fits$Stack$learner_fits$`Pipeline(INT->LASSO)`$learner_fits$LASSO$coefficients
  #vimp$qfit$learner_fits$Stack$learner_fits$OLS$fit_object$coefficients
  obj <- as.matrix(vimp$res)
  obj2 <- as.matrix(vimp$res)
  tr = dat$tr
  names(tr) <- colnames(dat$X)
  lmfit <- summary(lm(y~., data.frame(y=dat$y, dat$X/0.1)))
  # vimp$qfit$learner_fits$Stack$learner_fits$OLS$fit_object$R
  # chol(solve(lmfit$cov)) # close to cholesky decomposition of the Hessian, but varying in signs
  c(
    AIPWest = obj[1:2,1],
    AIPWse = obj[1:2,2],
    TMLEest = obj2[1:2,1],
    TMLEse = obj2[1:2,2],#apply(vimp2$boots,2,sd),
    lmest = lmfit$coefficients[2:3, 1],
    lmse = lmfit$coefficients[2:3, 2],
    tr = tr
    #TMLEseasymp = obj2[1:2,2]
  )
}


#dgm(n=1000000, delta = 0.1, beta = c(2,1, .0))$tr
t(res1 <- analyze(1231321, n=50, delta = 0.1, beta = c(1,2, -1.0, 1.0, -1.0, -1.0), degree=1, zk = c(-1.5, 0, 1.5)))

future::availableCores()
future::plan("multicore")
resL = future.apply::future_lapply(1:1000, analyze, n=50, delta = 0.1, beta = c(1,2, -1.0, 1.0, -1.0, -1.0), degree=1, zk = c(-1.5, 0, 1.5), future.seed=TRUE)
res = as.data.frame(do.call(rbind, resL))
rm <- apply(res, 2, mean)

cipow <- function(res, root="AIPW", exp="x", type="cover"){
  nm = paste0(root, type, ".", exp)
  est = res[,paste0(root, "est.", exp)]
  se = res[,paste0(root, "se.", exp)]
  tr = rm[paste0("tr.", exp)] # truth based on average across simulations
  res[,nm] <<- switch(type,
                   cover= as.numeric(((est + 1.96*se) > tr) & ((est - 1.96*se) < tr)),
                   power= as.numeric(((est - 1.96*se) > 0) | ((est + 1.96*se) < 0)),
                   bias= as.numeric(est-tr)
  )
}
cipow(res, "AIPW", "z", "cover")
cipow(res, "AIPW", "x", "cover")
cipow(res, "TMLE", "z", "cover")
cipow(res, "TMLE", "x", "cover")
cipow(res, "lm", "z", "cover")
cipow(res, "lm", "x", "cover")
#
cipow(res, "AIPW", "x", "power")
cipow(res, "AIPW", "z", "power")
cipow(res, "TMLE", "x", "power")
cipow(res, "TMLE", "z", "power")
cipow(res, "lm", "z", "power")
cipow(res, "lm", "x", "power")
#
cipow(res, "AIPW", "z", "bias")
cipow(res, "AIPW", "x", "bias")
cipow(res, "TMLE", "z", "bias")
cipow(res, "TMLE", "x", "bias")
cipow(res, "lm", "z", "bias")
cipow(res, "lm", "x", "bias")

print(apply(res[,  grep("est.", names(res))], 2, function(x) c(mean=mean(x), sd=sd(x))))
print(apply((res)[,grep("bias", names(res))], 2, function(x) c(bias=mean(x), mse=mean(x^2), sd.bias=sd(x))))
print(apply(res[,grep("se[a]*", names(res))], 2, function(x) c(mean=mean(x))))
print(apply(res[, grep("cover", names(res))], 2, function(x) c(mean=mean(x))))



#lnr  = Pipeline$new(Lrnr_define_interactions$new(name="INT", list(c(1, 2))), Lrnr_glm$new(name="OLS"))
#XX <- sl3_Task$new(
#  data=data.frame(dat$X, y=dat$y),
#  outcome="y",
#  covariates=names(data.frame(dat$X))
#)
#
#res = lnr$train(XX)
#res$fit_object$learner_fits$OLS$coefficients


