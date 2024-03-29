# sl3 like functions -----------------------------------------------------------

.call_with_args_vibr <- function (fun, args, other_valid = list(), keep_all = FALSE,
                                  silent = FALSE, ignore = c())
{
  args <- args[!(names(args) %in% ignore)]
  if (!keep_all) {
    formal_args <- names(formals(fun))
    all_valid <- c(formal_args, other_valid)
    invalid <- names(args)[which(!(names(args) %in% all_valid))]
    args <- args[which(names(args) %in% all_valid)]
    if (!silent & length(invalid) > 0) {
      message(sprintf("Learner called function %s with unknown args: %s. These will be dropped.\nCheck the params supported by this learner.",
                      as.character(substitute(fun)), paste(invalid,
                                                           collapse = ", ")))
    }
  }
  do.call(fun, args)
}


# convenience functions -------------------------------------------------------------

.logit <- function(p){log(p) - log1p(-p)}
.expit <- function(mu){1/(1+exp(-mu))}
.identity <- function(x){x}
.invidentity <- .identity

.shift <- function(X,Acol,shift){
  X[,Acol] = X[,Acol] + shift
  X
}


# Data handlers -------------------------------------------------------------

#' @importFrom sl3 sl3_Task variable_type
.create_tasks <- function(X,Y, V=NULL, delta=0.1, whichcols=NULL, ...){
  X <- as.data.frame(X)
  p <- ncol(X)
  if(is.null(whichcols)) whichcols = seq_len(p)
  nm <- names(X)
  slX <- list()
#  slXpdelta <- list() # X + delta
#  slXmdelta <- list() # X-delta
  if(!is.null(V[1])){
    X = data.frame(cbind(X, V)) # V can hold weights, offsets, etc.
  }
  for(j in whichcols){
    isbin = length(unique(X[,j]))==2
    tt = ifelse(isbin, "binomial", "continuous")
    slX[[j]] <- sl3_Task$new(outcome_type=variable_type(type=tt), data=X,
                             covariates=nm[-j],
                             outcome=nm[j], ...
    )
  }
  # TODO: give this a class
  list(
    slX = slX
  )
}

.scale_continuous <- function(X, isbin_vec=NULL){
  p <- ncol(X)
  if(is.null(isbin_vec[1])){
    if(p==1) isbin_vec = c(length(unique(X))==2)
    if(p>1) isbin_vec <- apply(X, 2, function(x) length(unique(x))==2)
  }
  for(j in seq_len(p)){
    if(!isbin_vec[j]){
      X[,j] <- X[,j]/(2*sd(X[,j]))
    }
  }
  X
}

#' @importFrom MASS ginv
.stdres <- function(resids, df, Z, train=FALSE){
  sse <- sum(resids^2)
  sig <- sqrt(sse/df) # sqrt of "usual" estimate of error variance if df = n-p
  if(train){
    # yields standard error of residuals, reflecting the fact that residuals at ends of domain
    # will be more influential
    Z <- as.matrix(Z)
    itzz = MASS::ginv(t(Z) %*% Z)
    # leverage
    Hi <-  apply(Z, 1, function(zi) ((zi %*% itzz) %*% zi)[1])
    ret = sqrt(1-Hi)*sig # standard error of the residual
  } else{
    ret = sig # for new data
  }
  ret
}


.bootsample <- function(n, strata=NULL, id=NULL){
  if(is.null(strata)){
    ridx <- sample(seq_len(n), n, replace=TRUE)
  }
  if(!is.null(strata)){
    idx <- seq_len(n)
    lidx <- tapply(idx,strata, function(x) sample(x, length(x), replace=TRUE))
    ridx <- as.numeric(do.call(c, lidx))
  }
  ridx
}


# default learners -------------------------------------------------------------



#' @import sl3
#' @export
.default_density_learners_big <- function(n_bins=c(7, 12), histtypes=c("equal.mass", "equal.length")){
  density_learners=list(
     Lrnr_density_gaussian$new(name="dens_gaussian_glm"),
     Lrnr_density_discretize$new(name="hist_multinom_10", categorical_learner = Lrnr_multinom$new(trace=FALSE), n_bins = 10, bin_method="equal.mass"),
     Lrnr_density_gaussian$new(name="dens_gaussian_glmlogy", transfun = log),
     Lrnr_density_hse$new(name="dens_hse_glm", mean_learner = Lrnr_glm$new(), var_learner = Lrnr_glm$new()),
     Lrnr_density_hse$new(name="dens_hse_lasso", mean_learner = Lrnr_glmnet$new(name="Lasso", alpha=1.0, family="gaussian"), var_learner = Lrnr_glmnet$new(name="Lasso", alpha=1.0, family="gaussian")),
     Lrnr_density_hse$new(name="dens_hse_sw", mean_learner = Lrnr_stepwise$new(name="Stepwise"), var_learner = Lrnr_stepwise$new(name="Stepwise")),
     Lrnr_density_hse$new(name="dens_hse_polymars", mean_learner = Lrnr_polspline$new(name="polymars", family=gaussian()), var_learner = Lrnr_glmnet$new(name="Lasso", alpha=1.0, family="gaussian")),
     Lrnr_density_hse$new(name="dens_hse_nnet", mean_learner = Lrnr_nnet$new(name="nnet",maxit=200, trace=FALSE), var_learner = Lrnr_glmnet$new(name="Lasso", alpha=1.0, family="gaussian")),
     Pipeline$new(Lrnr_pca$new(name="pca"), Lrnr_density_gaussian$new(name="dens_gaussian_glm")),
     Pipeline$new(Lrnr_pca$new(name="pca"),Lrnr_density_gaussian$new(name="dens_gaussian_glmlogtrans", transfun=log)),
     Lrnr_density_semiparametric$new(name="dens_sp_unadj", mean_learner = Lrnr_mean$new(), var_learner = NULL),
     Lrnr_density_semiparametric$new(name="dens_sp_glm", mean_learner = Lrnr_glm$new(), var_learner = NULL)
  )
  idx  = length(density_learners)+1
  for(nb in n_bins){
    for(histtype in histtypes){
      nm <- paste0("hist_rf_", nb, histtype)
      density_learners[[idx]] <- Lrnr_density_discretize$new(name=nm, categorical_learner = Lrnr_randomForest$new(), n_bins = nb, bin_method=histtype)
      idx  = idx+1
      nm = paste0("hist_unadj_", nb, histtype)
      density_learners[[idx]] <- Lrnr_density_discretize$new(name=nm, categorical_learner = Lrnr_mean$new(), n_bins = nb, bin_method=histtype)
      idx  = idx+1
    }
  }
  density_learners
}


#' @export
.default_density_learners <- function(...){
  vibr::.default_density_learners_big()[1:4]
}


#' @import sl3
#' @export
.default_continuous_learners_big <- function(){
  continuous_learners=list(
    Lrnr_mean$new(name="mean"),
    Lrnr_glm$new(name="ols", family=gaussian()),
    Lrnr_gam$new(name="gam"),
    Lrnr_polspline_quiet$new(name="polymars"),
    Lrnr_glmnet$new(name="cv_ridge", alpha=0.0, family="gaussian"), # formerly cv_elastic_net
    Lrnr_stepwise$new(name="stepwise", family=gaussian()),
    #Lrnr_xgboost$new(),
    Lrnr_nnet$new(name="nnet",maxit=200, trace=FALSE),
    Lrnr_glmnet$new(name="cv_lasso", alpha=1.0, family="gaussian"),
    Pipeline$new(Lrnr_pca$new(name="pca"), Lrnr_glm$new(name="ols", family=gaussian())), # PCA plus glm
    #Pipeline$new(Lrnr_screener_coefs$new(name="lassocreen", learner=Lrnr_glmnet$new(name="lasso", alpha=1.0, family="gaussian")), Lrnr_earth$new(name="MARS", family=gaussian())), # screen by coefficient size then OLS
    Pipeline$new(Lrnr_screener_coefs$new(name="coefscreen", learner=Lrnr_glm$new(name="ols", family=gaussian())), Lrnr_glm$new(name="OLS", family=gaussian())), # screen by coefficient size then OLS
    Pipeline$new(Lrnr_screener_importance$new(name="rfimpscreen", learner=Lrnr_randomForest$new()), Lrnr_glm$new(name="OLS", family=gaussian())), # screen by variable importance then OLS
    Lrnr_glmnet$new(name="cv_elastic_net", alpha=0.5, family="gaussian")
  )
  continuous_learners
}

#' @export
.default_continuous_learners <- function(){
  .default_continuous_learners_big()[1:4]
}

#' @import sl3
#' @export
.default_binary_learners_big <- function(){
  bin_learners=list(
    Lrnr_mean$new(name="mean"),
    Lrnr_glm$new(name="logit", family=binomial()),
    Lrnr_polspline_quiet$new(name="polymars", family=binomial()),
    Lrnr_stepwise$new(name="stepwise", family=binomial()),
    Lrnr_glmnet$new(name="cv_elastic_net", alpha=0.0, family="binomial"),
    #Lrnr_xgboost$new(),
    Lrnr_nnet$new(name="nnet",maxit=200, trace=FALSE),
    Lrnr_glmnet$new(name="cv_lasso", alpha=1.0, family="binomial"),
    Pipeline$new(Lrnr_pca$new(name="pca"), Lrnr_glm$new(name="logit", family=binomial())), # PCA plus glm
    Pipeline$new(Lrnr_screener_coefs$new(name="lassocreen", learner=Lrnr_glmnet$new(name="lasso", alpha=1.0, family="binomial")), Lrnr_polspline$new(name="polymars", family=binomial())), # screen by coefficient size then OLS
    Pipeline$new(Lrnr_screener_coefs$new(name="coefscreen", learner=Lrnr_glm$new(name="logit", family=binomial())), Lrnr_glm$new(name="logit", family=binomial())), # screen by coefficient size then OLS
    Pipeline$new(Lrnr_screener_importance$new(name="rfimpscreen", learner=Lrnr_randomForest$new()), Lrnr_glm$new(name="logit", family=binomial())) # screen by variable importance then LOGIT
  )
  bin_learners
}

#' @export
.default_binary_learners <- function(){
  .default_binary_learners_big()[1:4]
}

#' @import sl3
#' @export
.default_metalearner <- function(type){
  switch(substr(type[1],1,1),
                        d = Lrnr_solnp_density_quiet$new(trace=0),
                        h = Lrnr_solnp_density_quiet$new(trace=0),
                        e = Lrnr_nnls$new(convex=TRUE),
                        p = Lrnr_nnls$new(convex=TRUE)
  )
}



# training functions -------------------------------------------------------------





## Train super learner
##
## @param datatask sl3 task
## @param learners R list of learners
## @param type type of fit ("density", "probability", or "expectation")
##
## @description Train a super learner fit, given some data held in an sl3 task object and a list of learners
##
## @return a trained sl3 fit: defaults to super learner ensemble if a list of multiple learners are given (output from sl3::make_learner)
#' @export
#' @import sl3
.train_fullstack <- function(datatask, learners, metalearner=NULL, type=c("density", "probability", "expectation")){
  #sl3_Task$new()
  # train learners
  if(length(learners)>1){
    learner_stack <- Stack$new(learners)
    learner_fit <- learner_stack$train(datatask)
    #cross validation
    #cv_stack <- Lrnr_cv$new(learner_stack, full_fit = TRUE) # seems like it should be needed, but makes no difference
    cv_stack <- Lrnr_cv$new(learner_stack)
    cv_fit <- cv_stack$train(datatask)
    cv_task <- cv_fit$chain()
    # train super learner
    if(is.null(metalearner)) metalearner <- .default_metalearner(type)
    sl_fit <- metalearner$train(cv_task)
    sl_pipeline <- make_learner(Pipeline, learner_fit, sl_fit)
    ret <- sl_pipeline
  } else{
    ret <- learners[[1]]$train(datatask)
  }
  ret
}

.train_cvsuperlearner <- function(datatask, learners, metalearner=NULL, type=c("density", "probability", "expectation")){
  #sl3_Task$new()
  # train learners
  learner_stack <- Stack$new(learners) # alt is do.call(Stack.new, learners)
  learner_fit <- learner_stack$train(datatask)
  #cross validation
  cv_stack <- Lrnr_cv$new(learner_stack, full_fit=TRUE)
  cv_fit <- cv_stack$train(datatask)
  # train super learner
  if(is.null(metalearner)) metalearner <- .default_metalearner(type)
  cv_meta_task <- cv_fit$chain(datatask)
  #cv_meta_fit <- cv_meta_task$train()
  sl_fit <- metalearner$train(cv_meta_task) # trained sublearners

  full_stack_fit <- cv_fit$fit_object$full_fit
  #sl_pipeline <- make_learner(Pipeline, learner_fit, sl)fit

  sl_pipeline <- make_learner(Pipeline, full_stack_fit, sl_fit)
  sl_pipeline$train(datatask)
  sl_pipeline
}

.train_cvsuperlearner_delayed <- function(datatask, learners, metalearner=NULL, type=c("density", "probability", "expectation")){
  if(is.null(metalearner)) metalearner <- .default_metalearner(type)
  sl <- Lrnr_sl$new(learners=learners,
                    metalearner = metalearner,
  )
  obj <- delayed_learner_train(sl, datatask)
  cvsl_pipeline <- obj$compute()
  cvsl_pipeline
}




#' @import sl3
.train_Y <- function(X,Y, learners, metalearner=NULL, verbose=TRUE, isbin=FALSE, V=NULL, ...){
  df = data.frame(X, Y)
  nms = names(df)
  pp1 <- ncol(df)
  if(!is.null(V[1])) df = data.frame(cbind(df, V))
  XY <- sl3_Task$new(data=df,
                     outcome_type=variable_type(type=ifelse(isbin, "binomial", "continuous")),
                     covariates=nms[-pp1],
                     outcome=nms[pp1], ...
  )
  if(verbose) cat(paste0("Training: ", names(df)[pp1], "(", ifelse(isbin, "binomial", "continuous, predicted"), ")\n"))
  .train_fullstack(XY, learners, metalearner=NULL, type="expectation")
}


#' @import sl3
#' @importFrom future future value
#' @export
.train_allX <- function(X, tasks, bin_learners, density_learners, metalearner=NULL, verbose=TRUE, V=NULL, whichcols=NULL){
  # TODO: check for data frame X
  p = ncol(X)
  if(is.null(whichcols)) whichcols = seq_len(p)
  nvars = length(whichcols)
  #if(p==1) vartype = ifelse(length(unique(X[,1]))==2, "b", "c")
  #if(p>1) vartype = ifelse(apply(X, 2, function(x) length(unique(x))==2), "b", "c")
  if(nvars==1) vartype = ifelse(length(unique(X[,whichcols]))==2, "b", "c")
  if(nvars>1) vartype = ifelse(apply(X[,whichcols], 2, function(x) length(unique(x))==2), "b", "c")

  #trained_models <- list()
  nmseq <- names(X)
  ee <- new.env()
  #for(j in 1:p){
  for(j in whichcols){
    if(verbose) cat(paste0("Training: ", names(X)[j], "(", ifelse(vartype[j]=="b", "binomial", "continuous, density"), ")\n"))
    ee[[nmseq[j]]] <- future::future( {
      #trained_models[[j]] <-
      switch(vartype[j],
             b = .train_fullstack( tasks[[j]],bin_learners, metalearner=metalearner[[bin_learners]], type="prob"),
             c = .train_fullstack( tasks[[j]],density_learners, metalearner=metalearner[[density_learners]], type="density")
      )
    }, seed=TRUE, lazy=FALSE)
}
  trained_models <-  as.list(future::value(ee))
  trained_models <- trained_models[names(X)]
  names(trained_models) <- NULL
  trained_models

}

# prediction functions -------------------------------------------------------------

## covariate density/probability prediction ----
.gfunction_sl <- function(X=NULL,Acol=1,gfits=NULL, ...){
    if(!is.null(X)){
      XX <- sl3_Task$new(
        data=data.frame(X),
        covariates=names(X)[-Acol],
        outcome=names(X)[Acol],
        ...
      )
      pred <- gfits[[Acol]]$predict(XX)
  }
  if(is.null(X)){
    pred <- gfits[[Acol]]$predict()
  }
  if(typeof(pred)=="list") pred <- pred[[1]]
  pred
}

.gfunction <- function(X=NULL,Acol=1,gfits=NULL, ...){
  #.gfunction_glm(X,Acol,gfits, ...)
  #.gfunction_glm_density(X,Acol,gfits, ...)
  .gfunction_sl(X,Acol,gfits, ...)
}



## outcome prediction ----
.qfunction_sl <- function(X=NULL,Acol=1,qfit=NULL, ...){
    if(!is.null(X)){
    XX <- sl3_Task$new(
      data=data.frame(X),
      covariates=names(X),
      outcome_type = qfit$training_task$outcome_type,
      ...
    )
    pred <- qfit$predict(XX)
  }
  if(is.null(X))  pred <- qfit$predict()
  if(typeof(pred)=="list") pred <- pred[[1]]
  pred
}


.qfunction <- function(X=NULL,Acol=1,qfit=NULL, ...){
  #.qfunction_glm(X,1,qfit, ...)
  .qfunction_sl(X,1,qfit, ...)
}


# shared functions of esimators -------------------------------------------------------------

## clever covariate/weights  ----

.Haw <- function(gn, ga, gb){
  # if intervention would push exposure out of the support of A | W, then don't intervene (potential = observed outcome)
  # g(a-delta) = expected density after probabilistically shifting all individuals to a-delta (continuous)
  # i.e. if intervetion is g(a+delta), then observed individual represents g(a-delta) individuals in intervened world
  # I(A<u(w))g(a-delta|w)/g(a|w) + I(a>=u(w)-delta)
  Haw <- ifelse(gn>0, ga/gn, 0) + as.numeric(gb == 0)
  if(any(is.na(Haw))){
    stop(paste0("",sum(is.na(Haw))/length(Haw), "% of weights (continuous predictor) had missing values\n"))
  }
  if(any(Haw<0)) warning(paste0("",sum(Haw<0), " weights (continuous predictor) were < 0\n"))
  Haw
}

# clever covariate, binary (based on diaz and vdl 2018, translated to shift in propensity score)
.Hawb <- function(g0, shift, X, Acol, retcols=1){
  # if intervention would push exposure out of the support of A | W, then don't intervene (potential = observed outcome)
  # g(a-delta) = expected density after probabilistically shifting all individuals to a-delta (continuous)
  # i.e. if intervetion is g(a+delta), then observed individual represents 1/g(a-delta) individuals in intervened world
  # NOTE:
  #    g(a-shift) = I(A=1)*(g(1)-shift) + I(A=0)*(g(0)+shift)          since support is only on [0,1](binary)
  #    g(a-shift)/g(a)  =  (I(A=1)*(g(1)-shift) + I(A=0)*(g(0)+shift)) / (I(A=1)*g(1) + I(A=0)*g(0))
  #        =  (I(A=1)*(g(1)-shift) + I(A=0)*(g(0)+shift)) / (I(A=1)*g(1) + I(A=0)*g(0))
  #        =  (I(A=1)*(g(1)-shift)/ g(1) + I(A=0)*(g(0)+shift)/g(0)
  #        =  (I(A=1)*(1-shift/ g(1)) + I(A=0)*(1+shift/g(0))
  #        =  1 + (I(A=1)*(-shift/ g(n)) + I(A=0)*(+shift/g(n))
  #        =  1 + shift(-I(A=1)/g(n) + I(A=0)/g(n))
  #        =  1 + shift((I(A=0)-I(A=1))/g(n))
  #
  #
  #   1 + shift*(2A-1)/gn  = 1 + shift(I(A=1)/g1 - I(A=0)/g0)
  #    shift(I(A=1)/g1) = I(A=1)(g1+shift)/g1) -I(A=1)
  #   -I(A=1) + I(A=1)(g1+shift)/g1 -  (-I(A=0)  I(A=0)(g0- shift)/G0) = 1 + shift/G1 - shift/G0
  #   1 + I(A=1)(gn+shift)/gn + I(A=0)(gn- shift)/gn = 1 + I(A=1)(gn+shift)/gn - I(A=0)(gn + shift)/gn
  #   (2A-1)*(gn+shift)/gn = (2A-1)*(1+shift/gn)
  g1 <- 1-g0
  Haw1 <- ifelse(g1>0, 1 + shift/g1, 0) + as.numeric(g1 + shift > 1)
  Haw0 <- ifelse(g0>0, 1 - shift/g0, 0) + as.numeric(g0 - shift < 0)
  Haw <-  X[,Acol]*Haw1 + (1-X[,Acol])*Haw0
  if(any(is.na(Haw))){
    stop(paste0("",sum(is.na(Haw))/length(Haw), "% of weights (binary predictor) had missing values\n"))
  }
  if(any(Haw<0)) warning(paste0("",sum(Haw<0), " weights (binary predictor) were < 0\n"))
  cbind(Haw, Haw1, Haw0)[,1:retcols]
}


# prelim functions -------------------------------------------------------------


.prelims <- function(X, Y, V=NULL, whichcols=seq_len(ncol(X)), delta, Y_learners, Xbinary_learners, Xdensity_learners, verbose=verbose, isbin=NULL, ...){
  #.prelims(X, Y, delta, Y_learners, Xbinary_learners, Xdensity_learners, verbose=verbose, ...)
  # basic error checking
  stopifnot(is.data.frame(X))
  if(!is.null(V[1])) stopifnot(is.data.frame(V))

  doX = !(is.null(Xbinary_learners) & is.null(Xdensity_learners))
  doY = !(is.null(Y_learners))
  if(is.null(isbin)) isbin <- as.logical((length(unique(Y))==2))
  n = length(Y)
  args = list(...)
  # checking for weights
  nmargs = names(args)
  if("weights" %in% nmargs){
    wt <- V[,args$weights, drop=TRUE]
  } else{
    wt <- rep(1.0, n)
  }
  if(!is.logical(all.equal(sum(wt),as.double(n)))){
    if(verbose) cat("Normalizing weights to sum to length(Y)\n")
    wt <- wt/mean(wt)
  }
  tasklist <- NULL
  #if(doX) tasklist = .create_tasks(X[,whichcols, drop=FALSE],Y,V,delta, ...)
  if(doX) tasklist = .create_tasks(X,Y,V,delta,, whichcols=whichcols, ...)
  if(verbose) cat(paste0("delta = ", delta, "\n"))
  yb = .bound_zero_one(Y)
  Ybound = yb[[1]]
  sl.qfit <- sl.gfits <- NULL
  if(doY){
    #yb = .bound_zero_one(Y)
    sl.qfit <- .train_Y(X,Y, Y_learners, verbose=verbose, isbin=isbin, V=V, ...)
  }
  if(doX){
    sl.gfits <- .train_allX(X[,whichcols, drop=FALSE], tasklist$slX, Xbinary_learners, Xdensity_learners, verbose=verbose, V=V)
  }
  list(n=n, Ybound=Ybound, tasklist = tasklist, sl.qfit = sl.qfit, sl.gfits=sl.gfits, isbin = isbin, weights=wt, whichcols=whichcols)
}


.ChExstractFit <- function(obj){
  yfit <- xfit <- tasklist <- FALSE
  if(!is.null(obj$qfit)){
    yfit <- obj$qfit$is_trained
  }
  if(!is.null(obj$gfits)){
    xfit <- obj$gfits[[1]]$is_trained
    if(xfit){
      tasklist <- list()
      tasklist$slX <- lapply(1:length(obj$gfits), function(x) obj$gfits[[x]]$training_task)
    }
  }
  list(yfit=yfit,
       xfit=xfit,
       sl.qfit=obj$qfit,
       sl.gfits=obj$gfits,
       tasklist=tasklist,
       isbin=obj$binomial,
       oldtype=obj$type,
       scaled =obj$scaled,
       weights=obj$weights,
       n = length(obj$weights),
       rank = obj$rank,
       whichcols=obj$whichcols
       )
}


.attach_misc <- function(obj, scale_continuous, delta, B=NULL, whichcols=seq_len(length(obj$gfits))){
  obj$scaled = scale_continuous
  obj$delta = delta
  obj$whichcols = whichcols
  if(is.null(B)){
    obj$rank = rank(-abs(obj$res$est))
    obj$rankz = rank(-abs(obj$res$z))
  }
  if(!is.null(B)){
    obj$rank = rank(-abs(obj$est$res$est))
    obj$sds = apply(obj$boots,2,sd)
    obj$rankz = rank(-abs(obj$est$res$est/obj$sds))
  }
  obj
}

# cross-fitting functions -----------------------------------------------------

.checkeven <- function(val){
  !as.logical(val %% 2)
}

.permn <- function(n){
  if(n==1){
    return(matrix(1))
  } else {
    sp <- .permn(n-1)
    p <- nrow(sp)
    A <- matrix(nrow=n*p,ncol=n)
    for(i in 1:n){
      A[(i-1)*p+1:p,] <- cbind(i,sp+(sp>=i))
    }
    return(A)
  }
}

## double cross fitting ----

.make_xfit_folds <- function(fold, set1, set2, set3){
  fold <- list(fold = fold,
               set1 = set1,
               set2 = set2,
               set3 = set3
  )
  fold
}

.xfitfolds_from_foldvec <- function (r, folds, ordermat)
{
  nfolds <- length(unique(folds))
  remfolds = (nfolds-1)/2
  s1 <- which(folds %in% ordermat[1:remfolds,r])
  s2 <- which(folds %in% ordermat[(remfolds+1):(2*remfolds),r])
  if(.checkeven(r)){
    set1 = s1
    set2 = s2
  }
  if(!.checkeven(r)){
    set1 = s2
    set2 = s1
  }
  set3 <- which(folds == ordermat[nfolds,r])
  .make_xfit_folds(r, set1, set2, set3)
}

.xfitsplit <- function(r=1,n, V=5){
  if(V>1) {
    folds <- rep(seq_len(V), length = n)
    folds <- sample(folds)
    vm1 = max(1,V-1)
    combinations <- combn(V,vm1)
    combinations <- apply(combinations, 2, function(x) x[order(runif(vm1))])
    combinations <- rbind(combinations, apply(combinations, 2, function(x) setdiff(1:V,x)))
    foldobj = lapply(1:V, .xfitfolds_from_foldvec, folds=folds, ordermat=combinations)
  } else if(V==1){
    # degenerate case where we just average across multiple fits to the same data
    foldobj = list(.make_xfit_folds(fold=1, set1 = 1:n, set2 = 1:n, set3 = 1:n))
  }
  foldobj
}

# foldrepeats=3;xfitfolds=5; n=20
#lapply(seq_len(foldrepeats), .xfitsplit,n=n,V=xfitfolds)
#V= xfitfolds;fold=1



## single cross fitting ----
.make_xfit_folds_plugin <- function(fold, set1, set2){
  fold <- list(fold = fold,
               set1 = set1,
               set2 = set2
  )
  fold
}

.xfitfolds_from_foldvec_plugin <- function (r, folds, ordermat)
{
  nfolds <- length(unique(folds))
  set1 <- which(folds %in% ordermat[1:(nfolds-1),r])
  set2 <- which(folds == ordermat[nfolds,r])
  .make_xfit_folds_plugin(r, set1, set2)
}


.xfitsplit_plugin <- function(r=1,n, V=5){
  folds <- rep(seq_len(V), length = n)
  folds <- sample(folds)
  combinations <- combn(V,V-1)
  combinations <- rbind(combinations, apply(combinations, 2, function(x) setdiff(1:V,x)))
  if(V>1) foldobj = lapply(1:V, .xfitfolds_from_foldvec_plugin, folds=folds, ordermat=combinations)
  # degenerate case where we just average across multiple fits to the same data
  if(V==1) foldobj = list(.make_xfit_folds_plugin(fold=1, set1 = 1:n, set2 = 1:n))
  foldobj
}
# foldrepeats=3;xfitfolds=5; n=20
#lapply(seq_len(foldrepeats), .xfitsplit_plugin,n=n,V=xfitfolds)


## safe functions  ----
.safeapply <- function(X, MARGIN, FUN, ...){
  if(!is.null(dim(X))){
    return(apply(X, MARGIN, FUN, na.rm=TRUE, ...))
  } else{
    return(X)
  }
}

.safesweepminus <- function(x, MARGIN, STATS, FUN = "-", check.margin = TRUE, ...){
  if(!is.null(dim(x))){
    return(sweep(x, MARGIN, STATS, FUN, check.margin, ...))
  } else{
    return(0)# assumes subtraction
  }
}
