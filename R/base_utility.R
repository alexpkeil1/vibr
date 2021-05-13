.logit <- function(p) log(p) - log1p(-p)
.expit <- function(mu) 1/(1+exp(-mu))
.identity <- function(x) x
.invidentity <- .identity

.shift <- function(X,Acol,shift){
  X[,Acol] = X[,Acol] + shift
  X
}

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


################################
# Data handlers
################################
#' @importFrom sl3 sl3_Task variable_type
.create_tasks <- function(X,Y, V=NULL, delta=0.1, ...){
  X <- as.data.frame(X)
  p <- ncol(X)
  nm = names(X)
  slX <- list()
#  slXpdelta <- list() # X + delta
#  slXmdelta <- list() # X-delta
  if(!is.null(V[1])){
    X = data.frame(cbind(X, V)) # V can hold weights, offsets, etc.
  }
  for(j in seq_len(p)){
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

.scale_continuous <- function(X){
  p <- ncol(X)
  isbin_vec <- apply(X, 2, function(x) length(unique(x))==2)
  for(j in seq_len(p)){
    if(!isbin_vec[j]){
      X[,j] <- X[,j]/(2*sd(X[,j]))
    }
  }
  X
}


################################
# default learners
################################


#' @import sl3
#' @export
.default_density_learners_big <- function(n_bins=c(7, 12), histtypes=c("equal.mass", "equal.length")){
  density_learners=list()
  idx  = 1
  nm <- paste0("hist_multinom_10")
  density_learners[[idx]] <- Lrnr_density_discretize$new(name=nm, categorical_learner = Lrnr_multinom$new(trace=FALSE), n_bins = 10, bin_method="equal.mass")
  idx  = idx+1
  nm <- paste0("dens_hse_glm")
  density_learners[[idx]] <- Lrnr_density_hse$new(name=nm, mean_learner = Lrnr_glm$new(), var_learner = Lrnr_glm$new())
  idx  = idx+1
  nm <- paste0("dens_hse_lasso")
  density_learners[[idx]] <- Lrnr_density_hse$new(name=nm, mean_learner = Lrnr_glmnet$new(name="Lasso", alpha=1.0, family="gaussian"), var_learner = Lrnr_glmnet$new(name="Lasso", alpha=1.0, family="gaussian"))
  idx  = idx+1
  nm <- paste0("dens_hse_sw")
  density_learners[[idx]] <- Lrnr_density_hse$new(name=nm, mean_learner = Lrnr_stepwise$new(name="Stepwise"), var_learner = Lrnr_stepwise$new(name="Stepwise"))
  idx  = idx+1
  nm <- paste0("dens_hse_mars")
  density_learners[[idx]] <- Lrnr_density_hse$new(name=nm, mean_learner = Lrnr_earth$new(name="MARS", family=gaussian()), var_learner = Lrnr_glmnet$new(name="Lasso", alpha=1.0, family="gaussian"))
  idx = idx + 1
  nm <- paste0("dens_hse_nnet")
  density_learners[[idx]] <- Lrnr_density_hse$new(name=nm, mean_learner = Lrnr_nnet$new(name="nnet",maxit=200, trace=FALSE), var_learner = Lrnr_glmnet$new(name="Lasso", alpha=1.0, family="gaussian"))
  idx  = idx+1
  nm <- paste0("dens_sp_unadj")
  # uniform marginal with kernel density estimator
  density_learners[[idx]] <- Lrnr_density_semiparametric$new(name=nm, mean_learner = Lrnr_mean$new(), var_learner = NULL)
  idx  = idx+1
  nm <- paste0("dens_sp_glm")
  density_learners[[idx]] <- Lrnr_density_semiparametric$new(name=nm, mean_learner = Lrnr_glm$new(), var_learner = NULL)
  idx  = idx+1
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
  .default_density_learners_big(...)[1:4]
}


#' @import sl3 glmnet earth stats nnet xgboost stats
#' @export
.default_continuous_learners_big <- function(){
  continuous_learners=list(
    Lrnr_glm$new(name="ols", family=gaussian()),
    Lrnr_glmnet$new(name="lasso", alpha=1.0, family="gaussian"),
    Lrnr_stepwise$new(name="stepwise", family=gaussian()),
    Lrnr_earth$new(name="mars", family=gaussian()),
    Lrnr_xgboost$new(name="xgboost"),
    Lrnr_nnet$new(name="nnet",maxit=200, trace=FALSE),
    Lrnr_glmnet$new(name="enet", alpha=0.0, family="gaussian"),
    Pipeline$new(Lrnr_pca$new(name="pca"), Lrnr_glm$new(name="ols", family=gaussian())), # PCA plus glm
    Pipeline$new(Lrnr_screener_coefs$new(name="lassocreen", learner=Lrnr_glmnet$new(name="lasso", alpha=1.0, family="gaussian")), Lrnr_earth$new(name="MARS", family=gaussian())), # screen by coefficient size then OLS
    Pipeline$new(Lrnr_screener_coefs$new(name="coefscreen", learner=Lrnr_glm$new(name="ols", family=gaussian())), Lrnr_glm$new(name="OLS", family=gaussian())), # screen by coefficient size then OLS
    Pipeline$new(Lrnr_screener_importance$new(name="rfimpscreen", learner=Lrnr_randomForest$new()), Lrnr_glm$new(name="OLS", family=gaussian())) # screen by variable importance then OLS
  )
  continuous_learners
}

#' @export
.default_continuous_learners <- function(){
  .default_continuous_learners_big()[1:4]
}

#' @import sl3 glmnet earth stats nnet xgboost stats
#' @export
.default_binary_learners_big <- function(){
  bin_learners=list(
    Lrnr_glm$new(name="logit", family=binomial()),
    Lrnr_glmnet$new(name="lasso", alpha=1.0, family="binomial"),
    Lrnr_stepwise$new(name="stepwise", family=binomial()),
    Lrnr_earth$new(name="mars", family=binomial()),
    Lrnr_xgboost$new(),
    Lrnr_nnet$new(name="nnet",maxit=200, trace=FALSE),
    Lrnr_glmnet$new(name="enet", alpha=0.0, family="binomial"),
    Pipeline$new(Lrnr_pca$new(name="pca"), Lrnr_glm$new(name="logit", family=binomial())), # PCA plus glm
    Pipeline$new(Lrnr_screener_coefs$new(name="lassocreen", learner=Lrnr_glmnet$new(name="lasso", alpha=1.0, family="binomial")), Lrnr_earth$new(name="mars", family=binomial())), # screen by coefficient size then OLS
    Pipeline$new(Lrnr_screener_coefs$new(name="coefscreen", learner=Lrnr_glm$new(name="logit", family=binomial())), Lrnr_glm$new(name="logit", family=binomial())), # screen by coefficient size then OLS
    Pipeline$new(Lrnr_screener_importance$new(name="rfimpscreen", learner=Lrnr_randomForest$new()), Lrnr_glm$new(name="logit", family=binomial())) # screen by variable importance then LOGIT
  )
  bin_learners
}

#' @export
.default_binary_learners <- function(){
  .default_binary_learners_big()[1:4]
}

#' @import sl3 glmnet earth stats
#' @export
.default_metalearner <- function(type){
  switch(substr(type[1],1,1),
                        d = Lrnr_solnp_density_quiet$new(trace=0),
                        h = Lrnr_solnp_density_quiet$new(trace=0),
                        e = Lrnr_nnls$new(convex=TRUE),
                        p = Lrnr_nnls$new(convex=TRUE)
  )
}



################################
# training functions
################################




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
  isbin <- length(unique(Y))==2
  XY <- sl3_Task$new(data=df,
                     outcome_type=variable_type(type=ifelse(isbin, "binomial", "continuous")),
                     covariates=nms[-pp1],
                     outcome=nms[pp1], ...
  )
  if(verbose) cat(paste0("Training: ", names(df)[pp1], "(", ifelse(isbin, "binomial", "continuous, predicted"), ")\n"))
  .train_fullstack(XY, learners, metalearner=NULL, type="expectation")
}


#' @import sl3
.train_allX <- function(X, tasks, bin_learners, density_learners, metalearner=NULL, verbose=TRUE, V=NULL){
  # TODO: check for data frame X
  vartype = ifelse(apply(X, 2, function(x) length(unique(x))==2), "b", "c")
  p = ncol(X)
  trained_models <- list()
  for(j in 1:p){
    if(verbose) cat(paste0("Training: ", names(X)[j], "(", ifelse(vartype[j]=="b", "binomial", "continuous, density"), ")\n"))
    trained_models[[j]] <- switch(vartype[j],
                                  b = .train_fullstack( tasks[[j]],bin_learners, metalearner=metalearner[[bin_learners]], type="prob"),
                                  c = .train_fullstack( tasks[[j]],density_learners, metalearner=metalearner[[density_learners]], type="density")
    )
  }
  trained_models
}

################################
## prediction functions ##
################################


# now create some functions that can automate density prediction for every variable
#.gfunction <- function(X=NULL,Acol=1,gfits=sl.gfits, ...){
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

.gfunction_glm_density <- function(X=NULL,Acol=1,gfits=NULL, ...){
  # assume normal model
  yy = X[,-Acol,drop=TRUE]
  if(!is.null(X)){
    if(!is.data.frame(X)){
      X = as.data.frame(X[,-Acol, drop=FALSE])
    }
    preds <- predict(gfits[[Acol]], newdata = X)
  }
  if(is.null(X)){
    preds <- predict(gfits[[Acol]])
  }
  if(typeof(preds)=="list") preds <- preds[[1]]
  err <- y - preds
  sderr <- sd(err)
  dens <- dnorm(err, 0, sderr)
  dens
}

.gfunction_glm <- function(X=NULL,Acol=1,gfits=NULL, ...){
  if(!is.null(X)){
    if(!is.data.frame(XX)){
      X = as.data.frame(X[,-Acol, drop=FALSE])
    }
    pred <- predict(gfits[[Acol]], newdata = X)
  }
  if(is.null(X)){
    pred <- predict(gfits[[Acol]])
  }
  if(typeof(pred)=="list") pred <- pred[[1]]
  pred
}

# outcome prediction
#.qfunction <- function(X=NULL,Acol=1,qfit=sl.qfit, ...){
.qfunction_sl <- function(X=NULL,Acol=1,qfit=NULL, ...){
    if(!is.null(X)){
    XX <- sl3_Task$new(
      data=data.frame(X),
      covariates=names(X),
      ...
    )
    pred <- qfit$predict(XX)
  }
  if(is.null(X))  pred <- qfit$predict()
  if(typeof(pred)=="list") pred <- pred[[1]]
  pred
}


.qfunction_glm <- function(X=NULL,Acol=1,gfits=NULL, ...){
  if(!is.null(X)){
    if(!is.data.frame(XX)){
      X = as.data.frame(X[,, drop=FALSE])
    }
    pred <- predict(gfits[[Acol]], newdata = X)
  }
  if(is.null(X)){
    pred <- predict(gfits[[Acol]])
  }
  if(typeof(pred)=="list") pred <- pred[[1]]
  pred
}

.gfunction <- function(X=NULL,Acol=1,gfits=NULL, ...){
  #.gfunction_glm(X,Acol,gfits, ...)
  #.gfunction_glm_density(X,Acol,gfits, ...)
  .gfunction_sl(X,Acol,gfits, ...)
}

.qfunction <- function(X=NULL,Acol=1,qfit=NULL, ...){
  #.qfunction_glm(X,1,qfit, ...)
  .qfunction_sl(X,1,qfit, ...)
}


################################
## shared functions of esimators ##
################################

# clever covariate/weight continuous
.Haw <- function(gn, ga, gb){
  # if intervention would push exposure out of the support of A | W, then don't intervene (potential = observed outcome)
  # g(a-delta) = expected density after probabilistically shifting all individuals to a-delta (continuous)
  # i.e. if intervetion is g(a+delta), then observed individual represents g(a-delta) individuals in intervened world
  # I(A<u(w))g(a-delta|w)/g(a|w) + I(a>=u(w)-delta)
  Haw <- ifelse(gn>0, ga/gn, 0) + as.numeric(gb == 0)
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
  if(any(Haw<0)) warning(paste0("",sum(Haw<0), " weights (binary predictor) were < 0\n"))
  cbind(Haw, Haw1, Haw0)[,1:retcols]
}


##########
# prelim functions
#########
.prelims <- function(X, Y, V=NULL, delta, Y_learners, Xbinary_learners, Xdensity_learners, verbose=verbose, ...){
  #.prelims(X, Y, delta, Y_learners, Xbinary_learners, Xdensity_learners, verbose=verbose, ...)
  # basic error checking
  stopifnot(is.data.frame(X))
  if(!is.null(V[1])) stopifnot(is.data.frame(V))

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
    if(verbose) cat("Normalizing weights to sum to length(Y)")
    wt <- wt/mean(wt)
  }
  tasklist = .create_tasks(X,Y,V,delta, ...)
  isbin <- as.character((length(unique(Y))==2))
  if(verbose) cat(paste0("delta = ", delta, "\n")) # TODO: better interpretation
  yb = .bound_zero_one(Y)
  Ybound = yb[[1]]
  sl.qfit <- sl.gfits <- NULL
  if(!is.null(Y_learners)){
    #yb = .bound_zero_one(Y)
    sl.qfit <- .train_Y(X,Y, Y_learners, verbose=verbose, isbin, V=V, ...)
  }
  if(!is.null(Xbinary_learners) || !is.null(Xdensity_learners)){
    sl.gfits <- .train_allX(X, tasklist$slX, Xbinary_learners, Xdensity_learners, verbose=verbose, V=V)
  }
  list(n=n, Ybound=Ybound, tasklist = tasklist, sl.qfit = sl.qfit, sl.gfits=sl.gfits, isbin = isbin, weights=wt)
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
       sl.gfits=obj$gfit,
       tasklist=tasklist,
       binomial=obj$binomial,
       oldtype=obj$type,
       scaled =obj$scaled,
       weights=obj$weights,
       n = length(obj$weights)
       )
}


.attach_misc <- function(obj, scale_continuous, delta){
  obj$scaled = scale_continuous
  obj$delta = delta
  obj
}
