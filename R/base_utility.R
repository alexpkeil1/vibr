.logit <- function(p) log(p) - log1p(-p)
.expit <- function(mu) 1/(1+exp(-mu))

.shift <- function(X,Acol,shift){
  X[,Acol] = X[,Acol] + shift
  X
}

################################
# Data handlers
################################

.create_tasks <- function(X,Y, delta=0.1){
  slX <- list()
  slXpdelta <- list() # X + delta
  slXmdelta <- list() # X-delta
  p <- ncol(X)
  nm = names(X)
  for(j in seq_len(p)){
    isbin = length(unique(X))==2
    tt = ifelse(isbin, "binomial", "continuous")
    slX[[j]] <- sl3_Task$new(outcome_type=variable_type(type=tt), data=data.frame(X),
                             covariates=nm[-j],
                             outcome=nm[j]
    )
    Xt = X
    Xt[j] = ifelse(isbin, 1, X[j]+delta)
    slXpdelta[[j]] <- sl3_Task$new(outcome_type=variable_type(type=tt), data=data.frame(Xt),
                                   covariates=nm[-j],
                                   outcome=nm[j]
    )
    Xt[j] = ifelse(isbin, 0, X[j]-delta)
    slXmdelta[[j]] <- sl3_Task$new(outcome_type=variable_type(type=tt), data=data.frame(Xt),
                                   covariates=nm[-j],
                                   outcome=nm[j]
    )
  }
  # TODO: give this a class
  list(
    slX = slX,
    slXpdelta = slXpdelta, # X + delta
    slXmdelta = slXmdelta # X-delta
  )
}




################################
# default learners
################################


#' @import sl3
#' @export
.default_density_learners <- function(n_bins=c(10, 20, 30), histtypes=c("equal.mass", "equal.length", "dhist")){
  density_learners=list()
  idx  = idx+1
  nm <- paste0("dens_sp_mean_")
  # uniform marginal with kernel density estimator
  density_learners[[idx]] <- Lrnr_density_semiparametric$new(name=nm, mean_learner = Lrnr_mean$new(), var_learner = NULL)
  idx = 1
  for(nb in n_bins){
    for(histtype in histtypes){
      nm <- paste0("hist_RF_", nb, histtype)
      density_learners[[idx]] <- Lrnr_density_discretize$new(name=nm, categorical_learner = Lrnr_randomForest$new(), n_bins = nb, bin_method=histtype)
      idx  = idx+1
      nm = paste0("hist_Unadj_", nb, histtype)
      density_learners[[idx]] <- Lrnr_density_discretize$new(name=nm, categorical_learner = Lrnr_mean$new(), n_bins = nb, bin_method=histtype)
      idx  = idx+1
    }
  }
  density_learners
}

#' @import sl3 glmnet earth
#' @export
.default_continuous_learners_big <- function(){
  continuous_learners=list(
    Lrnr_glm$new(name="OLS", family="gaussian"),
    Lrnr_glmnet$new(name="Lasso", alpha=1.0),
    Lrnr_glmnet$new(name="Enet", alpha=0.0),
    Lrnr_earth$new(name="MARS", alpha=1.0),
    Pipeline$new(Lrnr_pca$new(name="PCA"), Lrnr_glm$new(name="OLS", family="gaussian")) # PCA plus glm
  )
  continuous_learners
}

#' @import sl3 glmnet earth
#' @export
.default_continuous_learners <- function(){
  continuous_learners=list(
    Lrnr_glm$new(name="OLS"),
    Lrnr_earth$new(name="MARS", alpha=1.0)
  )
  continuous_learners
}

#' @import sl3 glmnet earth stats
#' @export
.default_binary_learners_big <- function(){
  bin_learners=list(
    Lrnr_glm$new(name="LOGIT"),
    Lrnr_glmnet$new(name="Lasso", alpha=1.0),
    Lrnr_glmnet$new(name="Enet", alpha=0.0),
    Lrnr_earth$new(name="MARS", alpha=1.0),
    Pipeline$new(Lrnr_pca$new(name="PCA"), Lrnr_glm$new(name="LOGIT")) # PCA plus glm
  )
  bin_learners
}

#' @import sl3 glmnet earth stats
#' @export
.default_binary_learners <- function(){
  bin_learners=list(
    Lrnr_glm$new(name="LOGIT"),
    Lrnr_earth$new(name="MARS", alpha=1.0)
  )
  bin_learners
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




#' Train super learner
#'
#' @param datatask sl3 task
#' @param learners R list of learners
#' @param type type of fit ("density", "probability", or "expectation")
#'
#' @description Train a super learner fit, given some data held in an sl3 task object and a list of learners
#'
#' @return a trained super learner fit (output from sl3::make_learner)
#' @export
#' @import sl3
#'
.train_superlearner <- function(datatask, learners, metalearner=NULL, type=c("density", "probability", "expectation")){
  #sl3_Task$new()
  # train learners
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
  sl_pipeline
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
.train_Y <- function(X,Y, learners, metalearner=NULL, verbose=TRUE, isbin=FALSE){
  df = data.frame(X, Y)
  pp1 <- ncol(df)
  isbin <- length(unique(Y))==2
  XY <- sl3_Task$new(data=df,
                     outcome_type=variable_type(type=ifelse(isbin, "binomial", "continuous")),
                     covariates=names(df)[-pp1],
                     outcome=names(df)[pp1]
  )
  if(verbose) cat(paste0("Training: ", names(df)[pp1], "(", ifelse(isbin, "binomial", "continuous, predicted"), ")\n"))
  .train_superlearner(XY, learners, metalearner=NULL, type="expectation")
}


#' @import sl3
.train_allX <- function(X, tasks, bin_learners, density_learners, metalearner=NULL, verbose=TRUE){
  # TODO: check for data frame X
  vartype = ifelse(apply(X, 2, function(x) length(unique(x))==2), "b", "c")
  p = ncol(X)
  trained_models <- list()
  for(j in 1:p){
    if(verbose) cat(paste0("Training: ", names(X)[j], "(", ifelse(vartype[j]=="b", "binomial", "continuous, density"), ")\n"))
    trained_models[[j]] <- switch(vartype[j],
                                  b = .train_superlearner( tasks[[j]],bin_learners, metalearner=metalearner[[bin_learners]], type="prob"),
                                  c = .train_superlearner( tasks[[j]],density_learners, metalearner=metalearner[[density_learners]], type="density")
    )
  }
  trained_models
}

################################
# prediction functions
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
# shared functions of esimators
################################

# clever covariate/weight continuous
.Haw <- function(gn, ga, gb){
  # I(A<u(w))g(a-delta|w)/g(a|w) + I(a>=u(w)-delta)
  Haw <- ifelse(gn>0, ga/gn, 0) + as.numeric(gb == 0)
  Haw
}

# clever covariate, binary (based on diaz and vdl 2012)
.Hawb_old <- function(gn, ga, gb, shift, X, Acol){
  # if intervention would push exposure out of the support of A | W, then don't intervene
  # (delta*(I(A=1)-I(A=0)) + gn)/gn + gn/gn
  #
  Haw <- ifelse(gn>0, (shift*(2*X[,Acol] - 1)/gn + 1), 0)+ as.numeric(gb == 0)
  Haw
}

# clever covariate, binary (based on diaz and vdl 2018, translated to shift in propensity score)
.Hawb <- function(gn, shift, X, Acol){
  # if intervention would push exposure out of the support of A | W, then don't intervene
  # NOTE: (gn-shift)/gn = 1+(-shift)/gn
  Haw1 <- ifelse(gn>0, 1+(-shift)/gn, 0) + as.numeric(gn + shift > 1)
  Haw0 <- ifelse((1-gn)>0, 1+(shift)/(1-gn), 0) + as.numeric((1-gn) + shift > 1)
  Haw <- X[,Acol]*Haw1 + (1-X[,Acol])*Haw0
  Haw
}


