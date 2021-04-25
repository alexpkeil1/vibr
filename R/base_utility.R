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
.default_continuous_learners <- function(){
  continuous_learners=list(
    Lrnr_glm$new(name="OLS"),
    Lrnr_glmnet$new(name="Lasso", alpha=1.0),
    Lrnr_glmnet$new(name="Enet", alpha=0.0),
    Lrnr_earth$new(name="MARS", alpha=1.0),
    Pipeline$new(Lrnr_pca$new(name="PCA"), Lrnr_glm$new(name="OLS")) # PCA plus glm
  )
  continuous_learners
}

#' @import sl3 glmnet earth stats
#' @export
.default_binary_learners <- function(){
  bin_learners=list(
    Lrnr_glm$new(name="LOGIT"),
    Lrnr_glmnet$new(name="Lasso", alpha=1.0),
    Lrnr_glmnet$new(name="Enet", alpha=0.0),
    Lrnr_earth$new(name="MARS", alpha=1.0),
    Pipeline$new(Lrnr_pca$new(name="PCA"), Lrnr_glm$new(name="LOGIT")) # PCA plus glm
  )
  bin_learners
}

################################
# training functions
################################



#  train all SL model on continuous covariate
#' .train_superlearner
#'
#' @param datatask
#' @param learners
#' @param type
#'
#' @description Train a super learner fit, given some data held in an sl3 task object and a list of learners
#'
#' @return
#' @export
#'
#' @examples
.train_superlearner <- function(datatask, learners, type=c("density", "probability", "expectation")){
  #sl3_Task$new()
  # train learners
  learner_stack <- Stack$new(learners)
  learner_fit <- learner_stack$train(datatask)
  #cross validation
  cv_stack <- Lrnr_cv$new(learner_stack)
  cv_fit <- cv_stack$train(datatask)
  cv_task <- cv_fit$chain()
  # train super learner
  metalearner <- switch(substr(type[1],1,1),
                        d = make_learner(Lrnr_solnp_density),
                        h = make_learner(Lrnr_solnp_density),
                        e = make_learner(Lrnr_nnls),
                        p = make_learner(Lrnr_nnls)
  )
  sl_fit <- metalearner$train(cv_task)
  sl_pipeline <- make_learner(Pipeline, learner_fit, sl_fit)
  sl_pipeline
}


.train_Y <- function(X,Y, learners, verbose=TRUE, isbin=FALSE){
  df = data.frame(X, Y)
  pp1 <- ncol(df)
  isbin <- length(unique(Y))==2
  XY <- sl3_Task$new(data=df,
                     outcome_type=variable_type(type=ifelse(isbin, "binomial", "continuous")),
                     covariates=names(df)[-pp1],
                     outcome=names(df)[pp1]
  )
  if(verbose) cat(paste0("Training: ", names(df)[pp1], "(", ifelse(isbin, "binomial", "continuous, predicted"), ")\n"))
  .train_superlearner(XY,learners, type="expectation")
}


.train_allX <- function(X, tasks, bin_learners, density_learners, verbose=TRUE){
  # TODO: check for data frame X
  vartype = ifelse(apply(X, 2, function(x) length(unique(x))==2), "b", "c")
  p = ncol(X)
  trained_models <- list()
  for(j in 1:p){
    if(verbose) cat(paste0("Training: ", names(X)[j], "(", ifelse(vartype[j]=="b", "binomial", "continuous, density"), ")\n"))
    trained_models[[j]] <- switch(vartype[j],
                                  b = .train_superlearner( tasks[[j]],bin_learners, type="prob"),
                                  c = .train_superlearner( tasks[[j]],density_learners, type="density")
    )
  }
  trained_models
}

################################
# prediction functions
################################


# now create some functions that can automate density prediction for every variable
.gfunction <- function(X=NULL,Acol,gfits=sl.gfits){
  if(!is.null(X)){
    XX <- sl3_Task$new(
      data=data.frame(X),
      covariates=names(X)[-Acol],
      outcome=names(X)[Acol]
    )
    pred <- gfits[[Acol]]$predict(XX)
  }
  if(is.null(X)){
    pred <- gfits[[Acol]]$predict()
  }
  if(typeof(pred)=="list") pred <- pred[[1]]
  pred
}

# outcome prediction
.qfunction <- function(X=NULL,Acol,qfit=sl.qfit){
  if(!is.null(X)){
    XX <- sl3_Task$new(
      data=data.frame(X),
      covariates=names(X)
    )
    pred <- qfit$predict(XX)
  }
  if(is.null(X))  pred <- qfit$predict()
  if(typeof(pred)=="list") pred <- pred[[1]]
  pred
}

