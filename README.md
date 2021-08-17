# vibr

Estimating variable importance based on a stochastic intervention to increase each predictor by a small amount (continuous) or increase the probability of a predictor by the same small amount (binary). The underlying methodology is based on papers by Ivan Diaz and Mark van der Laan and some improvements suggested by Haneuse and Rotnitzky. 

The basic idea behind variable importance is that one wants to rank predictors based on how strongly they predict the target (e.g. health outcome), conditional on other predictors. Several variable importance algorithms are algorithm specific (e.g. random forest node-purity-based importance, linear model coefficients) and not generally applicable. Others are generally applicable. For example, we could compare the change in the predicted outcome from a population wide change from exposed to unexposed for some binary predictor, but such importance measures often map back to contrasts that simply are not well supported by the data. 

The stochastic intervention approach assesses variable importance via small changes to each continuous predictor or the propensity score of each binary predictor, which assesses variable importance in an area of the predictors that we expect to have good support from the data. Thus, this approach focuses on a measure of variable importance that is heavily tied to the population distribution of predictors at hand. The output looks much like output from a set of linear model coefficients, but represents a generalization of these coefficients. A standard linear model coefficient represents the expected change in the mean target value for a one unit change in the predictor <at all values of other predictors>. Here, we replace "one unit" with "delta", where delta is some small value for which we can reasonably assume support in the data, and we also replace "at all values of other predictors" with "at the observed levels of all other predictors." This formulation allows us to model the target flexibly.

Flexible modeling of the target is farmed out to the sl3 package, but some defaults are provided in the vibr package. This package allows assessment of variable importance via inverse probability weighting (which requires additional, possibly flexible, models for each predictor, conditional on all other predictors), g-computation, targeted maximum likelihood estimation and augmented inverse probability weighting (where the latter two require models for the target and the predictors)


note this package is under rapid development: please create an issue if you find something that looks like an error
	
## Current capabilities
- [x] continuous target variables
- [ ] multivariate target variables
- [x] binary target variables
- [x] continuous predictors
- [ ] categorical predictors
- [x] binary predictors
- [x] Super learner based prediction for target variables
- [x] Super learner based prediction for predictors
- [x] GLM based prediction for outcome
- [x] GLM based prediction for predictors
- [x] TMLE
- [x] AIPW
- [x] IPW
- [x] G COMPUTATION
- [x] DOUBLE CROSS FITTING
- [ ] GLM with individually specified models for each predictor (e.g. omitting some variables, specifying functional forms)
- [ ] GLM with individually specified models for the target  (e.g. specifying functional forms)


## quick start
## installation
    devtools::install_github("tlverse/sl3", build_vignettes = TRUE)
    devtools::install_github("alexpkeil1/vibr", build_vignettes = TRUE)

Alternatively, you can clone it locally, open the R project associated with this, and, in R studio, install through the "build" tab/menu item.

## variable importance for a set of mainly continuous predictors on a continuous target
	
    data(metals, package="qgcomp")
    XYlist = data.frame(metals[,1:23], Y=metals$y+10) # adding constant to Y purely for illustration
    # split 
    spldat <- qgcomp::split_data(XYlist)
    trdat <- spldat$traindata
    vdat <- spldat$validdata
    
    Y_learners = .default_continuous_learners()
    Xbinary_learners = list(sl3::Lrnr_glm$new(name="logit"))
    Xdensity_learners = .default_density_learners(n_bins=c(10))

    X = trdat[,1:23]
    Y = trdat$Y

    # basic variable importance using targeted maximum likelihood
    vi <- varimp(X=X, 
                 Y=Y, 
                 delta=0.1, 
                 Y_learners = Y_learners,
                 Xdensity_learners=Xdensity_learners[1:2],
                 Xbinary_learners=Xbinary_learners,
                 estimator="TMLE"
                 )
                 
    # this will generate error messages from sl3 when various learners do not work
    # with the data. These can be safely ignored because sl3 handles failed learners
    # gracefully.

    # now view vi object: if vi measure is > 0, then exposure increases the outcome
    vi

# variable importance/treatment effects for a subset of covariates

    # one can also look at variable importance for a basic subset of the data
    # Note: this uses the full data in X and W, but estimates are only created
    # for X, which can save considerable time if inverse probability weights
    # are required and one is only interested in a small number of specific
    # variables: here we look at "mage35"
    Xsub = X[,23, drop=FALSE]
    W = X[,-c(23), drop=FALSE]
    vi_sub <- varimp(
                 X=Xsub,   
                 W=W,   
                 Y=Y, 
                 delta=0.1, 
                 Y_learners = Y_learners,
                 Xdensity_learners=Xdensity_learners[1:2],
                 Xbinary_learners=Xbinary_learners,
                 estimator="TMLE"
                 )
    vi_sub
    vi$res[23,]
    
    # these might differ by quite a bit - machine learning based esimators
    # can be numerically unstable: we can address this via cross fitting

# Cross fitting, double cross fitting

    # double-cross fitted TMLE can be invoked with estimator = TMLEX
    # this estimator ought to perform asymptotically better
    # and has faster convergence rates than standard TMLE and so
    # may also perform better at moderate sample sizes.
    # estimates will generally be more stable. It will take much more time.
    vi_subXfit <- varimp(
                 X=Xsub,   
                 W=W,   
                 Y=Y, 
                 delta=0.1, 
                 Y_learners = c(Y_learners, list(Lrnr_xgboost$new(), Lrnr_randomForest$new())),
                 Xdensity_learners=Xdensity_learners[1:2],
                 Xbinary_learners=Xbinary_learners,
                 estimator="TMLEX",
                 xfitfolds=3,
                 foldrepeats = 10,
                 )
    vi_subXfit
    
# More stable estimates with machine learning based learners
    
    # note that if you just want a more stable estimator (e.g. if using superlearner
    # and you get substantially different answers between runs, you can average
    # TMLE over multiple runs and calculate a variance based on the average estimated
    # variance across runs plus the variance of the estimates between runs)
    # This is accomplished by setting xfitfolds to 1, so this routine will average
    # over 10 fits (this number is chosen arbitrarily and is not a recommendation)
    #
    vi_sub2avg <- varimp(
                 X=Xsub,   
                 W=W,   
                 Y=Y, 
                 delta=0.1, 
                 Y_learners = c(Y_learners, list(Lrnr_xgboost$new(), Lrnr_randomForest$new())),
                 Xdensity_learners=Xdensity_learners[1:2],
                 Xbinary_learners=Xbinary_learners,
                 estimator="TMLEX",
                 xfitfolds=1,
                 foldrepeats = 10,
                 )
                 
    vi_sub
    vi_subXfit                 
    vi_sub2avg
    
    # can look at estimates of each repeated run (partition) to get a sense of the 
    # stability of the estimate
    # here: the estimate is very stable (but cross fit estimate is likely lower bias)!
    vi_sub2avg$ests
    
    # now examine stability with a continuous variable, which ought to be more unstable

    Xsub2 = X[,4,drop=FALSE]
    W2 = X[,-4,drop=FALSE]
    # single run
    vi_sub2cont <- varimp(
                 X=Xsub2,   
                 W=W2,   
                 Y=Y, 
                 delta=0.1, 
                 Y_learners = c(Y_learners, list(Lrnr_xgboost$new(), Lrnr_randomForest$new())),
                 Xdensity_learners=Xdensity_learners,
                 Xbinary_learners=Xbinary_learners,
                 estimator="TMLE",
                 )
    vi_sub2cont
    
    # cross fit
    vi_sub2Xcont <- varimp(
                 X=Xsub2,   
                 W=W2,   
                 Y=Y, 
                 delta=0.1, 
                 Y_learners = c(Y_learners, list(Lrnr_xgboost$new(), Lrnr_randomForest$new())),
                 Xdensity_learners=Xdensity_learners,
                 Xbinary_learners=Xbinary_learners,
                 estimator="TMLEX",
                 xfitfolds=3,
                 foldrepeats = 5,
                )
    vi_sub2Xcont
    
    # multi-run average
    vi_sub2avgcont <- varimp(
                 X=Xsub2,   
                 W=W2,   
                 Y=Y, 
                 delta=0.1, 
                 Y_learners = c(Y_learners, list(Lrnr_xgboost$new(), Lrnr_randomForest$new())),
                 Xdensity_learners=Xdensity_learners,
                 Xbinary_learners=Xbinary_learners,
                 estimator="TMLEX",
                 xfitfolds=1,
                 foldrepeats = 5,
                 )

    vi_sub2cont
    vi_sub2avgcont
    vi_sub2Xcont

  
    # variation between runs (multirun vs. cross fit)
    # variation of in sample estimates
     vi_sub2avgcont$ests   
     sd(vi_sub2avgcont$ests)
    # variation of out of sample estimates
     vi_sub2Xcont$ests
     sd(vi_sub2Xcont$ests)

    vi_sub2avgcont$vars # variance at each run
    vi_sub2Xcont$vars # variance at each run, cross fit

    

    # and this can be done with the entire covariate set to get more stable
    # estimates of variable importance + ranking for all variables at once
    # this will take a while
    vi_avg <- varimp(
                 X=X,   
                 W=NULL,   
                 Y=Y, 
                 delta=0.1, 
                 Y_learners = Y_learners,
                 Xdensity_learners=Xdensity_learners[1:2],
                 Xbinary_learners=Xbinary_learners,
                 estimator="TMLEX",
                 xfitfolds=1,
                 foldrepeats = 5,
                 )
    # compare with original             
    vi_avg
    vi
    vi_avg$ests

# utilizing existing fits to efficiently (compute-wise) estimate new quantities

    # estimate mean outcome after each "intervention" (p value not very useful here)
    mean(Y)
    varimp_refit(vi,X=X, 
                 Y=Y, 
                 delta=0.1, 
                 estimator="TMLE",
                 estimand="mean"
                 )
    
    # change intervention strength
    vibr::plotshift_dens(vi, 1, delta=0.1)
    vibr::plotshift_scatter(vi, Acol=1, Bcol=2, delta=0.1)
    varimp_refit(vi,X=X, 
                 Y=Y, 
                 delta=0.01, 
                 estimator="TMLE",
                 estimand="diff"
                 )

    # estimate using outcome model only (NOTE STD ERROR IS NOT VALID)
    varimp_refit(vi,X=X, 
                 Y=Y, 
                 delta=0.1, 
                 estimator="GCOMP",
                 estimand="diff"
                 )

# bootstrapping for better standard errors (can be combined with stabilization/cross fitting)

    viB <- varimp(X=X, 
                 Y=Y, 
                 delta=0.1, 
                 Y_learners = Y_learners,
                 Xdensity_learners=NULL,
                 Xbinary_learners=NULL,
                 estimator="GCOMP",
                 B=10
                 )
    viB

# utilizing underlying (nuisance) models

    # getting predictions (not bootstrapped)
    outcome_model <- vi$qfit
    outcome_model$predict() # predictions at observed values
    plot(outcome_model$predict(), Y)
    
    # prediction in new data
    Xv = vdat[,1:23]
    Yv = vdat$Y
    
    #outcome_model$predict(valid_task) # = predictions in validation data based on model fit to training data
    
    valid_task = sl3::sl3_Task$new(data = data.frame(Xv,Y=Yv), outcome="Y", covariates=names(Xv))
    library(ggplot2)
    ggplot(data= data.frame(pred=outcome_model$predict(valid_task), y=Yv), aes(x=y, y=pred)) + 
      geom_point() + geom_smooth()
    



# using weights (under construction - use with caution)

    V = data.frame(wt=runif(nrow(metals)))
    viw <- varimp(X=XYlist$X,Y=XYlist$Y, V=V, delta=0.1, 
       Y_learners = Y_learners,
       Xdensity_learners=Xdensity_learners[1:2], 
       Xbinary_learners=Xbinary_learners,
       estimator="TMLE", weights="wt")
    viw



# binomial outcome
    Y = metals$disease_state[spldat$trainidx]
    Ybinary_learners = list(
      sl3::Lrnr_glm$new(name="logit"),
      Lrnr_stepwise$new(name="stepwise"),
      Lrnr_polspline_quiet$new(name="mars")
      )
    vib <- varimp(X=X, 
                 Y=Y, 
                 delta=0.1, 
                 Y_learners = Ybinary_learners,
                 Xdensity_learners=Xdensity_learners[1:2],
                 Xbinary_learners=Xbinary_learners,
                 estimator="TMLE", 
                 family="binomial"
                 )
    vib

    vib2 <- varimp(X=X[,1:3], 
                 Y=Y, 
                 delta=0.1, 
                 Y_learners = Ybinary_learners,
                 Xdensity_learners=Xdensity_learners[1:2],
                 Xbinary_learners=Xbinary_learners,
                 estimator="GCOMP", 
                 family="binomial"
                 )
   vib2

	
  - Diaz Muñoz I, van der Laan M. Population intervention causal effects based on stochastic interventions. Biometrics. 2012;68(2):541–549. 
  - Diaz I, Hubbard A, Decker A, et al. Variable importance and prediction methods for longitudinal problems with missing variables. PloS one. 2015;10(3):e0120031. 
  - Van der Laan MJ, Rose S. Targeted learning: causal inference for observational and experimental data. Springer Science & Business Media; 2011.
  - Haneuse S, Rotnitzky A. Estimation of the effect of interventions that modify the received treatment. Statistics in Medicine. 2013;32(30):5260–5277. 
