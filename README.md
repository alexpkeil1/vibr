# vibr

Estimating variable importance based on a stochastic intervention to increase each predictor by a small amount (continuous) or increase the probability of a predictor by the same small amount (binary). The underlying methodology is based on papers by Ivan Diaz and Mark van der Laan and some improvements suggested by Haneuse and Rotnitsky. 

The basic idea behind variable importance is that one wants to rank predictors based on how strongly they predict the outcome, conditional on other predictors. Several variable importance algorithms are algorithm specific (e.g. random forest node-purity-based importance, linear model coefficients) and not generally applicable. Others are generally applicable. For example, we could compare the change in the predicted outcome from a population wide change from exposed to unexposed for some binary predictor, but such importance measures often map back to contrasts that simply are not well supported by the data. 

The stochastic intervention approach assesses variable importance via small changes to each continuous predictor or the propensity score of each binary predictor, which assesses variable importance in an area of the predictors that we expect to have good support from the data. Thus, this approach focuses on a measure of variable importance that is heavily tied to the population distribution of predictors at hand. The output looks much like output from a set of linear model coefficients, but represents a generalization of these coefficients. A standard linear model coefficient represents the expected change in the mean target value for a one unit change in the predictor <at all values of other predictors>. Here, we replace "one unit" with "delta", where delta is some small value for which we can reasonably assume support in the data, and we also replace "at all values of other predictors" with "at the observed levels of all other predictors." This formulation allows us to model the target flexibly.

Flexible modeling of the target is farmed out to the sl3 package, but some defaults are provided in the vibr package. This package allows assessment of variable importance via inverse probability weighting (which requires additional, possibly flexible, models for each predictor, conditional on all other predictors), g-computation, targeted maximum likelihood estimation and augmented inverse probability weighting (where the latter two require models for the target and the predictors)


## note this package is under rapid development and is bound to be error prone
## Current capabilities
- [x] continuous target variables
- [ ] multivariate target variables
- [ ] binary target variables
- [x] continuous predictors
- [ ] categorical predictors
- [x] binary predictors
- [x] Super learner based prediction for target variables
- [x] Super learner based prediction for predictors
- [x] GLM based prediction for predictors
- [x] GLM based prediction for predictors
- [ ] GLM with individually specified models for each predictor (e.g. omitting some variables, specifying functional forms)
- [ ] GLM with individually specified models for the target  (e.g. specifying functional forms)


## quick start
## variable importance for a set of mainly continuous predictors on a continuous target
    data(metals, package="qgcomp")
    XYlist = list(X=metals[,1:23], Y=metals$y)
    Y_learners = .default_continuous_learners()
    Xbinary_learners = list(Lrnr_stepwise$new(name="SW"))
    Xdensity_learners = .default_density_learners(n_bins=c(10))


    vi <- varimp(X=XYlist$X,Y=XYlist$Y, delta=0.1, 
       Y_learners = Y_learners,
       Xdensity_learners=Xdensity_learners[1:2],
       Xbinary_learners=Xbinary_learners,
       estimator="TMLE")
    vi

# reusing model fits from another vibr object to obtain an alternative estimator
    vi2 <- varimp_refit(vi, X=XYlist$X,Y=XYlist$Y, delta=0.1,
                    estimator="AIPW")

    vi3 <- varimp_refit(vi, X=XYlist$X,Y=XYlist$Y, delta=0.1,
                    estimator="IPW")

    vi2
    vi3


# using weights (under construction - use with caution)
    V = data.frame(wt=runif(nrow(metals)))
    viw <- varimp(X=XYlist$X,Y=XYlist$Y, V=V, delta=0.1, 
       Y_learners = Y_learners,
       Xdensity_learners=Xdensity_learners[1:2], 
       Xbinary_learners=Xbinary_learners,
       estimator="TMLE", weights="wt")
    viw
