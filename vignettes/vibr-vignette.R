## ----invisibles, echo=FALSE, results='markup', message=FALSE------------------
library("knitr")

## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----Prelims, eval=FALSE, fig.show='hold'-------------------------------------
#  devtools::install_github("tlverse/sl3")
#  devtools::install_github("alexpkeil1/vibr")

## ----Prelims 2, results='markup', message=FALSE-------------------------------
library("vibr")
library("qgcomp")


## ----first example------------------------------------------------------------
    data(metals, package="qgcomp")
    XYlist = data.frame(metals[,1:23], Y=metals$y+10) # adding constant to Y purely for illustration
    # split 
    spldat <- qgcomp::split_data(XYlist)
    trdat <- spldat$traindata
    vdat <- spldat$validdata
    
    Y_learners = .default_continuous_learners()
    Xbinary_learners = list(sl3::Lrnr_glm$new(name="logit"))
    Xdensity_learners = .default_density_learners(n_bins=c(5))

    X = trdat[,c(1:5,23)] # selecting a subset to make the example run faster
    Y = trdat$Y

    # basic variable importance using targeted maximum likelihood
    vi <- varimp(X=X, 
                 Y=Y, 
                 delta=0.1, 
                 Y_learners = Y_learners,
                 Xdensity_learners=Xdensity_learners[1:2],
                 Xbinary_learners=Xbinary_learners,
                 estimator="TMLE" # targeted maximum likelihood (AIPW, GCOMP. IPW and cross-vit versions of these are all available)
                 )
  vi

## ---- fig.cap = paste("Visualizing the exposure distribution shift implied by the stochastic intervention:", names(X)[c(1,2)])----
vibr::plotshift_dens(vi, X, Acol=1, delta = 0.1)
vibr::plotshift_dens(vi, X, Acol=2, delta = 0.1)

## ---- fig.cap = paste("Visualizing bivariate exposure value shifts implied by the stochastic intervention:", names(X)[2], " with ", names(X)[c(3,6)])----
vibr::plotshift_scatter(vi, Acol=2, Bcol=3, delta = 0.1)
vibr::plotshift_scatter(vi, Acol=2, Bcol=6, delta = 0.1)

## ---- fig.cap = paste("Visualizing the inverse density of exposure value shift implied by the stochastic intervention:", names(X)[c(4,5)])----
vibr::plotshift_wt(vi, X, Acol=4, delta = 0.1)
vibr::plotshift_wt(vi, X, Acol=5, delta = 0.1)

