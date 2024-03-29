% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/base_identify.R
\name{plotshift_wt}
\alias{plotshift_wt}
\title{Weight-based graphical identifiability diagnostics and visualization (continuous predictors only)}
\usage{
plotshift_wt(vibr_fit, Acol = 1, delta = 0.01, ...)
}
\arguments{
\item{vibr_fit}{a vibr_fit object from varimp}

\item{Acol}{(integer) which column of predictors in call to varimp to diagnose (can only be continuous column of predictors in call to varimp)}

\item{delta}{(numeric, default=0.01) change in each column of predictors in call to varimp corresponding to stochastic intervention}

\item{...}{not used}
}
\description{
This produces a plot of estimated density (horizontal axis) by a ratio of estimated density under a stochastic intervention (intervention exposure = observed exposure -delta) over the estimated density of the observed data. This gives an indication of influential observations (high weights) and allows one to cross-check against estimated density in the observed data.
}
\details{
Note that the "weights" used in vibr use a more difficult calculation than the weight given here (the actual weights used by \code{varimp} depend on whether the stochastic intervention pushes exposure outside of the support of the data), so this function is useful for identifying potentially influential points that may turn out to be a non-issue in analysis.
}
\examples{
\dontrun{
data(metals, package="qgcomp")
# subset of predictors
XYlist = list(X=metals[,1:4], Y=metals$y)
Y_learners = .default_continuous_learners()
Xbinary_learners = list(Lrnr_stepwise$new(name="SW"))
Xdensity_learners = .default_density_learners(n_bins=c(10))
set.seed(1231)
vi_ipw <- varimp(X=XYlist$X,Y=XYlist$Y, delta=0.1, Y_learners = Y_learners,
       Xdensity_learners=Xdensity_learners[1:2], Xbinary_learners=Xbinary_learners,
       estimator="TMLE")
plotshift_wt(vi_ipw, Acol=3, delta=0.01)
# Note many weights are close to 1.0, meaning that the intervention effect
# (in terms of weighted estimators) is based on the relatively smaller number
# of observations with weights > 1.0. This could be reflected by a large
# overlap of the observed and "shifted" density functions, which is reflected
# by the call to plotshift_dens
plotshift_dens(vi_ipw, Acol=3, delta=0.01)
# There is not necessarily an ideal appearance for plotshift_wt plots, but
# a smooth and well estimated density will tend to lead to y-axes where the
# values concentrate away from 1.0 but do not get too large. Thus,
# this plot can diagnose issues in identifiability (really large weights),
# or density estimation (not very many non-zero weights even with large
# value of delta, or really large weights)
}
}
