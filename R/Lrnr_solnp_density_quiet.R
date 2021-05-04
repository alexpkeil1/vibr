#' Nonlinear Optimization via Augmented Lagrange
#'
#' This version is a copy of Lrnr_solnp with additional consideration
#' for users that want explicit control over printed output of solnp
#'
#' This meta-learner provides fitting procedures for density estimation, finding
#' convex combinations of candidate density estimators by minimizing the
#' cross-validated negative log-likelihood loss of each candidate density. The
#' optimization problem is solved by making use of \code{\link[Rsolnp]{solnp}},
#' using Lagrange multipliers. For further details, consult the documentation of
#' the \code{Rsolnp} package.
#'
#' @docType class
#'
#' @import sl3
#' @importFrom R6 R6Class
#' @importFrom stats runif
#' @importFrom data.table setnames
#'
#' @export
#'
#' @keywords data
#'
#' @return Learner object with methods for training and prediction. See
#'  \code{\link{Lrnr_base}} for documentation on learners.
#'
#' @format \code{\link{R6Class}} object.
#'
#' @family Learners
#'
#' @section Parameters:
#' \describe{
#'   \item{\code{trace=0}}{The value of the objective function and the
#'   parameters is printed at every major iteration (default 0).}
#'   \item{\code{tol=0}}{Relative tolerance on feasibility and optimality
#'   (default 1e-8, default in Rsolnp package is 1e-8).}
#'   \item{\code{...}}{Not currently used.}
#' }
#'
Lrnr_solnp_density_quiet <- R6Class(
  classname = "Lrnr_solnp_density_quiet",
  inherit = Lrnr_base, portable = TRUE,
  class = TRUE,
  public = list(
    initialize = function(tol = 1e-5, trace=0, ...) {
      params <- args_to_list()
      super$initialize(params = params, ...)
    }
  ),
  private = list(
    .covariates = NULL,
    .properties = "density",

    .train = function(task) {
      requireNamespace("sl3", quietly = TRUE)
      verbose <- getOption("sl3.verbose")
      params <- self$params
      eval_fun_loss <- function(alphas) {
        sum(-log(as.vector(as.matrix(task$X) %*% alphas)))
      }
      eq_fun <- function(alphas) {
        sum(alphas)
      }
      fit_object <- Rsolnp::solnp(
        stats::runif(ncol(task$X)), eval_fun_loss,
        eqfun = eq_fun, eqB = 1,
        LB = rep(0L, ncol(task$X)),
        control = list(trace = params$trace, tol = params$tol)
        )
      fit_object$coef <- fit_object$pars
      names(fit_object$coef) <- colnames(task$X)
      if (verbose) {
        cat("\ndensity meta-learner fit:\n")
        print(fit_object$coef)
      }
      fit_object$name <- "solnp"
      return(fit_object)
    },

    .predict = function(task = NULL) {
      requireNamespace("sl3", quietly = TRUE)
      requireNamespace("data.table", quietly = TRUE)
      verbose <- getOption("sl3.verbose")
      X <- task$X
      predictions <- rep.int(NA, nrow(X))
      if (nrow(X) > 0) {
        coef <- private$.fit_object$coef
        if (!all(is.na(coef))) {
          predictions <- data.table::data.table(as.matrix(X
                                                          [,
                                                            which(!is.na(coef)),
                                                            drop = FALSE, with = FALSE
                                                          ]) %*%
                                                  coef[!is.na(coef)])
        } else {
          stop("all SL model coefficients are NA.")
        }
        data.table::setnames(predictions, "likelihood")
      }
      return(predictions)
    },
    .required_packages = c("Rsolnp")
  )
)
