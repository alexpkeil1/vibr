#' sl3 extension: Nonlinear Optimization via Augmented Lagrange
#'
#' This version is a copy of sl3::Lrnr_solnp with additional consideration
#' for users that want explicit control over printed output of solnp
#'
#' This meta-learner provides fitting procedures for any pairing of loss
#' function and metalearner function, subject to constraints. The optimization
#' problem is solved by making use of \code{\link[Rsolnp]{solnp}}, using
#' Lagrange multipliers. For further details, consult the documentation of the
#' \code{Rsolnp} package.
#'
#' @docType class
#'
#' @importFrom R6 R6Class
#' @import sl3
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
#'   \item{\code{learner_function=metalearner_linear}}{A function(alpha, X) that
#'     takes a vector of covariates and a matrix of data and combines them into
#'     a vector of predictions. See \link{metalearners} for options.}
#'   \item{\code{loss_function=loss_squared_error}}{A function(pred, truth) that
#'     takes prediction and truth vectors and returns a loss vector. See
#'     \link{loss_functions} for options.}
#'   \item{\code{make_sparse=TRUE}}{If TRUE, zeros out small alpha values.}
#'   \item{\code{convex_combination=TRUE}}{If \code{TRUE}, constrain alpha to
#'     sum to 1.}
#'   \item{\code{init_0=FALSE}}{If TRUE, alpha is initialized to all 0's, useful
#'     for TMLE. Otherwise, it is initialized to equal weights summing to 1,
#'     useful for SuperLearner.}
#'   \item{\code{trace=0}}{The value of the objective function and the
#'   parameters is printed at every major iteration (default 0).}
#'   \item{\code{tol=0}}{Relative tolerance on feasibility and optimality
#'   (default 1e-5, default in Rsolnp package is 1e-8).}
#'   \item{\code{...}}{Not currently used.}
#' }
#'
Lrnr_solnp_quiet <- R6Class(
  classname = "Lrnr_solnp_quiet",
  inherit = Lrnr_base, portable = TRUE,
  class = TRUE,
  public = list(
    initialize = function(learner_function = metalearner_linear,
                          loss_function = loss_squared_error,
                          make_sparse = TRUE, convex_combination = TRUE,
                          init_0 = FALSE, tol = 1e-5, trace=0, ...) {
      params <- args_to_list()
      super$initialize(params = params, ...)
    }
  ),
  private = list(
    .properties = c(
      "continuous", "binomial", "categorical", "weights",
      "offset"
    ),

    .train = function(task) {
      requireNamespace("sl3", quietly = TRUE)
      verbose <- getOption("sl3.verbose")
      params <- self$params
      learner_function <- params$learner_function
      loss_function <- params$loss_function
      outcome_type <- self$get_outcome_type(task)

      # specify data
      X <- as.matrix(task$X)
      Y <- outcome_type$format(task$Y)

      if (task$has_node("offset")) {
        offset <- task$offset
      } else {
        offset <- NULL
      }

      weights <- task$weights
      risk <- function(alphas) {
        requireNamespace("sl3", quietly = TRUE)
        if (!is.null(offset)) {
          preds <- learner_function(alphas, X, offset)
        } else {
          preds <- learner_function(alphas, X)
        }
        losses <- loss_function(preds, Y)
        risk <- weighted.mean(losses, weights)
        return(risk)
      }
      if (params$convex_combination) {
        eq_fun <- function(alphas) {
          sum(alphas)
        }
        eqB <- 1
        LB <- rep(0L, ncol(task$X))
      } else {
        eq_fun <- NULL
        eqB <- NULL
        LB <- NULL
      }
      p <- ncol(X)

      if (params$init_0) {
        init_alphas <- rep(0, p)
      } else {
        init_alphas <- rep(1 / p, p)
      }
      fit_object <- Rsolnp::solnp(
        init_alphas, risk,
        eqfun = eq_fun, eqB = eqB,
        LB = LB,
        control = list(trace = params$trace, tol = params$tol)
      )
      coefs <- fit_object$pars
      names(coefs) <- colnames(task$X)

      if (params$make_sparse) {
        max_coef <- max(coefs)
        threshold <- max_coef / 1000
        coefs[coefs < threshold] <- 0
        coefs <- coefs / sum(coefs)
      }
      fit_object$coefficients <- coefs
      fit_object$training_offset <- task$has_node("offset")
      fit_object$name <- "solnp"
      return(fit_object)
    },

    .predict = function(task = NULL) {
      requireNamespace("sl3", quietly = TRUE)
      verbose <- getOption("sl3.verbose")
      X <- as.matrix(task$X)
      alphas <- self$fit_object$coefficients

      if (self$fit_object$training_offset) {
        predictions <- self$params$learner_function(alphas, X, task$offset)
      } else {
        predictions <- self$params$learner_function(alphas, X)
      }
      return(predictions)
    },
    .required_packages = c("Rsolnp")
  )
)
