

#' Density Estimation With Mean Model and Homoscedastic normal errors
#'
#' This learner assumes a mean model with homoscedastic errors: Y ~ E(Y|W) + epsilon. E(Y|W) is fit using a glm,
#' and then the errors are assumed normally distributed epsilon_i ~ Normal(0, sigma_i) where sigma_i is the estimated standard error of the residual.
#'
#' @docType class
#'
#' @importFrom R6 R6Class
#' @importFrom assertthat assert_that is.count is.flag
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
#'   \item{\code{mean_learner}}{The learner to wrap.}
#' }
#'
Lrnr_density_gaussian <- R6Class(
  classname = "Lrnr_density_gaussian",
  inherit = Lrnr_base, portable = TRUE,
  class = TRUE,
  public = list(
    initialize = function(intercept = TRUE, ...) {
      params <- args_to_list()
      super$initialize(params = params, ...)
    }
  ),

  private = list(
    .properties = c("density"),

    .train = function(task) {
      args <- self$params
      outcome_type <- self$get_outcome_type(task)

      if (is.null(args$family)) {
        args$family <- outcome_type$glm_family(return_object = TRUE)
      }
      family_name <- args$family$family
      linkinv_fun <- args$family$linkinv
      link_fun <- args$family$linkfun

      # specify data
      if (args$intercept) {
        args$x <- as.matrix(task$X_intercept)
      } else {
        args$x <- as.matrix(task$X)
      }
      args$y <- outcome_type$format(task$Y)

      if (task$has_node("weights")) {
        args$weights <- task$weights
      }

      if (task$has_node("offset")) {
        args$offset <- task$offset_transformed(link_fun)
      }

      args$control <- glm.control(trace = FALSE)
      fit_object <- .call_with_args_vibr(stats::glm.fit, args)

      fit_object$linear.predictors <- NULL
      fit_object$weights <- NULL
      fit_object$prior.weights <- NULL
      fit_object$y <- NULL
      fit_object$residuals <- NULL
      fit_object$fitted.values <- NULL
      fit_object$effects <- NULL
      fit_object$qr <- NULL
      fit_object$linkinv_fun <- linkinv_fun
      fit_object$link_fun <- link_fun
      fit_object$training_offset <- task$has_node("offset")
      return(fit_object)
      },

    .predict = function(task) {
      verbose <- getOption("sl3.verbose")
      if (self$params$intercept) {
        X <- task$X_intercept
      } else {
        X <- task$X
      }

      predictions <- rep.int(NA, nrow(X))
      if (nrow(X) > 0) {
        coef <- self$fit_object$coef
        if (!all(is.na(coef))) {
          eta <- as.matrix(X
                           [, which(!is.na(coef)),
                             drop = FALSE,
                             with = FALSE
                           ]) %*% coef[!is.na(coef)]

          if (self$fit_object$training_offset) {
            offset <- task$offset_transformed(self$fit_object$link_fun, for_prediction = TRUE)
            eta <- eta + offset
          }

          predictions <- as.vector(self$fit_object$linkinv_fun(eta))
        }
      }
      errors <- task$Y - predictions
      n <- task$nrow
      df <- n - length(coef)
      stdres <- .stdres(resids=errors, df=df, X)

      dens_preds <- dnorm(errors, 0, stdres)
      # dens_preds[is.na(dens_preds)] <- 0
      return(dens_preds)
    },
    .required_packages = c()
  )
)
