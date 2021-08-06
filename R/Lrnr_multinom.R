#' sl3 extension: Feed-Forward Neural Networks and Multinomial Log-Linear Models
#'
#' This learner provides feed-forward neural networks with a single hidden layer,
#' and for multinomial log-linear models.
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
#' @return Learner object with methods for both training and prediction. See
#'  \code{\link{Lrnr_base}} for documentation on learners.
#'
#' @format \code{\link{R6Class}} object.
#'
#' @family Learners
#'
#' @section Parameters:
#' \describe{
#'   \item{\code{formula}}{A formula of the form class ~ x1 + x2 + ...}
#'   \item{\code{weights}}{(case) weights for each example – if missing defaults to 1}
#'   \item{\code{size}}{number of units in the hidden layer. Can be zero if there are skip-layer units.}
#'   \item{\code{entropy}}{switch for entropy (= maximum conditional likelihood) fitting. Default by least-squares.}
#'   \item{\code{decay}}{parameter for weight decay. Default 0.}
#'   \item{\code{maxit}}{maximum number of iterations. Default 100.}
#'   \item{\code{linout}}{switch for linear output units. Default logistic output units.}
#'   \item{\code{...}}{Other parameters passed to
#'     \code{\link[nnet]{nnet}}.}
#' }
#'
#
Lrnr_multinom <- R6Class(
  classname = "Lrnr_multinom",
  inherit = Lrnr_base, portable = TRUE, class = TRUE,
  public = list(
    initialize = function(decay = 0, maxit = 100, linout = FALSE, ...) {
      super$initialize(params = args_to_list(), ...)
    }
  ),

  private = list(
    .properties = c("binomial", "categorical", "weights"),

    .train = function(task) {
      args <- self$params
      outcome_type <- self$get_outcome_type(task)
      # specify data
      #args$x <- as.data.frame(task$X)
      #args$y <- outcome_type$format(task$Y)
      args$data <- data.frame(outcome_type$format(task$Y), as.data.frame(task$X))

      args$formula <- paste0(names(args$data)[1],"~.")

      if (task$has_node("weights")) {
        args$weights <- task$weights
      }

      if (task$has_node("offset")) {
        args$offset <- task$offset
      }

      fit_object <- .call_with_args_vibr(nnet::multinom, args, keep_all = TRUE)

      # if (self$params$serializeable) {
      #  invisible(fit_object$fit$state)
      # }

      return(fit_object)
    },

    .predict = function(task) {
      outcome_type <- private$.training_outcome_type
      predictions <- predict(private$.fit_object,
                             newdata = data.frame(task$X),
                             type = "probs"
      )
      if (outcome_type$type == "binomial") {
        # extract p(Y=1)
        predictions <- predictions[, 2]
      } else if (outcome_type$type == "categorical") {
        # pack predictions in a single column
        predictions <- pack_predictions(predictions)
      }
      return(predictions)
    },
    .required_packages = c("nnet")
  )
)
