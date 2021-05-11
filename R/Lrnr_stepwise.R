

##' Stepwise generalized linear model regression.
##'
##' This learner provides fitting procedures for generalized linear models with
##' forward and backword stepwise regression using glm.fit and step.
##'
##' Documentation is copied directly from stats::step documentation
##'
##' @docType class
##' @importFrom R6 R6Class
##' @importFrom stats glm
##' @importFrom stats step
##' @import sl3
##' @export
##' @keywords data
##' @return Learner object with methods for training and prediction. See \code{\link{Lrnr_base}} for documentation on learners.
##' @format \code{\link{R6Class}} object.
##' @family Learners
##'
##' @section Parameters:
##' \describe{
##'   \item{\code{direction="both"}}{ Passed to stats::step. The mode of stepwise search, can be one of "both", "backward", or "forward", with a default of "both". If the scope argument is missing the default for direction is "backward".
##'   }
##'   \item{\code{trace=0}}{ Passed to stats::step. If positive, information is printed during the running of stepAIC. Larger values may give more information on the fitting process.
##'   }
##'   \item{\code{k=2}}{ Passed to stats::step. The multiple of the number of degrees of freedom used for the penalty. Only k = 2 gives the genuine AIC: k = log(n) is sometimes referred to as BIC or SBC.
##'   }
##'   \item{\code{family="gaussian"}}{ GLM family passed to stats::glm ("gaussian" or "binomial" only). Vibr will make a guess based on the task outcome type if this is not specified.
##'   }
##'   \item{\code{...}}{ Other parameters passed directly to \code{\link[my_package]{my_ml_fun}}. See its documentation for details.
##'   }
##' }
##'
##' @section Methods:
##' \describe{
##' \item{\code{special_function(arg_1)}}{
##'   My learner is special so it has a special function.
##'
##'   \itemize{
##'     \item{\code{arg_1}: A very special argument.
##'    }
##'   }
##'   }
##' }
Lrnr_stepwise <- R6Class(
  classname = "Lrnr_stepwise", inherit = Lrnr_base,
  portable = TRUE, class = TRUE,
  public = list(
    # you can define default parameter values here
    # if possible, your learner should define defaults for all required parameters
    initialize = function(
      direction = "both",
      trace = 0,
      k = 2,
      family = NULL,
      ...
    ) {
      params <- sl3::args_to_list()
      super$initialize(params = params, ...)
    }
  ),
  private = list(
    # list properties your learner supports here.
    # Use sl3_list_properties() for a list of options
    .properties = c(
      "binomial", "continuous", "weights", "offset"
      ),

    # list any packages required for your learner here.
    .required_packages = c("stats"),

    # .train takes task data and returns a fit object that can be used to generate predictions
    .train = function(task) {
      requireNamespace("sl3", quietly = TRUE)
      verbose <- getOption("sl3.verbose")
      args <- self$params

      outcome_type <- self$get_outcome_type(task)
      args$data <- as.data.frame(cbind(Y=task$Y, task$X))
      args$x <- TRUE#as.matrix(task$X)
      args$y <- TRUE#outcome_type$format(task$Y)
      # weights, if any
      if (task$has_node("weights")) {
        args$weights <- task$weights
      } else{
        args$weights <- NULL
      }
      # offset, if any
      if (task$has_node("offset")) {
        args$offset <- task$offset
      } else{
        args$offset <- NULL
      }
      if (is.null(args$family)) {
        args$family <- outcome_type$glm_family(return_object = TRUE)
      }

      args$formula <- as.formula(Y ~ .)
      .step.fun <- function(args){
        fit.glm <- .call_with_args_vibr(glm, args)
        args$object <- fit.glm
        fit.step <- .call_with_args_vibr(step, args)
        fit.step
      }

      # call a function that fits your algorithm
      # with the argument list you constructed
      #fit_object <- .step.fun(Y,X,outcome_type,params$direction, params$trace, params$k)

      suppressWarnings({
        fit_object <- .step.fun(args)
      })

      # return the fit object, which will be stored
      # in a learner object and returned from the call
      # to learner$predict
      return(fit_object)
    },

    # .predict takes a task and returns predictions from that task
    .predict = function(task = NULL) {
      self$training_task
      self$training_outcome_type
      self$fit_object

      predictions <- predict(self$fit_object,
                             newdata = task$X,
                             type="response"
                             )
      return(predictions)
    }
  )
)
