#' @include class_ITRSurvStep.R
setClass(Class = "ITRSurv",
         slots = c("call" = "call",
                   "phaseResults" = "list",
                   "value" = "ANY",
                   "params" = "ANY")) # change back to "Parameters" after fixing!!!

#-------------------------------------------------------------------------------
# method to print key results to screen
#-------------------------------------------------------------------------------
# method is exported
#-------------------------------------------------------------------------------
#' Print Analysis Results
#'
#' Prints the key results of the analysis.
#'
#' @param x A ITRSurv object. The value returned by itrSurv().
#'
#' @param ... Ignored.
#'
#' @export
#' @name print
#' @aliases print,ITRSurv-method
#' @returns No return value, called to display key results.
#' @examples
#'
#'
#' dt <- data.frame("Y.1" = sample(1:100,100,TRUE), "Y.2" = sample(101:200,100,TRUE),
#'                  "D.1" = rbinom(100, 1, 0.9), "D.2" = rbinom(100,1,0.9),
#'                  "A.1" = rbinom(100, 1, 0.5), "A.2" = rbinom(100,1,0.5),
#'                  "X.1" = rnorm(100), "X.2" = rnorm(100))
#'
#' result <- itrSurv(data = dt,
#'                   txName = c("A.1", "A.2"),
#'                   models = list(Surv(Y.1,D.1)~X.1+A.1,
#'                                 Surv(Y.2,D.2)~X.2+A.2+Y.1))
#'
#' print(x = result)
setMethod(f = "print",
          signature = c(x = "ITRSurv"),
          definition = function(x, ...) {

            # message('TESTINGGG0')

              cat("\nCall:\n")
              print(x = x@call)

              cat("\n")

              cat("TOTAL: ", "N = ", length(x@phaseResults[[1]]@optimal@Ratio_Stopping_Ind), "\n",
                  "  tx ", x@phaseResults[[ 1 ]]@txName, "\n",
                  "  tx options: ", x@phaseResults[[ 1 ]]@txLevels, "\n")
              cat("PHASE 1", "\n",
                  "  Surv only (N = ", sum(x@phaseResults[[1]]@optimal@Ratio_Stopping_Ind), ")", "\n")
              cat("PHASE 2", "\n",
                  "  Endpoint (N = ", sum(x@phaseResults[[1]]@optimal@Ratio_Stopping_Ind == 0),")", "\n")

              cat("\n")

              print(class(x@params))
              if (is(x = x@params, class2 = "Param_SurvivalMeanEndPointMean")) {
                cat("Criterion: Truncated Mean Survival Time\n")
                if (x@call[['endPoint']] == "CR"){
                  cat("Criterion: Truncated Mean CIF Time\n")
                }
              }
              else if (is(x = x@params, class2 = "Param_SurvivalProbabilityEndPointMean")){
                cat("Criterion: Survival Probability at T=",
                    x@params@survivalparam@survivalTime,
                    " surv.", x@params@survivalparam@type, "\n")
                if (x@call[['endPoint']] == "CR"){
                  cat("Criterion: Truncated Mean CIF Time\n")
                }
              } else if (is(x = x@params, class2 = "Param_SurvivalMeanEndPointProbability")){
                cat("Criterion: Truncated Mean Survival Time\n")
                if (x@call[['endPoint']] == "CR"){
                  cat("Criterion: CIF Probability at T=",
                      x@params@endpointparam@CIFTime,
                      " cif.", x@params@endpointparam@type, "\n")
                }
              } else if (is(x = x@params, class2 = "Param_SurvivalProbabilityEndPointProbability")){
                cat("Criterion: Survival Probability at T=",
                    x@params@survivalparam@survivalTime,
                    " surv.", x@params@survivalparam@type, "\n")
                if (x@call[['endPoint']] == "CR"){
                  cat("Criterion: CIF Probability at T=",
                      x@params@endpointparam@CIFTime,
                      " cif.", x@params@endpointparam@type, "\n")
                }
              } else{
                stop("params are not right class2")
              }

              for (val in 1L:length(x = x@value)){
                cat(sprintf("Estimated Value V%s: ",val))
                for (i in 1L:length(x = x@value[[val]])) {
                  cat(" ", names(x@value[[val]])[i], ": ",
                      round(x = x@value[[val]][[ i ]], digits = 4), "  ")
                }
                cat("\n")
              }
              cat("\n")

            })

#-------------------------------------------------------------------------------
# method to show key results to screen
#-------------------------------------------------------------------------------
# method is exported
#-------------------------------------------------------------------------------
#' Show Analysis Results
#'
#' Shows the key results of the analysis.
#'
#' @param object A ITRSurv object. The value returned by itrSurv().
#'
#' @export
#' @name show
#' @aliases show,ITRSurv-method
#' @returns No return value, called to display key results.
#' @examples
#'
#'
#' dt <- data.frame("Y.1" = sample(1:100,100,TRUE), "Y.2" = sample(101:200,100,TRUE),
#'                  "D.1" = rbinom(100, 1, 0.9), "D.2" = rbinom(100,1,0.9),
#'                  "A.1" = rbinom(100, 1, 0.5), "A.2" = rbinom(100,1,0.5),
#'                  "X.1" = rnorm(100), "X.2" = rnorm(100))
#'
#' result <- itrSurv(data = dt,
#'                   txName = c("A.1", "A.2"),
#'                   models = list(Surv(Y.1,D.1)~X.1+A.1,
#'                                 Surv(Y.2,D.2)~X.2+A.2+Y.1))
#'
#' show(object = result)
setMethod(f = "show",
          signature = c(object = "ITRSurv"),
          definition = function(object) {

              cat("\nCall:\n")
              print(x = object@call)

              cat("\n")

              cat("TOTAL: ", "N = ", length(object@phaseResults[[1]]@optimal@Ratio_Stopping_Ind), "\n",
                  "  tx ", object@phaseResults[[ 1 ]]@txName, "\n",
                  "  tx options: ", object@phaseResults[[ 1 ]]@txLevels, "\n")
              cat("PHASE 1", "\n",
                    "  Surv only (N = ", sum(object@phaseResults[[1]]@optimal@Ratio_Stopping_Ind), ")", "\n")
              cat("PHASE 2", "\n",
                    "  Endpoint (N = ", sum(object@phaseResults[[1]]@optimal@Ratio_Stopping_Ind == 0),")", "\n")

              cat("\n")

              print(class(object@params))
              if (is(object = object@params, class2 = "Param_SurvivalMeanEndPointMean")) {
                cat("Criterion: Truncated Mean Survival and CIF Time\n")
              }
              else if (is(object = object@params, class2 = "Param_SurvivalProbabilityEndPointMean")){
                cat("Criterion: Survival Probability at T=",
                          object@params@survivalparam@survivalTime,
                          " surv.", object@params@survivalparam@type, "\n")
                cat("Criterion: Truncated Mean CIF Time\n")
              } else if (is(object = object@params, class2 = "Param_SurvivalMeanEndPointProbability")){
                cat("Criterion: Truncated Mean Survival Time\n")
                cat("Criterion: CIF Probability at T=",
                      object@params@endpointparam@CIFTime,
                      " cif.", object@params@endpointparam@type, "\n")
              } else if (is(object = object@params, class2 = "Param_SurvivalProbabilityEndPointProbability")){
                cat("Criterion: Survival Probability at T=",
                    object@params@survivalparam@survivalTime,
                    " surv.", object@params@survivalparam@type, "\n")
                cat("Criterion: CIF Probability at T=",
                    object@params@endpointparam@CIFTime,
                    " cif.", object@params@endpointparam@type, "\n")
              } else{
                stop("params are not right class2")
              }


              cat("\n")

              for (val in 1L:length(x = object@value)){
                cat(sprintf("Estimated Value V%s: ",val))
                for (i in 1L:length(x = object@value[[val]])) {
                  cat(" ", names(object@value[[val]])[i], ": ",
                      round(x = object@value[[val]][[ i ]], digits = 4), "  ")
                }
                cat("\n")
              }
              cat("\n")
              })

#' #-------------------------------------------------------------------------------
#' # method to return key stage results as a list
#' #-------------------------------------------------------------------------------
#' # method is exported
#' #-------------------------------------------------------------------------------
#' #' Retrieve Stage Results as a List
#' #'
#' #' Returns the key results from all stages or one stage of the Q-learning algorithm.
#' #'
#' #' @param object A ITRSurv object. The value returned by itrSurv().
#' #'
#' #' @param ... Ignored. Used to require named inputs.
#' #'
#' #' @export
#' #' @name stage
#' #' @rdname stage
#' setGeneric(name = "stage",
#'            def = function(object, ...) { standardGeneric("stage") })
#'
#' #' Retrieve Stage Results as a List
#' #'
#' #' @rdname itrSurv-internal-api
#' #'
#' setMethod(f = "stage",
#'           signature = c(object = "ANY"),
#'           definition = function(object, ...) { stop("not allowed") })
#'
#' #-------------------------------------------------------------------------------
#' # method to return key stage results as a list
#' #-------------------------------------------------------------------------------
#' # method is exported
#' #-------------------------------------------------------------------------------
#' #' Retrieve Stage Results for Decision Point q as a List
#' #'
#' #' @param object A ITRSurv object. The value returned by itrSurv().
#' #'
#' #' @param ... Ignored. Used to require named inputs.
#' #'
#' #' @param q An integer object. (optional) The stage for which results are
#' #'   desired. If q is not provided, all stage results will be returned.
#' #'
#' #' @return A list object containing the key results for each requested stage.
#' #'   If q is not provided, a list of these results will be returned, where the
#' #'   ith element of that list corresponds to the ith decision point.
#' #'
#' #' @rdname stage
#' #' @examples
#' #'
#' #'
#' #' dt <- data.frame("Y.1" = sample(1:100,100,TRUE), "Y.2" = sample(101:200,100,TRUE),
#' #'                  "D.1" = rbinom(100, 1, 0.9), "D.2" = rbinom(100,1,0.9),
#' #'                  "A.1" = rbinom(100, 1, 0.5), "A.2" = rbinom(100,1,0.5),
#' #'                  "X.1" = rnorm(100), "X.2" = rnorm(100))
#' #'
#' #' result <- itrSurv(data = dt,
#' #'                   txName = c("A.1", "A.2"),
#' #'                   models = list(Surv(Y.1,D.1)~X.1+A.1,
#' #'                                 Surv(Y.2,D.2)~X.2+A.2+Y.1))
#' #'
#' #' tt <- stage(object = result)
#' setMethod(f = "stage",
#'           signature = c(object = "ITRSurv"),
#'           definition = function(object, ..., q) {
#'             message('TESTINGGG3')
#'               if (missing(x = q)) {
#'                 res <- list()
#'                 for (i in 1L:length(object@phaseResults)) {
#'                   res[[ i ]] <- .stage(object = object@phaseResults[[ i ]])
#'                 }
#'                 return( res )
#'               }
#'               if (q > length(x = object@phaseResults)) {
#'                 stop("q > # of decision points of analysis", call. = FALSE)
#'               }
#'               return( .stage(object = object@phaseResults[[ q ]]) )
#'             })


#' Prediction Method
#'
#' Method to estimate the value for new data or to retrieve estimated value for
#'  training data
#'
#' @param object A ITRSurv object. The object returned by a call to itrSurv().
#'
#' @param ... Ignored. Used to require named inputs.
#'
#' @param newdata NULL or a data.frame object. If NULL, this method retrieves
#'   the estimated value for the training data. If a data.frame, the
#'   value is estimated based on the data provided.
#'
#' @param stage An integer object. The stage for which predictions are desired.
#'
#' @param findOptimal A logical object. If TRUE, the value is estimated for
#'   all treatment options and that leading to the maximum value for each
#'   individual is used to estimate the value.
#'
#' @export
#' @name predict
#' @aliases predict,ITRSurv-method
#' @returns a list object containing a matrix of the predicted survival/cif function
#'   (Func), the estimated mean survival/cif (mean), and the estimated
#'   survival/cif probability (if critical value is mean.prob.combo or prob)
#' @examples
#'
#'
#' dt <- data.frame("Y.1" = sample(1:100,100,TRUE), "Y.2" = sample(101:200,100,TRUE),
#'                  "D.1" = rbinom(100, 1, 0.9), "D.2" = rbinom(100,1,0.9),
#'                  "A.1" = rbinom(100, 1, 0.5), "A.2" = rbinom(100,1,0.5),
#'                  "X.1" = rnorm(100), "X.2" = rnorm(100))
#'
#' result <- itrSurv(data = dt,
#'                   txName = c("A.1", "A.2"),
#'                   models = list(Surv(Y.1,D.1)~X.1+A.1,
#'                                 Surv(Y.2,D.2)~X.2+A.2+Y.1))
#'
#' tt <- predict(object = result)
#' tt <- predict(object = result, Phase = 1)
#' tt <- predict(object = result, findOptimal = FALSE)
#' tt <- predict(object = result, newdata = dt)
#' tt <- predict(object = result, newdata = dt, Phase = 1)
#' tt <- predict(object = result, newdaata = dt, findOptimal = FALSE)

setMethod(f = "predict",
          signature = c(object = "ITRSurv"),
          definition = function(object,
                                ...,
                                newdata,
                                Phase, # phase = 1 means predict survival; phase = 2 means predict CIF
                                findOptimal = TRUE) {
            # message("predict function from class_ITRSurv.R: Line 325")
            ob_itrsurv <<- object
            if (Phase > length(x = object@phaseResults)) {
              stop("requested Phase not present in analysis", call. = FALSE)
            }
            if (Phase == 1){
              params = object@params@survivalparam
            } else if (Phase == 2){
              params = object@params@endpointparam
            } else{
              stop("requested Phase not in analysis")
            }
            if (missing(x = newdata)) {
              print("class_ITRSurv.R: missing newdata")
              return( .Predict(Phase = Phase,
                               object = object@phaseResults[[ Phase ]],
                               newdata = NULL,
                               params = params,
                               findOptimal = findOptimal) )
            } else {
              # print(sprintf("Predicting for Phase %s", Phase))
              # print(dim(newdata))

              if (is.language(object@call[['tol1']])){
                # print(object@call[['tol1']])
                tol1 = c(object@call[["tol1"]][[2]],
                            object@call[["tol1"]][[3]])
              } else{
                tol1 <<- object@call[['tol1']]
              }
              # message("below calls class_survRF.R")
              # View(object)

              if (object@call[["pooled"]] == TRUE){
                # Phase, eps0, object, ..., newdata, model, txLevels, txName, params
                # Phase, eps0, object, ..., newdata, model, params, txLevels
                return( .Predict(object = object@phaseResults[[ Phase ]],
                                 newdata = newdata,
                                 Phase = Phase,
                                 txName = object@call[["txName"]],
                                 params = params,
                                 eps0 = tol1,
                                 findOptimal = findOptimal) )
              } else{
                return( .Predict(object = object@phaseResults[[ Phase ]],
                                 newdata = newdata,
                                 Phase = Phase,
                                 params = params,
                                 eps0 = tol1,
                                 findOptimal = findOptimal))
              }
           }
          })

