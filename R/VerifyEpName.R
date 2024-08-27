# Verify input 'epName'
#
# method is not exported and is for internal convenience only
#
# ensures that 'epName' is provided as a character or character vector and
# that the provided names are present in data.
# Has various purposes depending on endPoint.
# For "RE" endpoint: Refers to indicator if row is recurrent event (1) or not (0).
#   Later used to subset Phase 1 survival dataset (one row per subject).
#   This is used to subset the Phase 1 survival data from the recurrent event data.
# For "CR" endpoint: Refers to status indicator for causes: 0 = censored, 1 = priority cause, 2 = any other cause.
#'  Should be inputted as one row per person for CR data.
#'  Not needed if splitRule = csh_cr for endpoint 'CR'
# We can ignore decision points in the single stage setting.
#
# successful methods return the original input without modification.
#
setGeneric(name = ".VerifyEpName",
           def = function(epName, ...) { standardGeneric(".VerifyEpName") })

# the default method generates an error
setMethod(f = ".VerifyEpName",
          signature = c(epName = "ANY"),
          definition = function(epName, ...) {
            stop("epName must be a vector of character objects",
                 call. = FALSE)
          })

setMethod(f = ".VerifyEpName",
          signature = c(epName = "character"),
          definition = function(epName, ..., data, endPoint) {

              if (length(x = epName) == 0L) {
                stop("epName must be provided", call. = FALSE)
              }

              # epName = status (0,1,2) variable for CR (required for gray_cr test)
              # epName = recurrent event indicator for RE
              # if epName is in data, it is *the* variable name
              # if it is not, then we set an error because user must input for RE, and if rule2 = gray_cr or NULL for CR.

              test <- tryCatch(expr = data[,epName,drop = FALSE],
                               error = function(e) { return( NULL ) })

              if (is.null(x = test) && length(x = epName) == 1L) {
                if (endPoint == "RE"){
                  stop("For RE endpoints, you must input 'epName', the recurrent event indicator for the dataset.")
                } else{
                  stop("For CR endpoint using gray_cr test, you must input 'epName', the status indicator for dataset (0,1,2).")
                }
              }

              test <- tryCatch(expr = data[,epName,drop = FALSE],
                               error = function(e) {
                                 stop("unable to retrieve 'epName' from data",
                                      e$message, call. = FALSE)
                               })

              if (any(sapply(X = test, FUN = is.nan))) {
                stop("epName cannot include NaN values", call. = FALSE)
              }

              # ensure status is factor or integer-like
              for (i in 1L:ncol(x = test)) {
                if (!is.factor(x = test[,i])) {
                  if (is.numeric(x = test[,i])) {
                    if (!isTRUE(all.equal(target = test[,i],
                                          current = round(x = test[,i], digits = 0L)))) {

                      stop("EpName: variable must be integer or factor",
                           call. = FALSE)
                    }
                  } else {
                    stop("EpName: variable must be integer or factor",
                         call. = FALSE)
                  }
                }
              }

              if (endPoint == "RE"){
                message("IMPORTANT: RECURRENT EVENT INDICATOR MUST = 1 FOR RECURRENT EVENTS, AND = 0 FOR NON-RECURRENT EVENTS")
              }
            return( epName )
          })
