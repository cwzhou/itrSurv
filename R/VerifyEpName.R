# Verify input 'epName'
#
# method is not exported and is for internal convenience only
#
# ensures that 'epName' is provided as a character or character vector and
# that the provided names are present in data.
# This is used to subset the Phase 1 survival data from the recurrent event data.
# This parameter is ignored for CR endpoint.
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
          definition = function(epName, ..., data) {

            if (length(x = epName) == 0L) {
              stop("epName must be provided", call. = FALSE)
            }

            # if recurrent event indicator name is in data, it is *the* recurrent event indicator name
            # if it is not, then we set an error because user must input if doing RE endpoint.

            test <- tryCatch(expr = data[,epName,drop = FALSE],
                             error = function(e) { return( NULL ) })

            if (is.null(x = test) && length(x = epName) == 1L) {
              stop("For RE endpoints, you must input 'epName', the recurrent event indicator for the dataset.")
            }

            test <- tryCatch(expr = data[,epName,drop = FALSE],
                             error = function(e) {
                               stop("unable to retrieve 'epName' from data",
                                    e$message, call. = FALSE)
                             })

            if (any(sapply(X = test, FUN = is.nan))) {
              stop("epName cannot include NaN values", call. = FALSE)
            }

            # ensure tx is factor or integer-like
            for (i in 1L:ncol(x = test)) {
              if (!is.factor(x = test[,i])) {
                if (is.numeric(x = test[,i])) {
                  if (!isTRUE(all.equal(target = test[,i],
                                        current = round(x = test[,i], digits = 0L)))) {

                    stop("EpName: recurrent event indicator variable must be integer or factor",
                         call. = FALSE)
                  }
                } else {
                  stop("EpName: recurrent event indicator variable must be integer or factor",
                       call. = FALSE)
                }
              }
            }

            message("IMPORTANT: RECURRENT EVENT INDICATOR MUST = 1 FOR RECURRENT EVENTS, AND = 0 FOR NON-RECURRENT EVENTS")

            return( epName )
          })
