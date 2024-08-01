# Verify input 'idName'
#
# method is not exported and is for internal convenience only
#
# ensures that 'idName' is provided as a character or character vector and
# that the provided names are present in data.
# This is used to subset the Phase 1 survival data from the recurrent event data.
# This parameter is ignored for CR endpoint.
# We can ignore decision points in the single stage setting.
#
# successful methods return the original input without modification.
#
setGeneric(name = ".VerifyIdName",
           def = function(idName, ...) { standardGeneric(".VerifyIdName") })

# the default method generates an error
setMethod(f = ".VerifyIdName",
          signature = c(idName = "ANY"),
          definition = function(idName, ...) {
            stop("idName must be a vector of character objects",
                 call. = FALSE)
          })

setMethod(f = ".VerifyIdName",
          signature = c(idName = "character"),
          definition = function(idName, ..., data) {

            if (length(x = idName) == 0L) {
              stop("idName must be provided", call. = FALSE)
            }

            # if id name is in data, it is *the* id name
            # if it is not, then we set an error because user must input if doing RE endpoint.

            test <- tryCatch(expr = data[,idName,drop = FALSE],
                             error = function(e) { return( NULL ) })

            if (is.null(x = test) && length(x = idName) == 1L) {
              stop("For RE endpoints, you must input 'idName', the variable that identifies subjects.")
            }

            test <- tryCatch(expr = data[,idName,drop = FALSE],
                             error = function(e) {
                               stop("unable to retrieve 'idName' from data",
                                    e$message, call. = FALSE)
                             })

            if (any(sapply(X = test, FUN = is.nan))) {
              stop("idName cannot include NaN values", call. = FALSE)
            }

            # ensure tx is factor or integer-like
            for (i in 1L:ncol(x = test)) {
              if (!is.factor(x = test[,i])) {
                if (is.numeric(x = test[,i])) {
                  if (!isTRUE(all.equal(target = test[,i],
                                        current = round(x = test[,i], digits = 0L)))) {

                    stop("id variable must be integer or factor",
                         call. = FALSE)
                  }
                } else {
                  stop("id variable must be integer or factor",
                       call. = FALSE)
                }
              }
            }

            return( idName )
          })
