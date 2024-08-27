# Verify input 'yName'
#
# method is not exported and is for internal convenience only
#
# ensures that 'yName' is provided as a character or character vector and
# that the provided names are present in data.
# This is used to make sure the dataset is sorted in ascending order by the observed failure time for CR endpoint.
# This parameter is ignored for RE endpoint.
# We can ignore decision points in the single stage setting.
#
# successful methods return the original input without modification.
#
setGeneric(name = ".VerifyYName",
           def = function(yName, ...) { standardGeneric(".VerifyYName") })

# the default method generates an error
setMethod(f = ".VerifyYName",
          signature = c(yName = "ANY"),
          definition = function(yName, ...) {
            stop("yName must be a vector of character objects",
                 call. = FALSE)
          })

setMethod(f = ".VerifyYName",
          signature = c(yName = "character"),
          definition = function(yName, ..., data) {

            if (length(x = yName) == 0L) {
              stop("yName must be provided", call. = FALSE)
            }

            # if observed time variable name is in data, it is *the* observed time variable name
            # if it is not, then we set an error because user must input for CR endpoint.

            test <- tryCatch(expr = data[,yName,drop = FALSE],
                             error = function(e) { return( NULL ) })

            if (is.null(x = test) && length(x = yName) == 1L) {
              stop("For CR endpoints, you must input 'yName', the observed failure time variable.")
            }

            test <- tryCatch(expr = data[,yName,drop = FALSE],
                             error = function(e) {
                               stop("unable to retrieve 'yName' from data",
                                    e$message, call. = FALSE)
                             })

            if (any(sapply(X = test, FUN = is.nan))) {
              stop("yName cannot include NaN values", call. = FALSE)
            }



            # ensure y is numeric
            for (i in 1L:ncol(x = test)) {
              if (!is.numeric(x = test[,i])) {
                stop("observed time variable must be numeric", call. = FALSE)
              }
            }

            # for (i in 1L:ncol(x = test)) {
            #   if (!is.factor(x = test[,i])) {
            #     if (!is.numeric(x = test[,i])) {
            #       if (!isTRUE(all.equal(target = test[,i],
            #                             current = round(x = test[,i], digits = 0L)))) {
            #
            #         stop("observed time variable variable must be NUMERIC",
            #              call. = FALSE)
            #       }
            #     } else {
            #       stop("observed time variable variable must be NUMERIC",
            #            call. = FALSE)
            #     }
            #   }
            # }

            return( yName )
          })
