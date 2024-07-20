# Verify input 'criticalValue'
#
# methods are not exported and are only for internal convenience
#
# ensures that 'criticalValue' is one of {'mean', 'prob', 'area', 'mean.prob.combo'}.
#
# successful methods return the original character possibly modified to be
# all lower case.
#
setGeneric(name = ".VerifyCriticalValue",
           def = function(criticalValue, ...) {
                   standardGeneric(".VerifyCriticalValue")
                 })

# the default method generates an error
setMethod(f = ".VerifyCriticalValue",
          signature = c(criticalValue = "ANY"),
          definition = function(criticalValue, ...) {
              stop("criticalValue must be one of ",
                   "{'mean', 'prob', 'area', 'mean.prob.combo'}",
                   call. = FALSE)
            })

setMethod(f = ".VerifyCriticalValue",
          signature = c(criticalValue = "character"),
          definition = function(criticalValue, ...) {

              criticalValue <- tolower(x = criticalValue)

              if (criticalValue %in% c('mean', 'prob', 'area', 'mean.prob.combo'))
                return( criticalValue )

              stop("criticalValue must be one of ",
                   "{'mean', 'prob', 'area', 'mean.prob.combo'}",
                   call. = FALSE)

            })
