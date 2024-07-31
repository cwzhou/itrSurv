# Verify input 'endpointTime'
#
# methods are not exported and are for internal convenience only
#
# ensures that 'endpointTime' is provided if criticalValue is {"prob", "mean.prob.combo"'}.
#
# successful methods return the numeric endpointTime or NULL.
#
setGeneric(name = ".VerifyEndpointTime",
           def = function(endpointTime, ...) {
             standardGeneric(".VerifyEndpointTime")
           })

# the default method generates an error
setMethod(f = ".VerifyEndpointTime",
          signature = c(endpointTime = "ANY"),
          definition = function(endpointTime, ...) {
            stop("evalTime must be numeric or NULL",
                 call. = FALSE)
          })

setMethod(f = ".VerifyEndpointTime",
          signature = c(endpointTime = "numeric"),
          definition = function(endpointTime, ..., criticalValue, tau) {

            if (!{criticalValue %in% c("prob", "mean.prob.combo")}) {
              message("evalTime is ignored if critical value is mean")
              return( NULL )
            }

            if (length(x = endpointTime) > 1L) {
              stop("only 1 value for evalTime can be given",
                   call. = FALSE)
            }

            if (endpointTime <= 0.0 || endpointTime > tau) {
              stop("evalTime must be between 0 and tau", call. = FALSE)
            }

            message("evalTime ", endpointTime)

            return( endpointTime )
          })

setMethod(f = ".VerifyEndpointTime",
          signature = c(endpointTime = "NULL"),
          definition = function(endpointTime, ..., tau) {
            message("NULL evalTime so we use: ", .VerifyEndpointTime(endpointTime = tau/2.0,
                                                                     tau = tau, ... ) )
            return( .VerifyEndpointTime(endpointTime = tau/2.0,
                                        tau = tau, ... ) )
          })
