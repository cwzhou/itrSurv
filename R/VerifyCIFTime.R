# Verify input 'CIFTime'
#
# methods are not exported and are for internal convenience only
#
# ensures that 'CIFTime' is provided if criticalValue is {"prob", "mean.prob.combo"'}.
#
# successful methods return the numeric CIFTime or NULL.
#
setGeneric(name = ".VerifyCIFTime",
           def = function(CIFTime, ...) {
             standardGeneric(".VerifyCIFTime")
           })

# the default method generates an error
setMethod(f = ".VerifyCIFTime",
          signature = c(CIFTime = "ANY"),
          definition = function(CIFTime, ...) {
            stop("evalTime must be numeric or NULL",
                 call. = FALSE)
          })

setMethod(f = ".VerifyCIFTime",
          signature = c(CIFTime = "numeric"),
          definition = function(CIFTime, ..., criticalValue, tau) {

            if (!{criticalValue %in% c("prob", "mean.prob.combo")}) {
              message("evalTime is ignored if critical value is mean")
              return( NULL )
            }

            if (length(x = CIFTime) > 1L) {
              stop("only 1 value for evalTime can be given",
                   call. = FALSE)
            }

            if (CIFTime <= 0.0 || CIFTime > tau) {
              stop("evalTime must be between 0 and tau", call. = FALSE)
            }

            message("evalTime ", CIFTime)

            return( CIFTime )
          })

setMethod(f = ".VerifyCIFTime",
          signature = c(CIFTime = "NULL"),
          definition = function(CIFTime, ..., tau) {
            message("NULL evalTime so we use: ", .VerifyCIFTime(CIFTime = tau/2.0,
                                                                     tau = tau, ... ) )
            return( .VerifyCIFTime(CIFTime = tau/2.0,
                                        tau = tau, ... ) )
          })
