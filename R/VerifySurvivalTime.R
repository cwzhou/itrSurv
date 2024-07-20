# Verify input 'survivalTime'
#
# methods are not exported and are for internal convenience only
#
# ensures that 'survivalTime' is provided if criticalValue is one of
# {'prob', 'mean.prob.combo'}.
#
# successful methods return the numeric survivalTime or NULL.
#
setGeneric(name = ".VerifySurvivalTime",
           def = function(survivalTime, ...) {
             # message("standardGeneric .VerifySurvivalTime")
                   standardGeneric(".VerifySurvivalTime")
                 })

# the default method generates an error
setMethod(f = ".VerifySurvivalTime",
          signature = c(survivalTime = "ANY"),
          definition = function(survivalTime, ...) {
            # message("default method to produce error.")
              stop("evalTime must be numeric or NULL",
                   call. = FALSE)
            })

setMethod(f = ".VerifySurvivalTime",
          signature = c(survivalTime = "numeric"),
          definition = function(survivalTime, ..., criticalValue, tau) {

            # message(".VerifySurvivalTime for numeric survival time")
            if (!{criticalValue %in% c("prob", "mean.prob.combo")}) {
              message("evalTime is ignored if critical value is mean")
              return( NULL )
            }
              if (length(x = survivalTime) > 1L) {
                stop("only 1 value for evalTime can be given",
                     call. = FALSE)
              }

              if (survivalTime <= 0.0 || survivalTime > tau) {
                stop("evalTime must be between 0 and tau", call. = FALSE)
              }

              message("evalTime ", survivalTime)

              return( survivalTime )
            })

setMethod(f = ".VerifySurvivalTime",
          signature = c(survivalTime = "NULL"),
          definition = function(survivalTime, ..., tau) {

            message("NULL evalTime so we use: ", .VerifySurvivalTime(survivalTime = tau/2.0,
                                                                     tau = tau, ... ) )
              return( .VerifySurvivalTime(survivalTime = tau/2.0,
                                          tau = tau, ... ) )
            })
