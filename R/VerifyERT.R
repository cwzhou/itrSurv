# Verify input 'ERT'
#
# methods are not exported and are for internal convenience only
#
# ensures that 'ERT' is provided as a logical or NULL object.
# NOTE: this is only valid for RSF and RCIF for coding right now.
# NOT coded up for RMMF in Fortran as of Nov 2024 due to multiple records per person
# Nov 2024: Set ERT = FALSE if endPoint = RE and ERT is null. Added error for
#
# successful methods return a logical indicating if the
# Extremely Randomized Tree method is to be used
#
setGeneric(name = ".VerifyERT",
           def = function(ERT, ..., endPoint) { standardGeneric(".VerifyERT") })

# the default method generates an error
setMethod(f = ".VerifyERT",
          signature = c(ERT = "ANY"),
          definition = function(ERT, ..., endPoint) {
              stop("ERT must be logical or NULL. Not coded yet for RE.", call. = FALSE)
            })

setMethod(f = ".VerifyERT",
          signature = c(ERT = "logical"),
          definition = function(ERT, ..., endPoint) {
              if (endPoint == "RE" && ERT == TRUE) {
                stop("ERT not coded up yet for Endpoint: RE")
              }
              if (is.na(x = ERT)) {
                stop("ERT must be logical or NULL", call. = FALSE)
              }
              return( ERT )
            })

setMethod(f = ".VerifyERT",
          signature = c(ERT = "NULL"),
          definition = function(ERT, ..., endPoint) {
              if (endPoint == "RE") {
                print("ERT not coded up yet for Endpoint: RE")
                return ( FALSE )
              }
              return( TRUE )
            })
