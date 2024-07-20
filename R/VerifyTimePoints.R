# Verify inputs 'timePoints', 'tau', and 'nTimes'
#
# methods are not exported and are for internal convenience only
#
# ensures that 'timePoints' is numeric or character, generates time points if
# appropriate, ensures that 'nTimes' is appropriate.
#
# successful methods return a vector of unique time points that are sorted in
#   ascending order and the maximum time point, tau
#
setGeneric(name = ".VerifyTimePoints",
           def = function(timePoints, nTimes, ...) {
                   standardGeneric(".VerifyTimePoints")
                 })

# the default method generates an error
setMethod(f = ".VerifyTimePoints",
          signature = c(timePoints = "ANY",
                        nTimes = "ANY"),
          definition = function(timePoints, nTimes, ...) {
              stop("timePoints input must be one of {'quad', 'uni', 'exp'} ",
                   "or a numeric vector", call. = FALSE)
            })

# time points is character and the time points are then generated from the distribution used to obtain time points from response data
setMethod(f = ".VerifyTimePoints",
          signature = c(timePoints = "character",
                        nTimes = "numeric"),
          definition = function(timePoints, nTimes, ..., tau, response) {

              # if timePoints is provided as a character, the character indicates
              # which distribution should be used to obtain the time points from
              # the response data. Current options are {quad, uni, exp} for
              # quadratic, uniform, and exponential, respectively.

              # checks if "nTimes" is positive integer
              nTimes <- as.integer(x = nTimes)

              if (nTimes <= 0L) stop("nTimes must be positive", call. = FALSE)

              # if timePoints are generated internally, must use one
              # of "quad", "uni" or "exp"

              timePoints <- tolower(x = timePoints)
              # checks if "timePoints" is one of the values "quad," "uni," or "exp."
              # If not, it calls the generic method again with NULL arguments, effectively raising an error.
              if (!timePoints %in% c("quad", "uni", "exp")) {
                .VerifyTimePoints(timePoints = NULL, nTimes = NULL)
              }

              # identify the 75% percentile of "response" data (or "times" if "response" is matrix)
              if (is.matrix(x = response)) {
                times <- rowSums(x = response, na.rm = TRUE)
              } else {
                times <- response
              }

              if (is.null(x = tau)) {
                timeThirdQuartile <- stats::quantile(x = times, probs = 0.75)
                message("the study length (tau) was not provided ",
                        "and was set as the third quartile of the observed times: ", timeThirdQuartile)
              } else {
                timeThirdQuartile <- tau
              }

              if (timePoints == "quad") {

                # quadratically spaced times points between 0 and third quartile
                timePoints <- seq(from = 0.0,
                                  to = ceiling(x = sqrt(x = timeThirdQuartile)),
                                  length.out = nTimes)^2

              } else if (timePoints == "uni") {

                if (timeThirdQuartile > 1.0) {
                  # evenly spaced time points between 0 and ceiling third quartile
                  timePoints <- seq(from = 0.0,
                                    to = ceiling(x = timeThirdQuartile),
                                    length.out = nTimes)
                } else {
                  # evenly spaced time points between 0 and rounded third quartile
                  timePoints <- seq(from = 0.0,
                                    to = round(x = timeThirdQuartile, digits = 4L),
                                    length.out = nTimes)
                }

              } else if (timePoints == "exp") {

                if (timeThirdQuartile > 1.0) {
                  # exponentially spaced time points between 0 and third quartile
                  timePoints <- exp(x = seq(from = 0,
                                            to = log(x = ceiling(x = timeThirdQuartile) + 1L),
                                            length.out = nTimes)) - 1.0

                } else {
                  # evenly spaced time points between 0 and rounded third quartile
                  timePoints <- exp(x = seq(from = 0,
                                            to = log(x = round(x = timeThirdQuartile, digits = 4L) + 1L),
                                            length.out = nTimes))
                }
              }

              # ensures generated time points include 0 and maximum tau
              if (min(timePoints) > 1e-8) {
                timePoints <- c(0.0, timePoints)
              }

              if (is.null(x = tau)) {
                tau <- max(timePoints)
              } else if (tau > max(timePoints)) {
                timePoints <- c(timePoints, tau)
              } else if (tau < max(timePoints) && tau > min(timePoints)) {
                timePoints <- timePoints[timePoints < tau]
                timePoints <- c(timePoints, tau)
              } else if (tau < min(timePoints)) {
                stop("tau cannot be < all time points", call. = FALSE)
              }

              message('assigning tau: ', tau)

              # returns list containing generated time points and maximum time point "tau"
              return( list("timePoints" = timePoints, "tau" = tau) )

            })

# extends above for when time points is numeric
setMethod(f = ".VerifyTimePoints",
          signature = c(timePoints = "numeric",
                        nTimes = "ANY"),
          definition = function(timePoints, nTimes, ..., tau, response) {

              if (length(x = timePoints) == 0L) {
                stop("timePoints is of zero length", call. = FALSE)
              }

              # if timePoints provided, sort them and ensure uniqueness of values
              timePoints <- sort(x = unique(x = timePoints))

              # ensures generated time points include 0 and maximum tau
              if (min(timePoints) > 1e-8) {
                timePoints <- c(0.0, timePoints)
              }

              if (is.null(x = tau)) {
                tau <- max(timePoints)
              } else if (tau > max(timePoints)) {
                timePoints <- c(timePoints, tau)
              } else if (tau < max(timePoints) && tau > min(timePoints)) {
                timePoints <- timePoints[timePoints < tau]
                timePoints <- c(timePoints, tau)
              } else if (tau < min(timePoints)) {
                stop("tau cannot be < all time points", call. = FALSE)
              }
              # message('assigning tau: ', tau)

              return( list("timePoints" = timePoints, "tau" = tau) )
            })
