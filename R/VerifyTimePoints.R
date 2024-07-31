# Verify inputs 'timePointsPhase','timePoints', 'tau', and 'nTimes'
#
# methods are not exported and are for internal convenience only
#
# ensures that 'timePointsPhase' is numeric,
# ensures that 'timePoints' is numeric or character, generates time points if appropriate,
# ensures that 'nTimes' is appropriate.
#
# if both 'timePointsPhase' and 'timePoints'/'nTimes' are provided, then only 'timePointsPhase' is used
#
# successful methods return a vector of unique time points that are sorted in
#   ascending order and the maximum time point, tau
#
#
setGeneric(name = ".VerifyTimePoints",
           def = function(timePointsPhase, timePoints, nTimes, ...) {

                   standardGeneric(".VerifyTimePoints")
                 })

# the default method generates an error
setMethod(f = ".VerifyTimePoints",
          signature = c(timePointsPhase = "ANY",
                        # timePointsSurvival = "ANY",
                        # timePointsEndpoint = "ANY",
                        timePoints = "ANY",
                        nTimes = "ANY"),
          definition = function(timePointsPhase, timePoints, nTimes, ...) {
            stop("timePointsPhase must be a numeric vector",
                   # "timePointsEndpoint must be a numeric vector",
                   "if provided, timePoints input must be one of {'quad', 'uni', 'exp'} ",
                   "or a numeric vector but we strongly recommend using timePointsPhase",
                   call. = FALSE)
            })

# extends above for when timepointsphase is numeric
setMethod(f = ".VerifyTimePoints",
          signature = c(timePointsPhase = "numeric",
                        # timePointsSurvival = "numeric",
                        # timePointsEndpoint = "numeric",
                        timePoints = "ANY",
                        nTimes = "ANY"),
          definition = function(timePointsPhase, timePoints, nTimes, ..., tau, response) {

            # if (length(x = timePointsSurvival) == 0L) {
            #   stop("timePointsSurvival is of zero length", call. = FALSE)
            # }
            # if (length(x = timePointsEndpoint) == 0L) {
            #   stop("timePointsEndpoint is of zero length", call. = FALSE)
            # }
            # if timePoints provided, sort them and ensure uniqueness of values
            # timePointsSurvival <- sort(x = unique(x = timePointsSurvival))
            # timePointsEndpoint <- sort(x = unique(x = timePointsEndpoint))
            # ensures generated time points include 0 and maximum tau
            # if (min(timePointsSurvival) > 1e-8) {
              # timePointsSurvival <- c(0.0, timePointsSurvival)
            # }
            # if (min(timePointsEndpoint) > 1e-8) {
            #   timePointsEndpoint <- c(0.0, timePointsEndpoint)
            # }
            # if (is.null(tau)) {
            #   stop("tau should not be NULL. You can set it to the maximum of observed data for CR, or maximum of STOP column for RE.")
            # } else if (tau > max(timePointsSurvival) | tau > max(timePointsEndpoint)) {
            #   if (tau > max(timePointsSurvival)){
            #     timePointsPhase <- c(timePointsSurvival, tau)
            #   }
            #   if (tau > max(timePointsEndpoint)){
            #     timePointsEndpoint <- c(timePointsEndpoint, tau)
            #   }
            # } else if ((tau < max(timePointsSurvival) & tau > min(timePointsSurvival)) |
            #            (tau < max(timePointsEndpoint) & tau > min(timePointsEndpoint))) {
            #   if (tau < max(timePointsSurvival) & tau > min(timePointsSurvival)) {
            #     timePointsSurvival <- timePointsSurvival[timePointsSurvival < tau]
            #     timePointsSurvival <- c(timePointsSurvival, tau)
            #   }
            #   if (tau < max(timePointsEndpoint) & tau > min(timePointsEndpoint)) {
            #     timePointsEndpoint <- timePointsEndpoint[timePointsEndpoint < tau]
            #     timePointsEndpoint <- c(timePointsEndpoint, tau)
            #   }
            # } else if (tau < min(timePointsSurvival)){
            #   stop("tau cannot be less than all survival failure times", call. = FALSE)
            # } else if (tau < min(timePointsEndpoint)) {
            #   stop("tau cannot be less than all endpoint times", call. = FALSE)
            # }

            if (length(x = timePointsPhase) == 0L) {
              stop("timePointsPhase is of zero length", call. = FALSE)
            }

            # if timePointsPhase provided, sort them and ensure uniqueness of values
            timePointsPhase <- sort(x = unique(x = timePointsPhase))

            # ensures generated time points include 0 and maximum tau
            if (min(timePointsPhase) > 1e-8) {
              timePointsPhase <- c(0.0, timePointsPhase)
            }

            if (is.null(x = tau)) {
              stop("tau should not be NULL. You can set it to the maximum of observed data for CR, or maximum of STOP column for RE.")
            } else if (tau > max(timePointsPhase)) {
              timePointsPhase <- c(timePointsPhase, tau)
            } else if (tau < max(timePointsPhase) && tau > min(timePointsPhase)) {
              timePointsPhase <- timePointsPhase[timePointsPhase < tau]
              timePointsPhase <- c(timePointsPhase, tau)
            } else if (tau < min(timePointsPhase)) {
              stop("tau cannot be < all time points", call. = FALSE)
            }
            # message('assigning tau: ', tau)

            return( list("timePoints" = timePointsPhase,
                         # "timePointsSurvival" = timePointsSurvival,
                         # "timePointsEndpoint" = timePointsEndpoint,
                         "tau" = tau) )
          })

# extends above to prioritize timepointsurvival and timepointendpoint
setMethod(f = ".VerifyTimePoints",
          signature = c(timePointsPhase = "numeric",
                        # timePointsSurvival = "numeric",
                        # timePointsEndpoint = "numeric",
                        timePoints = "character",
                        nTimes = "numeric"),
          definition = function(timePointsPhase, timePoints, nTimes, ..., tau, response) {

            # if (length(x = timePointsSurvival) == 0L) {
            #   stop("timePointsSurvival is of zero length", call. = FALSE)
            # }
            # if (length(x = timePointsEndpoint) == 0L) {
            #   stop("timePointsEndpoint is of zero length", call. = FALSE)
            # }
            #
            # # if timePoints provided, sort them and ensure uniqueness of values
            # timePointsSurvival <- sort(x = unique(x = timePointsSurvival))
            # timePointsEndpoint <- sort(x = unique(x = timePointsEndpoint))
            #
            # # ensures generated time points include 0 and maximum tau
            # if (min(timePointsSurvival) > 1e-8) {
            #   timePointsSurvival <- c(0.0, timePointsSurvival)
            # }
            # if (min(timePointsEndpoint) > 1e-8) {
            #   timePointsEndpoint <- c(0.0, timePointsEndpoint)
            # }
            #
            # if (is.null(tau)) {
            #   stop("tau should not be NULL. You can set it to the maximum of observed data for CR, or maximum of STOP column for RE.")
            # } else if (tau > max(timePointsSurvival) | tau > max(timePointsEndpoint)) {
            #   if (tau > max(timePointsSurvival)){
            #     timePointsSurvival <- c(timePointsSurvival, tau)
            #   }
            #   if (tau > max(timePointsEndpoint)){
            #     timePointsEndpoint <- c(timePointsEndpoint, tau)
            #   }
            # } else if ((tau < max(timePointsSurvival) & tau > min(timePointsSurvival)) |
            #            (tau < max(timePointsEndpoint) & tau > min(timePointsEndpoint))) {
            #   if (tau < max(timePointsSurvival) & tau > min(timePointsSurvival)) {
            #     timePointsSurvival <- timePointsSurvival[timePointsSurvival < tau]
            #     timePointsSurvival <- c(timePointsSurvival, tau)
            #   }
            #   if (tau < max(timePointsEndpoint) & tau > min(timePointsEndpoint)) {
            #     timePointsEndpoint <- timePointsEndpoint[timePointsEndpoint < tau]
            #     timePointsEndpoint <- c(timePointsEndpoint, tau)
            #   }
            # } else if (tau < min(timePointsSurvival)){
            #   stop("tau cannot be less than all survival failure times", call. = FALSE)
            # } else if (tau < min(timePointsEndpoint)) {
            #   stop("tau cannot be less than all endpoint times", call. = FALSE)
            # }

            if (length(x = timePointsPhase) == 0L) {
              stop("timePointsPhase is of zero length", call. = FALSE)
            }

            # if timePointsPhase provided, sort them and ensure uniqueness of values
            timePointsPhase <- sort(x = unique(x = timePointsPhase))

            # ensures generated time points include 0 and maximum tau
            if (min(timePointsPhase) > 1e-8) {
              timePointsPhase <- c(0.0, timePointsPhase)
            }

            if (is.null(x = tau)) {
              stop("tau should not be NULL. You can set it to the maximum of observed data for CR, or maximum of STOP column for RE.")
            } else if (tau > max(timePointsPhase)) {
              timePointsPhase <- c(timePointsPhase, tau)
            } else if (tau < max(timePointsPhase) && tau > min(timePointsPhase)) {
              timePointsPhase <- timePointsPhase[timePointsPhase < tau]
              timePointsPhase <- c(timePointsPhase, tau)
            } else if (tau < min(timePointsPhase)) {
              stop("tau cannot be < all time points", call. = FALSE)
            }

            return( list("timePoints" = timePointsPhase,
                         # "timePointsEndpoint" = timePointsEndpoint,
                         "tau" = tau) )
          })

# extends above to prioritize timepointsurvival and timepointendpoint
setMethod(f = ".VerifyTimePoints",
          signature = c(timePointsPhase = "numeric",
                        # timePointsSurvival = "numeric",
                        # timePointsEndpoint = "numeric",
                        timePoints = "numeric",
                        nTimes = "ANY"),
          definition = function(timePointsPhase, timePoints, nTimes, ..., tau, response) {

            # if (length(x = timePointsSurvival) == 0L) {
            #   stop("timePointsSurvival is of zero length", call. = FALSE)
            # }
            # if (length(x = timePointsEndpoint) == 0L) {
            #   stop("timePointsEndpoint is of zero length", call. = FALSE)
            # }
            #
            # # if timePoints provided, sort them and ensure uniqueness of values
            # timePointsSurvival <- sort(x = unique(x = timePointsSurvival))
            # timePointsEndpoint <- sort(x = unique(x = timePointsEndpoint))
            #
            # # ensures generated time points include 0 and maximum tau
            # if (min(timePointsSurvival) > 1e-8) {
            #   timePointsSurvival <- c(0.0, timePointsSurvival)
            # }
            # if (min(timePointsEndpoint) > 1e-8) {
            #   timePointsEndpoint <- c(0.0, timePointsEndpoint)
            # }
            #
            # if (is.null(tau)) {
            #   stop("tau should not be NULL. You can set it to the maximum of observed data for CR, or maximum of STOP column for RE.")
            # } else if (tau > max(timePointsSurvival) | tau > max(timePointsEndpoint)) {
            #   if (tau > max(timePointsSurvival)){
            #     timePointsSurvival <- c(timePointsSurvival, tau)
            #   }
            #   if (tau > max(timePointsEndpoint)){
            #     timePointsEndpoint <- c(timePointsEndpoint, tau)
            #   }
            # } else if ((tau < max(timePointsSurvival) & tau > min(timePointsSurvival)) |
            #            (tau < max(timePointsEndpoint) & tau > min(timePointsEndpoint))) {
            #   if (tau < max(timePointsSurvival) & tau > min(timePointsSurvival)) {
            #     timePointsSurvival <- timePointsSurvival[timePointsSurvival < tau]
            #     timePointsSurvival <- c(timePointsSurvival, tau)
            #   }
            #   if (tau < max(timePointsEndpoint) & tau > min(timePointsEndpoint)) {
            #     timePointsEndpoint <- timePointsEndpoint[timePointsEndpoint < tau]
            #     timePointsEndpoint <- c(timePointsEndpoint, tau)
            #   }
            # } else if (tau < min(timePointsSurvival)){
            #   stop("tau cannot be less than all survival failure times", call. = FALSE)
            # } else if (tau < min(timePointsEndpoint)) {
            #   stop("tau cannot be less than all endpoint times", call. = FALSE)
            # }

            if (length(x = timePointsPhase) == 0L) {
              stop("timePointsPhase is of zero length", call. = FALSE)
            }

            # if timePointsPhase provided, sort them and ensure uniqueness of values
            timePointsPhase <- sort(x = unique(x = timePointsPhase))

            # ensures generated time points include 0 and maximum tau
            if (min(timePointsPhase) > 1e-8) {
              timePointsPhase <- c(0.0, timePointsPhase)
            }

            if (is.null(x = tau)) {
              stop("tau should not be NULL. You can set it to the maximum of observed data for CR, or maximum of STOP column for RE.")
            } else if (tau > max(timePointsPhase)) {
              timePointsPhase <- c(timePointsPhase, tau)
            } else if (tau < max(timePointsPhase) && tau > min(timePointsPhase)) {
              timePointsPhase <- timePointsPhase[timePointsPhase < tau]
              timePointsPhase <- c(timePointsPhase, tau)
            } else if (tau < min(timePointsPhase)) {
              stop("tau cannot be < all time points", call. = FALSE)
            }
            # message('assigning tau: ', tau)

            return( list("timePoints" = timePointsPhase,
                         # "timePointsEndpoint" = timePointsEndpoint,
                         "tau" = tau) )
          })


# time points is character and the time points are then generated from the distribution used to obtain time points from response data
setMethod(f = ".VerifyTimePoints",
          signature = c(timePointsPhase = "ANY",
                        # timePointsSurvival = "ANY",
                        # timePointsEndpoint = "ANY",
                        timePoints = "character",
                        nTimes = "numeric"),
          definition = function(timePointsPhase, timePoints, nTimes, ..., tau, response) {

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
              return( list("timePoints" = timePoints,
                           # "timePointsEndpoint" = timePoints,
                           "tau" = tau) )

            })

# extends above for when time points is numeric
setMethod(f = ".VerifyTimePoints",
          signature = c(timePointsPhase = "ANY",
                        # timePointsSurvival = "ANY",
                        # timePointsEndpoint = "ANY",
                        timePoints = "numeric",
                        nTimes = "ANY"),
          definition = function(timePointsPhase, timePoints, nTimes, ..., tau, response) {

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

              return( list("timePoints" = timePoints,
                           # "timePointsEndpoint" = timePoints,
                           "tau" = tau) )
            })


