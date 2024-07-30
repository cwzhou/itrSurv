# Class to store information regarding time points
#
# Class is not exported and is for internal convenience only
#
#  @slot timePointsPhase A numeric vector; the time points used for analysis Phase.
#  Takes priority over timePoints/nTimes if available, as seen in VerifyTimePoints.R
#
#  @slot timePoints A numeric vector; the timepoints used in the analysis
#
#  @slot timeDiff A numeric vector; the time differences between the timepoints
#
#  @slot tau A numeric object; maximum time
#
# Getters
#   .TimePoints(object, ...) {new; defined}
#   .NTimes(object, ...) {new; defined}
#   .Tau(object, ...) {new; defined}
#   .TimeDiff(object, ...) {new; defined}
#
# Methods
#   .TimeInfoAsList(object, ...) {new; defined}
#
# Functions
#   .timeInfo(timePoints, nTimes, response)
#
# setClass(Class = "TimeInfo",
#          slots = c(timePointsSurvival = "numeric",
#                    timePointsEndpoint = "numeric",
#                    timeDiffSurvival = "numeric",
#                    timeDiffEndpoint = "numeric",
#                    tau = "numeric"))
# setClass(Class = "TimeInfoSurvival",
#          slots = c(timePoints = "numeric",
#                    timeDiff = "numeric",
#                    tau = "numeric"))
#
# setClass(Class = "TimeInfoEndpoint",
#          slots = c(timePoints = "numeric",
#                    timeDiff = "numeric",
#                    tau = "numeric"))
# setClass("TimeInfo",
#          representation=representation(TimeInfoSurvival="TimeInfoSurvival",
#                                        TimeInfoEndpoint="TimeInfoEndpoint"))

setClass(Class = "TimeInfo",
         slots = c(timePoints = "numeric",
                   timeDiff = "numeric",
                   tau = "numeric"))

### Getters

# #-------------------------------------------------------------------------------
# # Method to retrieve timepointsSurvival (Phase 1)
# #-------------------------------------------------------------------------------
# # Method returns the vector of timepoints used for the Phase 1 analysis
# #-------------------------------------------------------------------------------
# setGeneric(name = ".TimePointsSurvival",
#            def = function(object, ...) { standardGeneric(".TimePointsSurvival") })
#
# setMethod(f = ".TimePointsSurvival",
#           signature = c(object = "ANY"),
#           definition = function(object, ...) { stop("not allowed") })
#
# setMethod(f = ".TimePointsSurvival",
#           signature = c(object = "TimeInfoSurvival"),
#           definition = function(object, ...) { return( object@timePointsSurvival ) })
#
# #-------------------------------------------------------------------------------
# # Method to retrieve the number of timepointsSurvival
# #-------------------------------------------------------------------------------
# # Method returns a numeric
# #-------------------------------------------------------------------------------
# setGeneric(name = ".NTimesSurvival",
#            def = function(object, ...) { standardGeneric(".NTimesSurvival") })
#
# setMethod(f = ".NTimesSurvival",
#           signature = c(object = "ANY"),
#           definition = function(object, ...) { stop("not allowed") })
#
# setMethod(f = ".NTimesSurvival",
#           signature = c(object = "TimeInfoSurvival"),
#           definition = function(object, ...) {
#             return( length(x = object@timePointsSurvival) )
#           })
#
# #-------------------------------------------------------------------------------
# # Method to retrieve the difference in timepointsSurvival
# #-------------------------------------------------------------------------------
# # Method returns the vector of the differences in timepoints used for the analysis
# #-------------------------------------------------------------------------------
# setGeneric(name = ".TimeDiffSurvival",
#            def = function(object, ...) { standardGeneric(".TimeDiffSurvival") })
#
# setMethod(f = ".TimeDiffSurvival",
#           signature = c(object = "ANY"),
#           definition = function(object, ...) { stop("not allowed") })
#
# setMethod(f = ".TimeDiffSurvival",
#           signature = c(object = "TimeInfoSurvival"),
#           definition = function(object, ...) { return( object@timeDiffSurvival ) })
#
# #-------------------------------------------------------------------------------
# #-------------------------------------------------------------------------------
# #-------------------------------------------------------------------------------
#
# #-------------------------------------------------------------------------------
# # Method to retrieve timepointsEndpoint (Phase 2)
# #-------------------------------------------------------------------------------
# # Method returns the vector of timepoints used for the Phase 2 analysis
# #-------------------------------------------------------------------------------
# setGeneric(name = ".TimePointsEndpoint",
#            def = function(object, ...) { standardGeneric(".TimePointsEndpoint") })
#
# setMethod(f = ".TimePointsEndpoint",
#           signature = c(object = "ANY"),
#           definition = function(object, ...) { stop("not allowed") })
#
# setMethod(f = ".TimePointsEndpoint",
#           signature = c(object = "TimeInfoEndpoint"),
#           definition = function(object, ...) { return( object@timePointsEndpoint ) })
#
# #-------------------------------------------------------------------------------
# # Method to retrieve the number of timepointsEndpoint
# #-------------------------------------------------------------------------------
# # Method returns a numeric
# #-------------------------------------------------------------------------------
# setGeneric(name = ".NTimesEndpoint",
#            def = function(object, ...) { standardGeneric(".NTimesEndpoint") })
#
# setMethod(f = ".NTimesEndpoint",
#           signature = c(object = "ANY"),
#           definition = function(object, ...) { stop("not allowed") })
#
# setMethod(f = ".NTimesEndpoint",
#           signature = c(object = "TimeInfoEndpoint"),
#           definition = function(object, ...) {
#             return( length(x = object@timePointsEndpoint) )
#           })
#
# #-------------------------------------------------------------------------------
# # Method to retrieve the difference in timepointsEndpoint
# #-------------------------------------------------------------------------------
# # Method returns the vector of the differences in timepoints used for the analysis
# #-------------------------------------------------------------------------------
# setGeneric(name = ".TimeDiffEndpoint",
#            def = function(object, ...) { standardGeneric(".TimeDiffEndpoint") })
#
# setMethod(f = ".TimeDiffEndpoint",
#           signature = c(object = "ANY"),
#           definition = function(object, ...) { stop("not allowed") })
#
# setMethod(f = ".TimeDiffEndpoint",
#           signature = c(object = "TimeInfoEndpoint"),
#           definition = function(object, ...) { return( object@timeDiffEndpoint ) })

#-------------------------------------------------------------------------------
# Method to retrieve timepoints
#-------------------------------------------------------------------------------
# Method returns the vector of timepoints used for the analysis
#-------------------------------------------------------------------------------
setGeneric(name = ".TimePoints",
           def = function(object, ...) { standardGeneric(".TimePoints") })

setMethod(f = ".TimePoints",
          signature = c(object = "ANY"),
          definition = function(object, ...) { stop("not allowed") })

setMethod(f = ".TimePoints",
          signature = c(object = "TimeInfo"),
          definition = function(object, ...) { return( object@timePoints ) })

# setMethod(f = ".TimePoints",
#           signature = c(object = "TimeInfoSurvival"),
#           definition = function(object, ...) { return( object@timePoints ) })
#
# setMethod(f = ".TimePoints",
#           signature = c(object = "TimeInfoEndpoint"),
#           definition = function(object, ...) { return( object@timePoints ) })

#-------------------------------------------------------------------------------
# Method to retrieve the number of timepoints
#-------------------------------------------------------------------------------
# Method returns a numeric
#-------------------------------------------------------------------------------
setGeneric(name = ".NTimes",
           def = function(object, ...) { standardGeneric(".NTimes") })

setMethod(f = ".NTimes",
          signature = c(object = "ANY"),
          definition = function(object, ...) { stop("not allowed") })

setMethod(f = ".NTimes",
          signature = c(object = "TimeInfo"),
          definition = function(object, ...) {
            return( length(x = object@timePoints) )
          })

# setMethod(f = ".NTimes",
#           signature = c(object = "TimeInfoSurvival"),
#           definition = function(object, ...) {
#             return( length(x = object@timePoints) )
#           })
#
# setMethod(f = ".NTimes",
#           signature = c(object = "TimeInfoEndpoint"),
#           definition = function(object, ...) {
#             return( length(x = object@timePoints) )
#           })

#-------------------------------------------------------------------------------
# Method to retrieve the difference in timepoints
#-------------------------------------------------------------------------------
# Method returns the vector of the differences in timepoints used for the analysis
#-------------------------------------------------------------------------------
setGeneric(name = ".TimeDiff",
           def = function(object, ...) { standardGeneric(".TimeDiff") })

setMethod(f = ".TimeDiff",
          signature = c(object = "ANY"),
          definition = function(object, ...) { stop("not allowed") })

setMethod(f = ".TimeDiff",
          signature = c(object = "TimeInfo"),
          definition = function(object, ...) { return( object@timeDiff ) })

# setMethod(f = ".TimeDiff",
#           signature = c(object = "TimeInfoSurvival"),
#           definition = function(object, ...) { return( object@timeDiff ) })
#
# setMethod(f = ".TimeDiff",
#           signature = c(object = "TimeInfoEndpoint"),
#           definition = function(object, ...) { return( object@timeDiff ) })

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Method to retrieve the maximum timepoint
#-------------------------------------------------------------------------------
# Method returns a numeric
#-------------------------------------------------------------------------------
setGeneric(name = ".Tau",
           def = function(object, ...) { standardGeneric(".Tau") })

setMethod(f = ".Tau",
          signature = c(object = "ANY"),
          definition = function(object, ...) { stop("not allowed") })

setMethod(f = ".Tau",
          signature = c(object = "TimeInfo"),
          definition = function(object, ...) {

            # below is only if we do "representation" which we are not anymore.
            # if (object@TimeInfoEndpoint@tau != object@TimeInfoSurvival@tau){
            #   stop("Endpoint and Survival Tau's are not the same. This is not correct.")
            # }
            # return(object@TimeInfoSurvival@tau ) })

            return( object@tau ) })

#-------------------------------------------------------------------------------
# Method to retrieve timepoints for printing
#-------------------------------------------------------------------------------
# Method returns a list containing 1 element "timePoints" containing the
#  vector of timepoints used for the analysis
#-------------------------------------------------------------------------------
setGeneric(name = ".TimeInfoAsList",
           def = function(object, ...) {
                   standardGeneric(".TimeInfoAsList")
                 })

setMethod(f = ".TimeInfoAsList",
          signature = c(object = "ANY"),
          definition = function(object, ...) { stop("not allowed") })

setMethod(f = ".TimeInfoAsList",
          signature = c(object = "TimeInfo"),
          definition = function(object, ...) {
            return( list("timePoints" = object@timePoints) )
          })

# setMethod(f = ".TimeInfoAsList",
#           signature = c(object = "TimeInfoSurvival"),
#           definition = function(object, ...) {
#               return( list("timePoints" = object@timePoints) )
#             })
#
# setMethod(f = ".TimeInfoAsList",
#           signature = c(object = "TimeInfoEndpoint"),
#           definition = function(object, ...) {
#             return( list("timePoints" = object@timePoints) )
#           })

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Function to verify inputs and create a TimeInfo object
#-------------------------------------------------------------------------------
# Function returns a TimeInfo object
#-------------------------------------------------------------------------------

# recurrent events:

#' @include VerifyTimePoints.R
.timeInfo <- function(timePointsPhase,
                      # timePointsSurvival, timePointsEndpoint,
                      timePoints, nTimes, tau){ #, response) {

  # ensure that timePoints and nTimes are appropriate. Methods return 2 vectors
  # of unique time points that are sorted in ascending order.
  # First is for observed failure times if timePointsSurvival was inputted
  # Second is for observed endpoint times (observed recurrent event times if RE)
  # For CIF, they are the same because subsets.
  timePoints0 <- .VerifyTimePoints(timePointsPhase = timePointsPhase,
                                  # timePointsSurvival = timePointsSurvival,
                                  # timePointsEndpoint = timePointsEndpoint,
                                  timePoints = timePoints,
                                  tau = tau,
                                  nTimes = nTimes#,
                                  # response = response
                                  )

  tau <- timePoints0$tau
  timePoints <- timePoints0$timePoints

  # the total number of times points
  nTimes <- length(x = timePoints)

  timeDiff <- timePoints[-1L] - timePoints[-nTimes]
  timeDiff[nTimes] <- 0.0

  # deltaT; T_{i} - T_{i-1} i = 1:nTimes with T_0 = 0
  # timeDiff <- c(timePoints[-1L] - timePoints[-nTimes], 0)
  # this is used for the truncated mean calculation in which it is
  # assumed that all t > tau equals tau to t(nTimes+1)-t(nTimes) = 0

  # timePointsSurvival <- timePoints0$timePointsSurvival
  # timePointsEndpoint <- timePoints0$timePointsEndpoint
  # nTimesSurvival <- length(x = timePointsSurvival) #times points for Phase 1
  # nTimesEndpoint <- length(x = timePointsEndpoint) # time points for Phase 2
  # timeDiffSurvival <- timePointsSurvival[-1L] - timePointsSurvival[-nTimesSurvival]
  # timeDiffSurvival[nTimesSurvival] <- 0.0
  #
  # timeDiffEndpoint <- timePointsEndpoint[-1L] - timePointsEndpoint[-nTimesEndpoint]
  # timeDiffEndpoint[nTimesEndpoint] <- 0.0

  # survtime = new(Class = "TimeInfoSurvival",
  #                "tau" = tau,
  #                "timePoints" = timePointsSurvival,
  #                "timeDiff" = timeDiffSurvival
  #                )
  #
  # endpointtime = new(Class = "TimeInfoEndpoint",
  #               "tau" = tau,
  #               "timePoints" = timePointsEndpoint,
  #               "timeDiff" = timeDiffEndpoint)
  #
  # time = new("TimeInfo",
  #            TimeInfoSurvival=survtime,
  #            TimeInfoEndpoint = endpointtime)

  # new(Class = "TimeInfo",
  #     "tau" = tau,
  #     "timePointsSurvival" = timePointsSurvival,
  #     "timePointsEndpoint" = timePointsEndpoint,
  #     "timeDiffSurvival" = timeDiffSurvival,
  #     "timeDiffEndpoint" = timeDiffEndpoint)

  time = new(Class = "TimeInfo",
             "tau" = tau,
             "timePoints" = timePoints,
             "timeDiff" = timeDiff)

  return( time )
}

