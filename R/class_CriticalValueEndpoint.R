# Class extends CriticalValue to indicate that critical value is endpoint
#
# Class is not exported and is for internal convenience only
#
# @slot EndpointTime A numeric object. The time at which the endpoint
#   probability is to be estimated
#
# @slot eIndex An integer object. The index of the timePoint vector above
#   which the endpoint time lies (and it is below the eIndex + 1 element)
#
# @slot eFraction A numeric object. The fractional location of EndpointTime in
#   t[eIndex] and t[eIndex+1]
#
# @slot type A character object. Indicates of mean is to be used to break ties.
#
# Methods
#  .CriticalValueCriterion(object, ...) {not allowed}
#  .CreateValueObject(object, ...) {defined}
#  .IsEndpoint(object, ...) {new; defined}
#
# Functions
# .criticalValueEndpoint(EndpointTime, timePoints)
#
#' @include class_CriticalValue.R
#'
setClass(Class = "CriticalValueEndpoint",
         slots = c(endpointTime = "ANY",
                   eIndex = "integer",
                   eFraction = "numeric",
                   type = "character"),
         contains = c("CriticalValueBase"))

#-------------------------------------------------------------------------------
# method to identify if critical value is a mean or a probability
#-------------------------------------------------------------------------------
# method returns a character (specifically "mean")
#-------------------------------------------------------------------------------
setMethod(f = ".CriticalValueCriterion",
          signature = c(object = "CriticalValueEndpoint"),
          definition = function(object, ...) {

            # message("class_CriticalValueEndpoint.R: LINE 41")
            # message("sfdkljksld: object@type:", object@type)
            if (object@type == "end.mean") return( "mean.prob.combo" )
            if (object@type == "end.prob") return( "prob" )
            if (object@type == "end.area") return( "area" )
          })

setMethod(f = "initialize",
          signature = c(.Object = "CriticalValueEndpoint"),
          def = function(.Object, ..., endpointTime, eIndex, eFraction, type) {

            # message("~~~ begin initialize ~~~ 1")

            obj <- list(...)
            tst <- sapply(X = obj,
                          FUN = function(x){
                            is(object = x,
                               class2 = "CriticalValueEndpoint")
                          })
            # message("~~~ begin initialize ~~~ 2")

            if (any(tst)) {
              .Object <- obj[[ which(tst) ]]
            } else if (missing(x = endpointTime) && missing(x = eIndex) &&
                       missing(x = eFraction) && missing(x = type)) {
              .Object@endpointTime <- Inf
              .Object@eIndex <- -1L
              .Object@eFraction <- 0
              .Object@type <- "none"
            } else {
              # message("~~~ begin initialize ~~~ 3")

              if (missing(x = endpointTime) || missing(x = eIndex) ||
                  missing(x = eFraction) || missing(x = type)) {
                gn <- unlist(lapply(list(...),is))
                if( "CriticalValueEndpoint" %in% gn ) return(.Object)
                stop("insufficient inputs provided")
              }
              if (type %in% c("prob")) {
                .Object@type <- "end.prob"
              } else if (type %in% c("area")) {
                .Object@type <- "end.area"
              } else if (type %in% c("mean","mean.prob.combo")) {
                .Object@type <- "end.mean"
              }
              # message("~~~ begin initialize ~~~ 4")

              .Object@endpointTime <- endpointTime
              .Object@eIndex <- eIndex
              .Object@eFraction <- eFraction
            }
            # message("~~~ begin initialize ~~~ 5")
            # print(.Object)
            # print("!!!")
            return( .Object )
          })

#-------------------------------------------------------------------------------
# method to identify if critical value is of endpoint type
#-------------------------------------------------------------------------------
# method returns a logical
#-------------------------------------------------------------------------------
setGeneric(name = ".IsEndpoint",
           def = function(object, ...) { standardGeneric(".IsEndpoint") })

setMethod(f = ".IsEndpoint",
          signature = c(object = "ANY"),
          definition = function(object, ...) { return( FALSE ) })

setMethod(f = ".IsEndpoint",
          signature = c(object = "CriticalValueEndpoint"),
          definition = function(object, ...) { return( TRUE ) })

#-------------------------------------------------------------------------------
# internal function to create an object of class CriticalValueEndpoint
#-------------------------------------------------------------------------------
# function returns a CriticalValueEndpoint object
#-------------------------------------------------------------------------------
.criticalValueEndpoint <- function(endpointTime, timePoints, type) {

  message("-------starting .criticalValueEndpoint ---------")

  nTimes <- length(x = timePoints)
  # message("nTimes: ", nTimes)

  # index of last time point <= EndpointTime
  eIndex <- sum(timePoints <= endpointTime)
  # message("eIndex: ", eIndex)

  if (eIndex < nTimes) {
    # message("if EndpointTime is below tau, determine the fraction")
    eFraction <- {endpointTime - timePoints[eIndex]} /
      {timePoints[eIndex + 1L] - timePoints[eIndex]}
    # message("eFraction: ", eFraction)
  } else if (eIndex == 0L) {
    # message("if it is below the minimum, stop -- cannot extrapolate")
    stop("endpoint time is below minimum timepoint")
  } else {
    # message("if it is above tau, use tau -- cannot extrapolate")
    eFraction <- 0.0
    endpointTime <- max(timePoints)
    message("evalTime > tau --> endpoint time reset to tau")
  }

  return( new("CriticalValueEndpoint",
              "endpointTime" = endpointTime,
              "eIndex" = eIndex,
              "eFraction" = eFraction,
              "type" = type) )

}

