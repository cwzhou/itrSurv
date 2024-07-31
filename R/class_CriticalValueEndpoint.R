# Class extends CriticalValue to indicate that critical value is for endpoint
#
# Class is not exported and is for internal convenience only
#
# @slot endpointTime A numeric object. The time at which the Endpoint
#   probability is to be estimated
#
# @slot cIndex An integer object. The index of the timePoint vector above
#   which the Endpoint time lies (and it is below the cIndex + 1 element)
#
# @slot cFraction A numeric object. The fractional location of endpointTime in
#   t[cIndex] and t[cIndex+1]
#
# @slot type A character object. Indicates of mean is to be used to break ties.
#
# Methods
#  .CriticalValueCriterion(object, ...) {not allowed}
#  .CreateValueObject(object, ...) {defined}
#  .IsEndpoint(object, ...) {new; defined}
#
# Functions
# .criticalValueCR(endpointTime, timePoints)
#
#' @include class_CriticalValue.R
setClass(Class = "CriticalValueEndpoint",
         slots = c(endpointTime = "ANY",
                   cIndex = "integer",
                   cFraction = "numeric",
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
            if (object@type == "cif.mean") return( "mean.prob.combo" )
            if (object@type == "cif.prob") return( "prob" )
            if (object@type == "cif.area") return( "area" )
          })

setMethod(f = "initialize",
          signature = c(.Object = "CriticalValueEndpoint"),
          def = function(.Object, ..., endpointTime, cIndex, cFraction, type) {

            obj <- list(...)
            tst <- sapply(X = obj,
                          FUN = function(x){
                            is(object = x,
                               class2 = "CriticalValueEndpoint")
                          })

            if (any(tst)) {
              # print("testing")
              # print(which(tst))
              # print("obj")
              # print(obj)
              .Object <- obj[[ which(tst) ]]
            } else if (missing(x = endpointTime) && missing(x = cIndex) &&
                       missing(x = cFraction) && missing(x = type)) {
              .Object@endpointTime <- Inf
              .Object@cIndex <- -1L
              .Object@cFraction <- 0
              .Object@type <- "none"
            } else {
              if (missing(x = endpointTime) || missing(x = cIndex) ||
                  missing(x = cFraction) || missing(x = type)) {
                gn <- unlist(lapply(list(...),is))
                if( "CriticalValueEndpoint" %in% gn ) return(.Object)
                stop("insufficient inputs provided")
              }
              if (type %in% c("prob")) {
                .Object@type <- "cif.prob"
              } else if (type %in% c("area")) {
                .Object@type <- "cif.area"
              } else if (type %in% c("mean","mean.prob.combo")) {
                .Object@type <- "cif.mean"
              }
              .Object@endpointTime <- endpointTime
              .Object@cIndex <- cIndex
              .Object@cFraction <- cFraction
            }
            return( .Object )
          })

#-------------------------------------------------------------------------------
# method to identify if critical value is of Endpoint type
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
.criticalValueCR <- function(endpointTime, timePoints, type) {

  # message("-------starting .criticalValueCR ---------")
  # print(type)
  nTimes <- length(x = timePoints)

  # index of last time point <= endpointTime
  cIndex <- sum(timePoints <= endpointTime)

  if (cIndex < nTimes) {
    # if endpointTime is below tau, determine the fraction
    cFraction <- {endpointTime - timePoints[cIndex]} /
      {timePoints[cIndex + 1L] - timePoints[cIndex]}
  } else if (cIndex == 0L) {
    # if it is below the minimum, stop -- cannot extrapolate
    stop("Endpoint time is below minimum timepoint")
  } else {
    # if it is above tau, use tau -- cannot extrapolate
    cFraction <- 0.0
    endpointTime <- max(timePoints)
    message("Endpoint time reset to tau")
  }
  # print(type)
  # print("endpointTime")
  # print(endpointTime)
  return( new("CriticalValueEndpoint",
              "endpointTime" = endpointTime,
              "cIndex" = cIndex,
              "cFraction" = cFraction,
              "type" = type) )

}

