# Class extends CriticalValue to indicate that critical value is CIF
#
# Class is not exported and is for internal convenience only
#
# @slot CIFTime A numeric object. The time at which the CIF
#   probability is to be estimated
#
# @slot cIndex An integer object. The index of the timePoint vector above
#   which the CIF time lies (and it is below the cIndex + 1 element)
#
# @slot cFraction A numeric object. The fractional location of CIFTime in
#   t[cIndex] and t[cIndex+1]
#
# @slot type A character object. Indicates of mean is to be used to break ties.
#
# Methods
#  .CriticalValueCriterion(object, ...) {not allowed}
#  .CreateValueObject(object, ...) {defined}
#  .IsCIF(object, ...) {new; defined}
#
# Functions
# .criticalValueCR(CIFTime, timePoints)
#
#' @include class_CriticalValue.R
setClass(Class = "CriticalValueCR",
         slots = c(CIFTime = "ANY",
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
          signature = c(object = "CriticalValueCR"),
          definition = function(object, ...) {
            if (object@type == "cif.mean") return( "mean.prob.combo" )
            if (object@type == "cif.prob") return( "prob" )
            if (object@type == "cif.area") return( "area" )
          })

setMethod(f = "initialize",
          signature = c(.Object = "CriticalValueCR"),
          def = function(.Object, ..., CIFTime, cIndex, cFraction, type) {

            obj <- list(...)
            tst <- sapply(X = obj,
                          FUN = function(x){
                            is(object = x,
                               class2 = "CriticalValueCR")
                          })

            if (any(tst)) {
              # print("testing")
              # print(which(tst))
              # print("obj")
              # print(obj)
              .Object <- obj[[ which(tst) ]]
            } else if (missing(x = CIFTime) && missing(x = cIndex) &&
                       missing(x = cFraction) && missing(x = type)) {
              .Object@CIFTime <- Inf
              .Object@cIndex <- -1L
              .Object@cFraction <- 0
              .Object@type <- "none"
            } else {
              if (missing(x = CIFTime) || missing(x = cIndex) ||
                  missing(x = cFraction) || missing(x = type)) {
                gn <- unlist(lapply(list(...),is))
                if( "CriticalValueCR" %in% gn ) return(.Object)
                stop("insufficient inputs provided")
              }
              if (type %in% c("prob")) {
                .Object@type <- "cif.prob"
              } else if (type %in% c("area")) {
                .Object@type <- "cif.area"
              } else if (type %in% c("mean","mean.prob.combo")) {
                .Object@type <- "cif.mean"
              }
              .Object@CIFTime <- CIFTime
              .Object@cIndex <- cIndex
              .Object@cFraction <- cFraction
            }
            return( .Object )
          })

#-------------------------------------------------------------------------------
# method to identify if critical value is of CIF type
#-------------------------------------------------------------------------------
# method returns a logical
#-------------------------------------------------------------------------------
setGeneric(name = ".IsCIF",
           def = function(object, ...) { standardGeneric(".IsCIF") })

setMethod(f = ".IsCIF",
          signature = c(object = "ANY"),
          definition = function(object, ...) { return( FALSE ) })

setMethod(f = ".IsCIF",
          signature = c(object = "CriticalValueCR"),
          definition = function(object, ...) { return( TRUE ) })

#-------------------------------------------------------------------------------
# internal function to create an object of class CriticalValueCR
#-------------------------------------------------------------------------------
# function returns a CriticalValueCR object
#-------------------------------------------------------------------------------
.criticalValueCR <- function(CIFTime, timePoints, type) {

  # message("-------starting .criticalValueCR ---------")
  # print(type)
  nTimes <- length(x = timePoints)

  # index of last time point <= CIFTime
  cIndex <- sum(timePoints <= CIFTime)

  if (cIndex < nTimes) {
    # if CIFTime is below tau, determine the fraction
    cFraction <- {CIFTime - timePoints[cIndex]} /
      {timePoints[cIndex + 1L] - timePoints[cIndex]}
  } else if (cIndex == 0L) {
    # if it is below the minimum, stop -- cannot extrapolate
    stop("CIF time is below minimum timepoint")
  } else {
    # if it is above tau, use tau -- cannot extrapolate
    cFraction <- 0.0
    CIFTime <- max(timePoints)
    message("CIF time reset to tau")
  }
  # print(type)
  # print("CIFTime")
  # print(CIFTime)
  return( new("CriticalValueCR",
              "CIFTime" = CIFTime,
              "cIndex" = cIndex,
              "cFraction" = cFraction,
              "type" = type) )

}

