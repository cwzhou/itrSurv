# Class extends CriticalValue to indicate that critical value is survival
#
# Class is not exported and is for internal convenience only
#
# @slot SurvivalTime A numeric object. The time at which the survival
#   probability is to be estimated
#
# @slot sIndex An integer object. The index of the timePoint vector above
#   which the survival time lies (and it is below the sIndex + 1 element)
#
# @slot sFraction A numeric object. The fractional location of SurvivalTime in
#   t[sIndex] and t[sIndex+1]
#
# @slot type A character object. Indicates of mean is to be used to break ties.
#
# Methods
#  .CriticalValueCriterion(object, ...) {not allowed}
#  .CreateValueObject(object, ...) {defined}
#  .IsSurvival(object, ...) {new; defined}
#
# Functions
# .criticalValueSurv(SurvivalTime, timePoints)
#
#' @include class_CriticalValue.R
#'
setClass(Class = "CriticalValueSurv",
         slots = c(survivalTime = "ANY",
                   sIndex = "integer",
                   sFraction = "numeric",
                   type = "character"),
         contains = c("CriticalValueBase"))

#-------------------------------------------------------------------------------
# method to identify if critical value is a mean or a probability
#-------------------------------------------------------------------------------
# method returns a character (specifically "mean")
#-------------------------------------------------------------------------------
setMethod(f = ".CriticalValueCriterion",
          signature = c(object = "CriticalValueSurv"),
          definition = function(object, ...) {

            # message("class_CriticalValueSurv.R: LINE 40")
            # message("sfdkljksld: object@type:", object@type)
              if (object@type == "surv.mean") return( "mean.prob.combo" )
              if (object@type == "surv.prob") return( "prob" )
              if (object@type == "surv.area") return( "area" )
            })

setMethod(f = "initialize",
         signature = c(.Object = "CriticalValueSurv"),
         def = function(.Object, ..., survivalTime, sIndex, sFraction, type) {

           # message("~~~ begin initialize ~~~ 1")

                   obj <- list(...)
                   tst <- sapply(X = obj,
                                 FUN = function(x){
                                         is(object = x,
                                            class2 = "CriticalValueSurv")
                                        })
                   # message("~~~ begin initialize ~~~ 2")

                   if (any(tst)) {
                     .Object <- obj[[ which(tst) ]]
                   } else if (missing(x = survivalTime) && missing(x = sIndex) &&
                         missing(x = sFraction) && missing(x = type)) {
                     .Object@survivalTime <- Inf
                     .Object@sIndex <- -1L
                     .Object@sFraction <- 0
                     .Object@type <- "none"
                   } else {
                     # message("~~~ begin initialize ~~~ 3")

                     if (missing(x = survivalTime) || missing(x = sIndex) ||
                         missing(x = sFraction) || missing(x = type)) {
                         gn <- unlist(lapply(list(...),is))
                         if( "CriticalValueSurv" %in% gn ) return(.Object)
                         stop("insufficient inputs provided")
                     }
                     if (type %in% c("prob")) {
                       .Object@type <- "surv.prob"
                     } else if (type %in% c("area")) {
                       .Object@type <- "surv.area"
                     } else if (type %in% c("mean","mean.prob.combo")) {
                       .Object@type <- "surv.mean"
                     }
                     # message("~~~ begin initialize ~~~ 4")

                     .Object@survivalTime <- survivalTime
                     .Object@sIndex <- sIndex
                     .Object@sFraction <- sFraction
                   }
                   # message("~~~ begin initialize ~~~ 5")
                   # print(.Object)
                   # print("!!!")
                   return( .Object )
                 })

#-------------------------------------------------------------------------------
# method to identify if critical value is of survival type
#-------------------------------------------------------------------------------
# method returns a logical
#-------------------------------------------------------------------------------
setGeneric(name = ".IsSurvival",
           def = function(object, ...) { standardGeneric(".IsSurvival") })

setMethod(f = ".IsSurvival",
          signature = c(object = "ANY"),
          definition = function(object, ...) { return( FALSE ) })

setMethod(f = ".IsSurvival",
          signature = c(object = "CriticalValueSurv"),
          definition = function(object, ...) { return( TRUE ) })

#-------------------------------------------------------------------------------
# internal function to create an object of class CriticalValueSurv
#-------------------------------------------------------------------------------
# function returns a CriticalValueSurv object
#-------------------------------------------------------------------------------
.criticalValueSurv <- function(survivalTime, timePoints, type) {

  # message("class_CriticalValueSurv.R")
  # message("-------starting .criticalValueSurv ---------")

  nTimes <- length(x = timePoints)
  # message("nTimes: ", nTimes)

  # index of last time point <= SurvivalTime
  sIndex <- sum(timePoints <= survivalTime)
  # message("sIndex: ", sIndex)

  if (sIndex < nTimes) {
    # message("if SurvivalTime is below tau, determine the fraction")
    sFraction <- {survivalTime - timePoints[sIndex]} /
      {timePoints[sIndex + 1L] - timePoints[sIndex]}
    # message("sFraction: ", sFraction)
  } else if (sIndex == 0L) {
    # message("if it is below the minimum, stop -- cannot extrapolate")
    stop("survival time is below minimum timepoint")
  } else {
    # message("if it is above tau, use tau -- cannot extrapolate")
    sFraction <- 0.0
    survivalTime <- max(timePoints)
    message("evalTime > tau --> survival time reset to tau")
  }

  return( new("CriticalValueSurv",
              "survivalTime" = survivalTime,
              "sIndex" = sIndex,
              "sFraction" = sFraction,
              "type" = type) )

}

