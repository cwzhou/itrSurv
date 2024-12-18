# Class extends CriticalValue to indicate that critical value is non-survival
#  area
#
# Class is not exported and is for internal convenience only
#
# Methods
#  .CriticalValueCriterion(object, ...) {defined}
#  .CreateValueObject(object, ...) {defined}
#
# Functions
#  .criticalValueMean(...)
#
#' @include class_CriticalValue.R
setClass(Class = "CriticalValueArea",
         contains = c("CriticalValueBase"))

#-------------------------------------------------------------------------------
# method to identify if critical value is a mean or a probability
#-------------------------------------------------------------------------------
# method returns a character (specifically "area")
#-------------------------------------------------------------------------------
setMethod(f = ".CriticalValueCriterion",
          signature = c(object = "CriticalValueArea"),
          definition = function(object, ...) {

            # print("class_CriticalValueArea.R: LINE 26")
            return("area")
          }
)

#-------------------------------------------------------------------------------
# internal function to create an object of class CriticalValueArea
#-------------------------------------------------------------------------------
# function returns a CriticalValueArea object
#-------------------------------------------------------------------------------
.criticalValueArea <- function(...) { return( new("CriticalValueArea") ) }
