# Virtual class to store information regarding critical value selection
#
# Class is not exported and is for internal convenience only
#
# Methods
#  .CriticalValueAsList(object, ...) {not allowed}
#  .CriticalValueCriterion(object, ...) {not allowed}
#  .CreateValueObject(object, ...) {not allowed}
#
setClass(Class = "CriticalValueBase", # defines the virtual class CriticalValueBase. Objects of this class cannot be directly instantiated but serve as a base for other classes.
         contains = c("VIRTUAL"))

#-------------------------------------------------------------------------------
# method to identify if critical value is a mean or a probability
#-------------------------------------------------------------------------------
# method is not defined for a general CriticalValue object
#-------------------------------------------------------------------------------
setGeneric(name = ".CriticalValueCriterion",
           def = function(object, ...) { standardGeneric(".CriticalValueCriterion") }) #not intended to be called directly on objects of the CriticalValueBase class. Instead, it's meant to be overridden in subclasses to identify if the critical value is a mean or a probability.

setMethod(f = ".CriticalValueCriterion", # .CriticalValueCriterion Default Method:
          signature = c(object = "ANY"),
          definition = function(object, ...) {

            print("class_CriticalValue.R")
            stop("class_CriticalValue.R: not allowed")

            }) # Generates an error message that it's not allowed to be called directly on objects of the CriticalValueBase class.



resample <- function(x, ...) x[sample.int(n = length(x = x), ...)] #Takes a vector x and performs a random resampling of its elements using sample.int.


# Defines a virtual class which serves as a base class for objects that store information related to critical value selection.
# This class contains methods that are not allowed to be directly called for any object of this class. Instead, these methods are meant to be overridden in subclasses of CriticalValueBase:
# In summary, the provided code establishes a virtual base class for objects related to critical value selection and includes a placeholder method for identifying the type of critical value (mean or probability).
# Subclasses are expected to override this method to provide specific implementations.
