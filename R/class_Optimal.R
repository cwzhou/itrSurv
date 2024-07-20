# Class for storing estimated optimal treatment and value
#
# Class is not exported and is for internal convenience only
#
#  @slot optimalTx A vector object. The index of the estimated optimal tx
#
#  @slot optimalY A vector. The estimated value
#
#  @slot type A character. One of "mean" or "prob"
#
# Getters
#  .OptimalY(object, ...) {new; defined}
#  .OptimalAsList(object, ...) {new; defined}
#
# Define the Optimal class with 3 slots (attributes)
setClass(Class = "Optimal",
         slots = c("optimalTx" = "vector", #vector object to store the index of the estimated optimal treatment
                   "optimalY" = "matrix", #matrix for storing estimated value
                   "type" = "character", #character that can have value "mean" or "prob"
                   "NumTrts" = "vector", # number of possible optimal treatments based on survival Phase 1
                   # "NonOpt_Opt_Ratio" = "vector",
                   "Ratio_Stopping_Ind" = "vector")) # indicator to stop at P1 ( = 1) or continue to P2 (= 0).
# if NumTrts > 1 then we go to P2.

#-------------------------------------------------------------------------------
# Method to retrieve the estimated optimal value
#-------------------------------------------------------------------------------
# Method returns a vector object
#-------------------------------------------------------------------------------
setGeneric(name = ".OptimalY", # define generic function to retrieve estimated optimal value from "Optimal" object
           def = function(object, ...) { standardGeneric(".OptimalY") })
# catch-all method that throws an error if you call this function on an object of any class other than "Optimal"
setMethod(f = ".OptimalY",
          signature = c(object = "ANY"),
          definition = function(object, ...) { stop("not allowed") })
# method to retrieve "optimalY" (matrix of estimated optimal values) slot/attribute from "Optimal" object
setMethod(f = ".OptimalY",
          signature = c(object = "Optimal"),
          definition = function(object, ...) { return( object@optimalY ) })

#-------------------------------------------------------------------------------
# Method to return an Optimal object as a list (used for printing)
#-------------------------------------------------------------------------------
# Method returns a list object
#-------------------------------------------------------------------------------
setGeneric(name = ".OptimalAsList", # define generic function to return an "Optimal" object as a list (could be useful for printing or extracting information)
           def = function(object, ...) { standardGeneric(".OptimalAsList") })
# catch-all method that throws error if calling function on an object of any class other than "Optimal"
setMethod(f = ".OptimalAsList",
          signature = c(object = "ANY"),
          definition = function(object, ...) { stop("not allowed") })
# method returns a list containing the three slots of "Optimal" object
setMethod(f = ".OptimalAsList",
          signature = c(object = "Optimal"),
          definition = function(object, ...) {
              return( list("optimalTx" = object@optimalTx,
                           "optimalY" = object@optimalY,
                           "type" = object@type,
                           # "NonOpt_Opt_Ratio" = object@Ratio,
                           "NumTrts" = object@NumTrts,
                           "Ratio_Stopping_Ind" = object@Ratio_Stopping_Ind) )
            })
