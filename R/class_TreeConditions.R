# Class to store parameters that regulate tree and specify analysis preferences
#
# Class is not exported and is for internal convenience only
#
#  @slot nTree An integer object. The number of trees to be generated in forest.
#
#  @slot nodeSizeEnd An integer object. The minimum number of cases allowed in a
#    node

#  @slot nodeSizeSurv An integer object. The minimum number of cases allowed in a
#    node
#
#  @slot minEventEnd An integer object. The minimum number of events allowed in a
#    node

#  @slot minEventSurv An integer object. The minimum number of events allowed in a
#    node
#
#  @slot pooled A logical object. TRUE indicates that treatment groups are to
#    be considered in combination.
#
#  @slot stratifiedSplit A number object. The coefficient phi for stratified
#    random spliting
#
# Getters
#   .NTree(object, ...) {new; defined}
#   .NodeSize(object, ...) {new; defined}
#   .MinEvent(object, ...) {new; defined}
#   .Pooled(object, ...) {new; defined}
#   .StratifiedSplit(object, ...) {new; defined}
#
# Methods
#   .TreeConditionsAsList(object, ...) {new; defined}
#
# Functions
# .treeConditions(..., nTree, nodeSizeEnd, nodeSizeSurv, minEventEnd, minEventSurv,
#                 pooled, stratifiedSplit)
#
setClass(Class = "TreeConditions",
         slots = c(nTree = "integer",
                   nodeSizeEnd = "integer",
                   nodeSizeSurv = "integer",
                   minEventEnd = "integer",
                   minEventSurv = "integer",
                   pooled = "logical",
                   stratifiedSplit = "numeric"))

## Getters

#-------------------------------------------------------------------------------
# Method to retrieve number of trees in the forest
#-------------------------------------------------------------------------------
# Method returns an integer
#-------------------------------------------------------------------------------
setGeneric(name = ".NTree",
           def = function(object, ...) {
                   standardGeneric(".NTree")
                 })

setMethod(f = ".NTree",
          signature = c(object = "ANY"),
          definition = function(object, ...) { stop("not allowed") })

setMethod(f = ".NTree",
          signature = c(object = "TreeConditions"),
          definition = function(object, ...) { return( object@nTree ) })

#-------------------------------------------------------------------------------
# Method to retrieve minimum number of endpoint cases allowed in a node
#-------------------------------------------------------------------------------
# Method returns an integer
#-------------------------------------------------------------------------------
setGeneric(name = ".NodeSizeEnd",
           def = function(object, ...) { standardGeneric(".NodeSizeEnd") })

setMethod(f = ".NodeSizeEnd",
          signature = c(object = "ANY"),
          definition = function(object, ...) { stop("not allowed") })

setMethod(f = ".NodeSizeEnd",
          signature = c(object = "TreeConditions"),
          definition = function(object, ...) { return( object@nodeSizeEnd ) })

#-------------------------------------------------------------------------------
# Method to retrieve minimum number of phase 1 cases allowed in a node
#-------------------------------------------------------------------------------
# Method returns an integer
#-------------------------------------------------------------------------------
setGeneric(name = ".NodeSizeSurv",
           def = function(object, ...) { standardGeneric(".NodeSizeSurv") })

setMethod(f = ".NodeSizeSurv",
          signature = c(object = "ANY"),
          definition = function(object, ...) { stop("not allowed") })

setMethod(f = ".NodeSizeSurv",
          signature = c(object = "TreeConditions"),
          definition = function(object, ...) { return( object@nodeSizeSurv ) })

#-------------------------------------------------------------------------------
# Method to retrieve minimum number of Phase 2 events allowed in a node
#-------------------------------------------------------------------------------
# Method returns an integer
#-------------------------------------------------------------------------------
setGeneric(name = ".MinEventEnd",
           def = function(object, ...) { standardGeneric(".MinEventEnd") })

setMethod(f = ".MinEventEnd",
          signature = c(object = "ANY"),
          definition = function(object, ...) { stop("not allowed") })

setMethod(f = ".MinEventEnd",
          signature = c(object = "TreeConditions"),
          definition = function(object, ...) { return( object@minEventEnd ) })

#-------------------------------------------------------------------------------
# Method to retrieve minimum number of Phase 1 events allowed in a node
#-------------------------------------------------------------------------------
# Method returns an integer
#-------------------------------------------------------------------------------
setGeneric(name = ".MinEventSurv",
           def = function(object, ...) { standardGeneric(".MinEventSurv") })

setMethod(f = ".MinEventSurv",
          signature = c(object = "ANY"),
          definition = function(object, ...) { stop("not allowed") })

setMethod(f = ".MinEventSurv",
          signature = c(object = "TreeConditions"),
          definition = function(object, ...) { return( object@minEventSurv ) })


#-------------------------------------------------------------------------------
# Method to retrieve flag for pooled analysis
#-------------------------------------------------------------------------------
# Method returns a logical
#-------------------------------------------------------------------------------
setGeneric(name = ".Pooled",
           def = function(object, ...) { standardGeneric(".Pooled") })

setMethod(f = ".Pooled",
          signature = c(object = "ANY"),
          definition = function(object, ...) { stop("not allowed") })

setMethod(f = ".Pooled",
          signature = c(object = "TreeConditions"),
          definition = function(object, ...) { return( object@pooled ) })

#-------------------------------------------------------------------------------
# Method to retrieve coefficient for stratified analysis
#-------------------------------------------------------------------------------
# Method returns a numeric
#-------------------------------------------------------------------------------
setGeneric(name = ".StratifiedSplit",
           def = function(object, ...) { standardGeneric(".StratifiedSplit") })

setMethod(f = ".StratifiedSplit",
          signature = c(object = "ANY"),
          definition = function(object, ...) { stop("not allowed") })

setMethod(f = ".StratifiedSplit",
          signature = c(object = "TreeConditions"),
          definition = function(object, ...) { return( object@stratifiedSplit ) })

#-------------------------------------------------------------------------------
# Method to retrieve primary slots for printing
#-------------------------------------------------------------------------------
# Method returns a list containing 6 elements
#   "nTree" an integer, the total number of trees in forest
#   "nodeSizeEnd" an integer, the minimum number of Phase 2 cases in a node
#   "nodeSizeSurv" an integer, the minimum number of Phase 1 cases in a node
#   "minEventEnd" an integer, the minimum number of Phase 2 events in a node
#   "minEventSurv" an integer, the minimum number of Phase 1 events in a node
#   "pooled" a logical indicating of treatments were pooled
#   "stratifiedSplit" the coefficient phi for stratified
#-------------------------------------------------------------------------------
setGeneric(name = ".TreeConditionsAsList",
           def = function(object, ...) {
                   standardGeneric(".TreeConditionsAsList")
                 })

setMethod(f = ".TreeConditionsAsList",
          signature = c(object = "ANY"),
          definition = function(object, ...) { stop("not allowed") })

setMethod(f = ".TreeConditionsAsList",
          signature = c(object = "TreeConditions"),
          definition = function(object, ...) {
              return( list("nTree" = object@nTree,
                           "nodeSizeEnd" = object@nodeSizeEnd,
                           "nodeSizeSurv" = object@nodeSizeSurv,
                           "minEventEnd" = object@minEventEnd,
                           "minEventSurv" = object@minEventSurv,
                           "pooled" = object@pooled,
                           "stratifiedSplit" = object@stratifiedSplit) )
            })

.treeConditions <- function(...,
                            nTree,
                            nodeSizeEnd,
                            nodeSizeSurv,
                            minEventEnd,
                            minEventSurv,
                            pooled,
                            stratifiedSplit) {
  # Validate and adjust the values of input parameters

  # Ensure minEventEnd is a positive integer
  # minimum number of End events must be integer and > 0
  if (!is.numeric(x = minEventEnd)) stop("minEventEnd must be integer", call. = FALSE)
  minEventEnd <- as.integer(x = minEventEnd)
  if (minEventEnd < 1L) stop("minEventEnd must be non-zero positive", call. = FALSE)

  # Ensure minEventSurv is a positive integer
  # minimum number of Surv events must be integer and > 0
  if (!is.numeric(x = minEventSurv)) stop("minEventSurv must be integer", call. = FALSE)
  minEventSurv <- as.integer(x = minEventSurv)
  if (minEventSurv < 1L) stop("minEventSurv must be non-zero positive", call. = FALSE)

  # Ensure nodeSizeEnd is a positive integer
  # minimum number of End cases in each node, must be integer and > 0
  if (!is.numeric(x = nodeSizeEnd)) stop("nodeSizeEnd must be integer", call. = FALSE)
  nodeSizeEnd <- as.integer(x = nodeSizeEnd)
  if (nodeSizeEnd < 1L) stop("nodeSizeEnd must be non-zero positive", call. = FALSE)

  # Ensure nodeSizeSurv is a positive integer
  # minimum number of Surv cases in each node, must be integer and > 0
  if (!is.numeric(x = nodeSizeSurv)) stop("nodeSizeSurv must be integer", call. = FALSE)
  nodeSizeSurv <- as.integer(x = nodeSizeSurv)
  if (nodeSizeSurv < 1L) stop("nodeSizeSurv must be non-zero positive", call. = FALSE)

  # Ensure nTree is a positive integer
  # number of trees to grow in forest, must be integer and > 0
  if (!is.numeric(x = nTree)) stop("nTree must be integer", call. = FALSE)
  nTree <- as.integer(x = nTree)
  if (nTree < 1L) stop("nTree must be non-zero positive", call. = FALSE)

  # Ensure pooled is a logical value
  if (!is.logical(x = pooled) || is.na(x = pooled)) {
     stop("pooled must be logical", call. = FALSE)
  }

  # Adjust and validate the value of stratifiedSplit
  if (is.null(x = stratifiedSplit) || stratifiedSplit <= 1e-8) {
    stratifiedSplit <- 0.0
  } else {
    if (stratifiedSplit > 1.0) {
      stop("stratifiedSplit must be [0,1]", call. = FALSE)
    }
  }

  # Create and return a new TreeConditions object with the validated parameters
  return( new(Class = "TreeConditions",
              "nTree" = nTree,
              "nodeSizeSurv" = nodeSizeSurv,
              "nodeSizeEnd" = nodeSizeEnd,
              "minEventSurv" = minEventSurv,
              "minEventEnd" = minEventEnd,
              "pooled" = pooled,
              "stratifiedSplit" = stratifiedSplit) )

}
