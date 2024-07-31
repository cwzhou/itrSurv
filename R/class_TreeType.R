# Class to store information regarding tree type ERT or Breiman
#
# Class is not exported and is for internal convenience only
#
#  @slot replace A logical object; TRUE indicates that samples can include
#    duplicate records
#
#  @slot randomSplit A numeric object; The probability of a random split
#
#  @slot ERT A logical object; TRUE indicates that extremely randomized trees
#    methods are to be used
#
#  @slot uniformSplit A logical object; TRUE indicates that cutoffs are to
#    be selected based on a uniformed random number
#
#  @slot splitRule A character object; must be one of {'logrank', 'mean', 'logrankcr', 'logrankre'}
#
#  @slot tieMethod A character object; must be one of {'first', 'random', 'NA'}
#
#  @slot nSamples A numeric object; the sample size of the dataset
#
# Getters
#
# Methods
#
# Functions
#
#' @include VerifyERT.R VerifyRandomSplit.R
#' @include VerifyUniformSplit.R VerifyReplace.R VerifySplitRule.R
#' @include VerifyTieMethod.R

setClass(Class = "TreeType",
         slots = c("replace" = "logical",
                   "randomSplit" = "numeric",
                   "ERT" = "logical",
                   "uniformSplit" = "logical",
                   "splitRule" = "character",
                   "tieMethod" = "character",
                   "nSamples" = "numeric"))

## Getters

# initializer
.treeType <- function(endPoint,
                      ERT,
                      nSamples,
                      uniformSplit,
                      replace,
                      splitRule,
                      tieMethod,
                      randomSplit,
                      criticalValue) {

  # ensure that ERT is logical or NULL. Methods return a logical.
  ERT <- .VerifyERT(ERT = ERT)

  # ensure that randomSplit is 0 <= rs < 1. Methods return a numeric.
  randomSplit <- .VerifyRandomSplit(randomSplit = randomSplit)

  # ensure that uniformSplit is logical or NULL. Methods return a logical.
  uniformSplit <- .VerifyUniformSplit(uniformSplit = uniformSplit, ERT = ERT)

  # ensure that replace is logical or NULL. Methods return a logical.
  replace <- .VerifyReplace(replace = replace, ERT = ERT)

  # # if endPoint = "CR" then splitRule1 = splitRule (user input) and splitRule2 = "logrankcr" or "meancr"
  # splitRule1 = splitRule
  # if (endPoint == "CR"){
  #   splitRule2 = paste0(splitRule1, "CR")
  # } else if (endPoint == "RE"){
  #   splitRule2 = paste0(splitRule1, "RE")
  # } else if (endPoint == "MC"){
  #   splitRule2 = paste0(splitRule1, "MC")
  # } else{
  #   stop("need to define a splitrule for step2 endpoint.")
  # }
  # verify splitRule. methods return the original character object with possible
  # modification to make all lower case
  # step1 survival: splitRule1
  splitRule <- .VerifySplitRule(splitRule = splitRule,
                                criticalValue = criticalValue)
  # splitRule2 <- .VerifySplitRule(splitRule = splitRule2,
  #                                criticalValue = criticalValue2)

  # successful methods return the original character input possibly modified to
  # lower case
  tieMethod <- .VerifyTieMethod(tieMethod = tieMethod)

  return( new(Class = "TreeType",
              "replace" = replace,
              "randomSplit" = randomSplit,
              "ERT" = ERT,
              "uniformSplit" = uniformSplit,
              "splitRule" = splitRule,
              "tieMethod" = tieMethod,
              "nSamples" = nSamples) )

}
