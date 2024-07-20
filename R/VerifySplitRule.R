# Verify input 'splitRule'
#
# methods are not exported and are for internal convenience only
#
# ensures that 'splitRule' is one of {'logrank', 'mean', 'logrankcr'}.
#
# successful methods return the original character object with possible
#  modification to make all lower case
#

setGeneric(name = ".VerifySplitRule", # Define a generic function with different methods for different input types
           def = function(splitRule, ...) {
                   standardGeneric(".VerifySplitRule")
                 })

# the default method generates an error
setMethod(f = ".VerifySplitRule",
          signature = c(splitRule = "ANY"), # Defines default method for ANY type that generates error if splitRule is not logrank or mean.
          definition = function(splitRule, ...) {
            print("default error")
              stop("splitRule must be one of {'logrank', 'mean', 'logrankcr'}", call. = FALSE) # logrankcr = gray's test
            })

setMethod(f = ".VerifySplitRule",
          signature = c(splitRule = "NULL"), # Defines method where splitRule is NULL input type.
          # Checks the value of endPoint and criticalValue and sets splitRule to either "mean" or "logrank" accordingly.
          definition = function(splitRule, ..., endPoint, criticalValue) {

            # 'mean', 'prob', 'area', 'mean.prob.combo'

              if (criticalValue == "mean" | criticalValue == "area") { # Basically sets a default value for splitRule based on the value of criticalValue.
                # step1: survival
                splitRule = "mean"
                # step2: endPoint
              } else{
                # step1: survival
                splitRule = "logrank"
                # step2: endPoint
                # if (endPoint == "CR"){
                #   splitRule = "logrankcr"
                # } else if (endPoint == "RE"){
                #   splitRule = "logrankRE"
                # } else if (endPoint == "MC"){
                #   splitRule = "logrankMC"
                # }
              }
            # print(splitRule)

              return( .VerifySplitRule(splitRule = splitRule, ...) ) # Recursively calls .VerifySplitRule to verify the value.
            })

setMethod(f = ".VerifySplitRule",
          signature = c(splitRule = "character"), # Defines method where splitRule is character input type.
          # Converts the splitRule to lowercase and checks if it is one of {'logrank', 'mean'}. If not, it generates an error message.
          definition = function(splitRule, ...) {

              splitRule <- tolower(x = splitRule)
              # print(splitRule)
              if (!(splitRule %in% c("logrank", "mean", "logrankcr", "meancr", "logrankRE", "meanRE", "logrankMC", "meanMC"))) {
                stop("splitRule must be one of {'logrank', 'mean', 'logrankcr', 'meancr', 'logrankRE', 'meanRE', 'logrankMC', 'meanMC'}")
              }
              return( splitRule )
            })

# In summary, this code ensures that the splitRule parameter is either "logrank" or "mean" and converts it to lowercase if needed.
# It provides default values for splitRule based on the value of criticalValue if splitRule is NULL.
# These methods are used to standardize and validate the input for the splitRule parameter within the package's internal functions.
