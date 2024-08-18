# Verify input 'splitRule'
#
# methods are not exported and are for internal convenience only
#
# ensures that 'splitRule' is one of \{'logrank_surv', 'mean_surv', 'gray_cr', 'csh_cr', 'gray_re'\}.
#
# successful methods return the original character object with possible
#  modification to make all lower case

#'   Must be one of \{'logrank_surv', 'mean_surv', 'gray_cr', 'csh_cr', 'gray_re'\}
#'   indicating the test used to determine an optimal split.
#'   If NULL and Phase = 'Survival' and 'criticalValue' = 'mean', then takes value 'mean'.
#'   If NULL and Phase = 'Survival' and 'criticalValue' = 'prob' or 'mean.prob.combo', then takes value 'logrank'.
#'   If NULL and Phase = 'CR' and endPoint = 'CR', takes value 'graycr'. # default over cshcr
#'   If NULL and Phase = 'RE' and endPoint = 'RE', takes value 'grayre'.
#'   ## survival endpoint ##
#'   logrank = 1
#'   mean = 2
#'   ## other endpoints (CR, RE) ##
#'   cr_gray = 3
#'   cr_csh = 4
#'   re_qlr = 5
#'   if splitRule <=2 then Phase 1 survival
#'   if splitRule 3-4 then Phase 2 endpoint (CR)
#'   if splitRule >4 then Phase 2 endpoint (RE)
#'
setGeneric(name = ".VerifySplitRule", # Define a generic function with different methods for different input types
           def = function(splitRule, ...) {
                   standardGeneric(".VerifySplitRule")
                 })

# Defines default method for ANY type that generates error if splitRule is not one of \{'logrank_surv', 'mean_surv', 'gray_cr', 'csh_cr', 'gray_re'\}.
# the default method generates an error
setMethod(f = ".VerifySplitRule",
          signature = c(splitRule = "ANY"),
          definition = function(splitRule, ...) {
            print("default error")
              stop("splitRule must be one of {'logrank_surv', 'mean_surv', 'gray_cr', 'csh_cr', 'gray_re'}", call. = FALSE) # logrankcr = gray's test
            })

setMethod(f = ".VerifySplitRule",
          signature = c(splitRule = "NULL"), # Defines method where splitRule is NULL input type.
          # Checks the value of endPoint and criticalValue and sets splitRule to either 'mean_surv' or 'logrank_surv' for Phase 1; or 'gray_cr','csh_cr','gray_re' for Phase2, accordingly.
          definition = function(splitRule, ..., Phase, endPoint, criticalValue) {

            # 'mean', 'prob', 'area', 'mean.prob.combo'
            #'   If NULL and Phase = 1 and 'criticalValue' = 'mean', then takes value 'mean'.
            #'   If NULL and Phase = 1 and 'criticalValue' = 'prob' or 'mean.prob.combo', then takes value 'logrank'.
            #'   If NULL and Phase = 2 and endPoint = 'CR', takes value 'gray_cr'. # default over csh_cr
            #'   If NULL and Phase = 2 and endPoint = 'RE', takes value 'gray_re'.

            message("splitRule is NULL -- read documentation for how defaults are set")
            if (Phase == 1){
              # Phase 1: Survival ('logrank_surv' or 'mean_surv')
              # Sets a default value for splitRule based on the value of criticalValue.
              if (criticalValue == "mean" | criticalValue == "area") {
                splitRule = "mean_surv"
              } else{
                # criticalValue is 'prob' or 'mean.prob.combo'
                splitRule = "logrank_surv"
              }
            } else{
              # Phase 2: Endpoint (CR: 'gray_cr' (default) and RE: 'gray_re')
              # splitRule 'csh_cr' must be directly specified
              splitRule = paste("gray", tolower(x = endPoint), sep = "_")
            }
              return( .VerifySplitRule(splitRule = splitRule, ...) ) # Recursively calls .VerifySplitRule to verify the value.
            })

setMethod(f = ".VerifySplitRule",
          signature = c(splitRule = "character"), # Defines method where splitRule is character input type.
          # Converts the splitRule to lowercase and checks if it is one of \{'logrank_surv', 'mean_surv', 'gray_cr', 'csh_cr', 'gray_re'\}.
          # If not, it generates an error message.
          definition = function(splitRule, ...) {

              splitRule <- tolower(x = splitRule)
              # print(splitRule)
              if (!(splitRule %in% c('logrank_surv', 'mean_surv', 'gray_cr', 'csh_cr', 'gray_re'))) {
                stop("splitRule must be one of {'logrank_surv', 'mean_surv', 'gray_cr', 'csh_cr', 'gray_re'}")
              }
              return( splitRule )
            })

# In summary, this code ensures that the splitRule parameter is either 'logrank_surv', 'mean_surv', 'gray_cr', 'csh_cr', 'gray_re'
# and converts it to lowercase if needed.
# It provides default values for splitRule based on the value of criticalValue if splitRule is NULL.
# These methods are used to standardize and validate the input for the splitRule parameter within the package's internal functions.
