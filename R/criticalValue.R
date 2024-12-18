# Function to verify inputs and create appropriate CriticalValue object
#
# Function is not exported and for internal convenience only
#
# Function returns an object of class CriticalValueMean or CriticalValueSurv
# class_CriticalValueSurv.R is for surv prob #OLD: class_CriticalValueSurvival.R
#'
#' @include class_CriticalValue.R
#' @include class_CriticalValueMean.R
#' @include class_CriticalValueArea.R
#' @include class_CriticalValueSurv.R
#' @include class_CriticalValueEndpoint.R
#' @include VerifyCriticalValue.R
#' @include VerifySurvivalTime.R VerifyEndpointTime.R
#'
.criticalValue <- function(criticalValue,
                           Time,
                           Step,
                           tau,
                           timePoints) {
  # message("criticalValue.R")
  # message("-- Starting .criticalValue Function in criticalValue.R--")

  # ensure criticalValue is one of {'mean', 'prob', 'area', 'mean.prob.combo'}.
  # Methods return the original character possibly modified to be lower case.
  # print(criticalValue)
  criticalValue <- .VerifyCriticalValue(criticalValue = criticalValue)
  # message("end of .VerifyCriticalValue")
  # message("criticalValue: ", criticalValue)
  message("Step: ", Step)
  # message("tau is: ", tau)
  if (grepl("surv", Step, ignore.case = TRUE)){
    # print(".VerifySurvivalTime")
    # ensure that survivalTime is provided if criticalValue is one of
    # {'surv.prob', 'surv.man'}. Methods return the numeric survivalTime or NULL.
    survivalTime <- .VerifySurvivalTime(survivalTime = Time, #this has to be input Time
                                        criticalValue = criticalValue,
                                        tau = tau)
    # message("end of .VerifySurvivalTime")
    # message("survivalTime: ", survivalTime)
    if (!is.null(x = survivalTime)) {
      message("survivalTime: ", survivalTime)
      # if survivalTime is given as input, verify it and the timePoints input and
      # create a CriticalValueSurvivalProb object (used to be .criticalValueSurvival)
      return( .criticalValueSurv(survivalTime = survivalTime, #this has tobe survivalTime
                                     timePoints = timePoints,
                                     type = criticalValue) )
    } else{
      # message("Survival Mean - no evalTime needed")
      # if survivalTime nor endpointTime are not given as input, create a CriticalValueMean object
      if (criticalValue == "area"){
        # print("area")
        return( .criticalValueArea() )
      } else{
        # print("mean")
        return( .criticalValueMean() )
    }
    }
  } else if (Step == "CR" | Step == "RE") {
    # print(".VerifyEndpointTime")

    # ensure that endpointTime is provided if criticalValue is
    # {'cif'}. Methods return the numeric endpointTime or NULL.
    endpointTime <- .VerifyEndpointTime(endpointTime = Time, # this has to be Time input
                              criticalValue = criticalValue,
                              tau = tau)
    # message("endpointTime:", endpointTime)
    if (!is.null(x = endpointTime)){
      # message("endpointTime:",endpointTime)
      # if endpointTime is given as input, verify it and the timePoints input and
      # create a CriticalValueEndpoint object
      return( .criticalValueEndpoint(endpointTime = endpointTime, # this haas to be endpointTime, not Time
                               timePoints = timePoints,
                               type = criticalValue) )

    } else{
      # message("Mean Cumulative Incidence - no evalTime needed")
      # if survivalTime nor endpointTime are not given as input, create a CriticalValueMean object
      # same as above for if not survival prob
      if (criticalValue == "area"){
        return( .criticalValueArea() )
      } else{
        return( .criticalValueMean() )
      }
    }

  }
}
