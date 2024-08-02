# Class to store parameters that regulate tree and specify analysis preferences
#
# Class is not exported and is for internal convenience only
#
# Methods
#   .ParametersAsList(object, ...) {new; defined}
#
# Functions
# .parameters(endPoint, timePointsSurvival, timePointsEndpoint, timePoints, nTimes, response, response_endpoint, nTree, ERT, uniformSplit,
#                      randomSplit, splitRule, replace, nodeSize,
#                      minEvent, tieMethod, criticalValue,
#                      survivalTime, nSamples, pooled, stratifiedSplit)
#
#' @include class_TimeInfo.R criticalValue.R class_TreeType.R
#' @include class_TreeConditions.R
#'
setClass(Class = "SurvivalParameters_Mean", # before: Parameters_Mean
         contains = c("TimeInfo", "CriticalValueMean", "TreeType", "TreeConditions"))

setClass(Class = "SurvivalParameters_Area", # before: Parameters_Mean
         contains = c("TimeInfo", "CriticalValueArea", "TreeType", "TreeConditions"))

setClass(Class = "SurvivalParameters_Probability", # before: Parameters_Survival
         contains = c("TimeInfo", "CriticalValueSurv", "TreeType", "TreeConditions"))

setClass(Class = "EndPointParameters_Mean",
         contains = c("TimeInfo", "CriticalValueMean", "TreeType", "TreeConditions"))

setClass(Class = "EndPointParameters_Area", # before: Parameters_Mean
         contains = c("TimeInfo", "CriticalValueArea", "TreeType", "TreeConditions"))

setClass(Class = "EndPointParameters_Probability",
         contains = c("TimeInfo", "CriticalValueEndpoint", "TreeType", "TreeConditions"))

# Define a class union based on the selected classes
setClassUnion("AllParameters",
              members = c("SurvivalParameters_Mean",
                          "SurvivalParameters_Area",
                          "SurvivalParameters_Probability",
                          "EndPointParameters_Mean",
                          "EndPointParameters_Area",
                          "EndPointParameters_Probability"))

# Define a unique class name based on the combination of classes
ss_list_mean = c("SurvivalParameters_Probability", "SurvivalParameters_Mean")
ee_list_mean = c("EndPointParameters_Probability", "EndPointParameters_Mean")
ss_list_area = c("SurvivalParameters_Probability", "SurvivalParameters_Area")
ee_list_area = c("EndPointParameters_Probability", "EndPointParameters_Area")
ss_list = c(ss_list_mean, ss_list_area)
ee_list = c(ee_list_mean, ee_list_area)
class_name_mat = matrix(0, nrow = length(ss_list), ncol = length(ee_list))
for (ss in 1:length(ss_list)){
  ss_class = ss_list[ss]
  # print(ss_class)
  # Remove the common parts of the class names
  class_name1 <- gsub("Parameters_", "", ss_class)
  for (ee in 1:length(ee_list)){
    ee_class = ee_list[ee]
    # print(ee_class)
    class_name2 <- gsub("Parameters_", "", ee_class)
    # Concatenate the modified class names
    class_name <- paste0(class_name1, class_name2)
    # print(class_name)
    class_name_mat[ss,ee] = class_name
    if (!isS4(class_name)) {
      # print("class does not exist so creating new one.")
      A_name = sprintf("Param_%s", class_name)
      # print(A_name)
      setClass(A_name,
                   representation=representation(survivalparam=ss_class,
                                                 endpointparam=ee_class))
      # print("class created.")
    }
  }
}
# print(class_name_mat)

# setClassUnion(name = "SurvivalParameters",
#               members = c("SurvivalParameters_Mean", "SurvivalParameters_Probability"))

# setClassUnion(name = "EndPointParameters",
#               members = c("EndPointParameters_Mean", "EndPointParameters_Probability"))

# setClass("finalParams", contains = c("Param_SurvivalMeanEndPointMean",
#                                          "Param_SurvivalMeanEndPointProbability",
#                                          "Param_SurvivalProbabilityEndPointMean",
#                                          "Param_SurvivalProbabilityEndPointProbability"))

# SurvivalParameters_Mean
setMethod(f = "initialize",
          signature = c(.Object = "SurvivalParameters_Mean"),
          def = function(.Object, ...) {

            obj <- list(...)

            for (i in 1L:length(x = obj)) {
              if (is(object = obj[[ i ]], class2 = "TimeInfo")) {
                as(.Object, "TimeInfo") <- obj[[ i ]]
              } else if (is(object = obj[[ i ]], class2 = "CriticalValueBase")) {
                as(.Object, is(object = obj[[ i ]])[1L]) <- obj[[ i ]]
              } else if (is(object = obj[[ i ]], class2 = "TreeType")) {
                as(.Object, "TreeType") <- obj[[ i ]]
              } else if (is(object = obj[[ i ]], class2 = "TreeConditions")) {
                as(.Object, "TreeConditions") <- obj[[ i ]]
              } else {
                stop("SurvivalParameters_Mean: unrecognized object sent to Paremeters object")
              }
            }
            return( .Object )

          })

# SurvivalParameters_Area
setMethod(f = "initialize",
          signature = c(.Object = "SurvivalParameters_Area"),
          def = function(.Object, ...) {

            obj <- list(...)

            for (i in 1L:length(x = obj)) {
              if (is(object = obj[[ i ]], class2 = "TimeInfo")) {
                as(.Object, "TimeInfo") <- obj[[ i ]]
              } else if (is(object = obj[[ i ]], class2 = "CriticalValueBase")) {
                as(.Object, is(object = obj[[ i ]])[1L]) <- obj[[ i ]]
              } else if (is(object = obj[[ i ]], class2 = "TreeType")) {
                as(.Object, "TreeType") <- obj[[ i ]]
              } else if (is(object = obj[[ i ]], class2 = "TreeConditions")) {
                as(.Object, "TreeConditions") <- obj[[ i ]]
              } else {
                stop("SurvivalParameters_Area: unrecognized object sent to Paremeters object")
              }
            }
            return( .Object )

          })


#SurvivalParameters_Probability
setMethod(f = "initialize",
          signature = c(.Object = "SurvivalParameters_Probability"),
          def = function(.Object, ...) {

            obj <- list(...)
            # print(obj)

            for (i in 1L:length(x = obj)) {
              if (is(object = obj[[ i ]], class2 = "TimeInfo")) {
                as(.Object, "TimeInfo") <- obj[[ i ]]
              } else if (is(object = obj[[ i ]], class2 = "CriticalValueBase")) {
                as(.Object, is(object = obj[[ i ]])[1L]) <- obj[[ i ]]
              } else if (is(object = obj[[ i ]], class2 = "TreeType")) {
                as(.Object, "TreeType") <- obj[[ i ]]
              } else if (is(object = obj[[ i ]], class2 = "TreeConditions")) {
                as(.Object, "TreeConditions") <- obj[[ i ]]
              } else {
                stop("SurvivalParameters_Probability: unrecognized object sent to Paremeters object")
              }
            }

            return( .Object )

          })

# EndPointParameters_Mean
setMethod(f = "initialize",
          signature = c(.Object = "EndPointParameters_Mean"),
          def = function(.Object, ...) {

            obj <- list(...)

            for (i in 1L:length(x = obj)) {
              if (is(object = obj[[ i ]], class2 = "TimeInfo")) {
                as(.Object, "TimeInfo") <- obj[[ i ]]
              } else if (is(object = obj[[ i ]], class2 = "CriticalValueBase")) {
                as(.Object, is(object = obj[[ i ]])[1L]) <- obj[[ i ]]
              } else if (is(object = obj[[ i ]], class2 = "TreeType")) {
                as(.Object, "TreeType") <- obj[[ i ]]
              } else if (is(object = obj[[ i ]], class2 = "TreeConditions")) {
                as(.Object, "TreeConditions") <- obj[[ i ]]
              } else {
                stop("EndPointParameters_Mean: unrecognized object sent to Paremeters object")
              }
            }
            return( .Object )

          })

# EndPointParameters_Area
setMethod(f = "initialize",
          signature = c(.Object = "EndPointParameters_Area"),
          def = function(.Object, ...) {

            obj <- list(...)

            for (i in 1L:length(x = obj)) {
              if (is(object = obj[[ i ]], class2 = "TimeInfo")) {
                as(.Object, "TimeInfo") <- obj[[ i ]]
              } else if (is(object = obj[[ i ]], class2 = "CriticalValueBase")) {
                as(.Object, is(object = obj[[ i ]])[1L]) <- obj[[ i ]]
              } else if (is(object = obj[[ i ]], class2 = "TreeType")) {
                as(.Object, "TreeType") <- obj[[ i ]]
              } else if (is(object = obj[[ i ]], class2 = "TreeConditions")) {
                as(.Object, "TreeConditions") <- obj[[ i ]]
              } else {
                stop("EndPointParameters_Area: unrecognized object sent to Paremeters object")
              }
            }
            return( .Object )

          })

#EndPointParameters_Probability
setMethod(f = "initialize",
          signature = c(.Object = "EndPointParameters_Probability"),
          def = function(.Object, ...) {

            obj <- list(...)

            for (i in 1L:length(x = obj)) {
              if (is(object = obj[[ i ]], class2 = "TimeInfo")) {
                as(.Object, "TimeInfo") <- obj[[ i ]]
              } else if (is(object = obj[[ i ]], class2 = "CriticalValueBase")) {
                as(.Object, is(object = obj[[ i ]])[1L]) <- obj[[ i ]]
              } else if (is(object = obj[[ i ]], class2 = "TreeType")) {
                as(.Object, "TreeType") <- obj[[ i ]]
              } else if (is(object = obj[[ i ]], class2 = "TreeConditions")) {
                as(.Object, "TreeConditions") <- obj[[ i ]]
              } else {
                stop("EndPointParameters_Probability: unrecognized object sent to Paremeters object")
              }
            }

            return( .Object )

          })

#-------------------------------------------------------------------------------
# Function to verify inputs and create a Parameters object
#-------------------------------------------------------------------------------
# Function returns a Parameters object
#-------------------------------------------------------------------------------
.parameters <- function(endPoint,
                        timePointsSurvival,
                        timePointsEndpoint,
                        timePoints,
                        tau,
                        nTimes,
                        response,
                        response_endpoint,
                        nTree,
                        ERT,
                        uniformSplit,
                        randomSplit,
                        splitRule,
                        replace,
                        nodeSize,
                        minEvent,
                        tieMethod,
                        criticalValue1,
                        criticalValue2,
                        survivalTime,
                        endpointTime,
                        nSamples,
                        pooled,
                        stratifiedSplit) {

  # message("==================== Starting .Parameters Function ====================")

  # initialize TimeInfo
  # function returns a TimeInfo object
  timeInfo1 <- .timeInfo(timePointsPhase = timePointsSurvival, #timePointsPhase is priotized over timePoints/nTimes
                        timePoints = timePoints,
                        tau = tau,
                        nTimes = nTimes,
                        response = response
                        )
  timeInfo2 <- .timeInfo(timePointsPhase = timePointsEndpoint,
                         timePoints = timePoints,
                         tau = tau,
                         nTimes = nTimes,
                         response = response_endpoint
  )

  # timeInfo <- .timeInfo(timePointsSurvival = timePointsSurvival,
  #                        timePointsEndpoint = timePointsEndpoint,
  #                        timePoints = timePointsSurvival,
  #                        tau = tau,
  #                        nTimes = nTimes,
  #                        response = response
  # )
  # timeInfo1 = timeInfo@TimeInfoSurvival
  # timeInfo2 = timeInfo@TimeInfoEndpoint

  cv1 <- tolower(criticalValue1)
  cv2 <- tolower(criticalValue2)
  # message("Line215: cv1: ", cv1, "cv2: ", cv2)
  # message("Line216: criticalValue1: ", criticalValue1)
  # message("Line217: criticalValue2: ", criticalValue2)

  # initialize CriticalValue
  # function returns an object of class CriticalValueMean or
  #   CriticalValueSurv depending on input survivalTime
  # print("criticalValue1")
  # print(criticalValue1)

  # PHASE 1
  criticalValue1 <- .criticalValue(criticalValue = criticalValue1,
                                  Time = survivalTime,
                                  Step = "survival",
                                  tau = .Tau(object = timeInfo1),
                                  timePoints = .TimePoints(object = timeInfo1))
  # print("criticalValue1")
  # View(criticalValue1)

  # PHASE 2
  # NEED TO UPDATE ALL endpointTIME TO "ENDPOINTTIME FOR ALL SCRIPTS
  criticalValue2 <- .criticalValue(criticalValue = criticalValue2,
                                   Time = endpointTime, # change all scripts later to endpointTime
                                   Step = toString(endPoint),
                                   tau = .Tau(object = timeInfo2),
                                   timePoints = .TimePoints(object = timeInfo2))

# we dont need below because we just change endpointTime to a generic endpointTime
  # if (endPoint == "CR"){
  #   # print("CR")
  #   criticalValue2 <- .criticalValue(criticalValue = criticalValue2,
  #                                    Time = endpointTime,
  #                                    Step = toString(endPoint),
  #                                    tau = .Tau(object = timeInfo),
  #                                    timePoints = .TimePoints(object = timeInfo))
  #   # print("criticalValue2")
  #   # View(criticalValue2)
  # } else if (endPoint == "RE"){
  #   print("RE")
  #   criticalValue2 <- .criticalValue(criticalValue = criticalValue2,
  #                                    Time = RETime,
  #                                    Step = endPoint,
  #                                    tau = .Tau(object = timeInfo),
  #                                    timePoints = .TimePoints(object = timeInfo))
  # } else if (endPoint == "MC"){
  #   print("MC")
  #   criticalValue2 <- .criticalValue(criticalValue = criticalValue2,
  #                                    Time = MCTime,
  #                                    Step = endPoint,
  #                                    tau = .Tau(object = timeInfo),
  #                                    timePoints = .TimePoints(object = timeInfo))
  # } else{
  #   stop(".criticalValue within class_Parameters.R: EndPoint isn't right.")
  # }

  # print("end of .criticalValue")

  # initialize tree type info
  # function returns a TreeType object
  # print("class_Parameters.R: Line263")
  # print(sprintf("cv1: %s", cv1))
  treeType1 <- .treeType(endPoint = endPoint,
                        ERT = ERT,
                        nSamples = nSamples,
                        uniformSplit = uniformSplit,
                        replace = replace,
                        randomSplit = randomSplit,
                        splitRule = splitRule,
                        tieMethod = tieMethod,
                        criticalValue = cv1)
  # View(treeType1)
  # print(treeType1)
  # print("end of .treeType")

  splitRule2 = paste0(treeType1@splitRule,endPoint)
  message("class_Parameters.R: Line 354: Use the same rule type of splitrule for both Phase 1 and Phase 2 (", splitRule2, ")")
  # print(splitRule2)
  # print("Line279: cv2")
  # print(cv2)
  treeType2 <- .treeType(endPoint = endPoint,
                         ERT = ERT,
                         nSamples = nSamples,
                         uniformSplit = uniformSplit,
                         replace = replace,
                         randomSplit = randomSplit,
                         splitRule = splitRule2,
                         tieMethod = tieMethod,
                         criticalValue = cv2)

  # initialize tree conditions info
  # function returns a TreeConditions object
  treeConditions <- .treeConditions(nTree = nTree,
                                    nodeSize = nodeSize,
                                    minEvent = minEvent,
                                    pooled = pooled,
                                    stratifiedSplit = stratifiedSplit)
  # print("end of .treeConditions")

  # message("saving Parameters based on endpoints and critical values")
  # print("Line302:criticalValue1")
  # print(criticalValue1)
  # Step1 Survival
  if (is(object = criticalValue1, class2 = "CriticalValueSurv")) {
    print("Survival Parameters for Survival Probability")
    survparam = new(Class = "SurvivalParameters_Probability",
                # endPoint,
                timeInfo1,
                criticalValue1,
                treeType1,
                treeConditions)
  } else if (
    is(object = criticalValue1, class2 = "CriticalValueMean")) {
    print("Survival Parameters for Survival Mean")
    survparam = new(Class = "SurvivalParameters_Mean",
                # endPoint,
                timeInfo1,
                criticalValue1,
                treeType1,
                treeConditions)
  } else if (
    is(object = criticalValue1, class2 = "CriticalValueArea")) {
    print("Survival Parameters for Survival Area Under Curve")
    survparam = new(Class = "SurvivalParameters_Area",
                    # endPoint,
                    timeInfo1,
                    criticalValue1,
                    treeType1,
                    treeConditions)
  } else{
      stop("class_Parameters.R Line 393: we have object = criticalValue1 but class2 is neither CVSurv nor CVMean nor CVArea")
  }

  # Step2 EndPoint
  # print("Line328:criticalValue2")
  # print(criticalValue2)
  if (is(object = criticalValue2, class2 = "CriticalValueEndpoint")){
    print(sprintf("EndPoint %s Parameters for CIF Probability", endPoint))
    endpointparam = new(Class = "EndPointParameters_Probability",
                # endPoint,
                timeInfo2,
                criticalValue2,
                treeType2,
                treeConditions)
  } else if (is(object = criticalValue2, class2 = "CriticalValueMean")) {
    print("EndPoint Parameters for CIF Mean")
    endpointparam = new(Class = "EndPointParameters_Mean",
                # endPoint,
                timeInfo2,
                criticalValue2,
                treeType2,
                treeConditions)
  } else if (is(object = criticalValue2, class2 = "CriticalValueArea")) {
    print("EndPoint Parameters for CIF Area Under Curve")
    endpointparam = new(Class = "EndPointParameters_Area",
                        # endPoint,
                        timeInfo2,
                        criticalValue2,
                        treeType2,
                        treeConditions)
  } else{
    stop("class_Parameters.R Line 424: we have object = criticalValue2 BUT class2 is neither CVCR nor CVM")
  }

  # message("stage params done")

  # message("creating param")
  # User input (replace these with actual user input)
  survival_class <- class(survparam)[1]
  A1 <- gsub("Parameters_", "", survival_class)
  endpoint_class <- class(endpointparam)[1]
  A2 <- gsub("Parameters_", "", endpoint_class)

  # Create objects based on user input
  A = paste0("Param_",A1,A2);#print(A)
  params <- new(A,
                survivalparam = survparam,
                endpointparam = endpointparam)

  # params <- .createCombinedObject(survival_class, endpoint_class, survparam, endpointparam)
  # print(params)
  # print("returning params")
  return(params)
}


# # Define a function to create objects based on user input
# .createCombinedObject <- function(survival_class, endpoint_class, survparam, endpointparam) {
#
#   print("starting createCombinedObject")
#   print(survival_class)
#   print(endpoint_class)
#   print(survparam)
#   print(endpointparam)
#
#   print("creating object")
#   # Create an object of the combined class
#
#
#   return(params)
# }
