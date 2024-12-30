# Class for storing primary results from survival analysis
#
# Class is not exported and is for internal convenience only
#
#  @slot txName A character object. The name of the treatment variable
#
#  @slot txLevels A vector. The treatment options available
#
#  @slot model A formula object. The model for extracting the covariates to be
#    considered for splitting
#
#  @slot survRF A SurvRF object. The primary results of the tree building
#    algorithm
#
#  @slot eligibility A logical vector object. TRUE indicates that the case was
#    eligible to be included in analysis
#
#  @slot valueAllTx A list object. The value of the tree for each tx level.
#
#  @slot optimal An Optimal object. The estimated optimal tx and optimal value.

# Methods
#
# Functions
#
#' @include class_Optimal.R class_SurvRF.R
#' @include class_Parameters.R

# Define new class "ITRSurvStep"
# Slots (Attributes):
# 1. txName (Slot Type: character): Stores the name of the treatment variable.
# 2. txLevels (Slot Type: vector): Holds a vector containing the treatment options available.
# 3. model (Slot Type: formula): Store a formula object representing the model for extracting the covariates to be considered for splitting in the analysis.
# 4. survRF (Slot Type: ANY): This slot is of type "ANY," which means it can store objects of any class. SurvRF object (primary results of the tree building algorithm).
# 5. eligibility (Slot Type: logical): Stores a logical vector that indicates whether each case was eligible to be included in the analysis.
# 6. valueAllTx (Slot Type: list): A list that stores values related to all treatment options.
# 7. optimal (Slot Type: Optimal): Stores an object of class "Optimal," from class_Optimal.R

setClass(Class = "ITRSurvStep",
         slots = c("txName" = "character",
                   "txLevels" = "vector",
                   "model" = "formula",
                   "survRF" = "ANY",
                   "eligibility" = "logical",
                   "valueAllTx" = "list",
                   "optimal" = "Optimal"))

#-------------------------------------------------------------------------------
# function to return key stage results as a list
#-------------------------------------------------------------------------------
# function is not exported
#-------------------------------------------------------------------------------
# .stage <- function(object, ...) {
#   result <- list()
#
#   result[[ "txName" ]] <- object@txName
#   result[[ "txLevels" ]] <- object@txLevels
#   result[[ "sampleSize" ]] <- sum(object@eligibility)
#   result[[ "valueAllTx" ]] <- object@valueAllTx
#   result[[ "optimal" ]] <- .OptimalAsList(object = object)
#   result[[ "stages" ]] <- .stageSurvRF(object = object@survRF)
#
#   return( result )
# }


.meanValue <- function(object,Phase,endPoint, ...) {
  res <- list()
  # message("OK1")

  # print(endPoint)
  object_surv <<- object[[1]]
  object_ep <<- object[[2]]
  # message("OK2")

  if (grepl("surv", Phase, ignore.case = TRUE)){
    # object = object_surv
    optfunc = max

    AU_name = "AUS"
    mean_name = "Et_survival"
    prob_name = "S(t)"

    # message("OK3")
    area_dat = cbind(object_surv@valueAllTx$AUS[[1]],
                     object_surv@valueAllTx$AUS[[2]])

    # message("OK3")
    mean_dat = cbind(object_surv@valueAllTx$mean[[1]],
                     object_surv@valueAllTx$mean[[2]])

    mdat1 <<- mean_dat
    # message("OK3")

    if (!is.null(object_surv@valueAllTx$Prob[[1]])) {
      # message("Calculate mean Survival probability curve")
      # View(object_surv@valueAllTx$Prob)
      prob_dat0 <<- object_surv@valueAllTx$Prob[[1]]
      prob_dat1 <<- object_surv@valueAllTx$Prob[[2]]
      prob_dat <<- cbind(prob_dat0, prob_dat1)
      res[[ prob_name ]] <-  mean(x = apply(X = prob_dat,
                                            MARGIN = 1L,
                                            FUN = optfunc))
    }
  }

  if (!is.null(endPoint)){
    # print(endPoint)
    if (endPoint == Phase){
      object = object_ep
      objep <<- object_ep
      objs <<- object_surv
      optfunc = min

      AU_name = "AUC"

      if (endPoint == "CR"){
        mean_name = "Et_cif" #"Et_cif"
        prob_name = "F_j(t)"
      } else if (endPoint == "RE"){
        mean_name = "Et_mff"
        prob_name = "mu(t)"
      }

      # message("OK4")

      indicator_vector <- object_surv@optimal@Ratio_Stopping_Ind
      dat_mean <- cbind(object_ep@valueAllTx$mean[[1]],
                        object_ep@valueAllTx$mean[[2]])
      # message("Obtain mean CIF times for those who continued to Phase 2 with mean critical values.")
      mean_dat <- dat_mean[indicator_vector == 0,]
      if (is.null(nrow(mean_dat))){
        mean_dat = rbind(mean_dat)
      }
      mdat2 <<- mean_dat
      # print("hi")
      if (endPoint == "CR"){
        # print('hi2')
        dat_area <- cbind(object_ep@valueAllTx$AUS[[1]],
                          object_ep@valueAllTx$AUS[[2]])
        # message("Obtain mean area under CIF curve for those who continued to Phase 2 with area critical values.")
        area_dat <- dat_area[indicator_vector == 0,]
        if (is.null(nrow(area_dat))){
          area_dat = rbind(area_dat)
        }
        # print("hi3")
      }
      if (!is.null(object_ep@valueAllTx$Prob[[1]])) {
        # message("Obtain mean CIF curve for those who continued to Phase 2 with prob critical values.")
        dat_prob0 <- object_ep@valueAllTx$Prob[[1]]
        dat_prob1 <- object_ep@valueAllTx$Prob[[2]]
        prob_dat <- cbind(dat_prob0, dat_prob1)[indicator_vector == 0,]
        # print(typeof(prob_dat))
        # print(prob_dat)
        if (nrow(prob_dat) == 1){
          res[[prob_name]] = optfunc(prob_dat)
        } else{
          res[[ prob_name ]] <-  mean(x = apply(X = prob_dat,
                                                MARGIN = 1L,
                                                FUN = optfunc))}
      }
    }
  }
  # print(endPoint)
  # message("OK6")
  # print(optfunc)
  mdat3 <<- mean_dat
  of <<- optfunc
  # returns the mean of the expected times
  # first take the minimum across rows (trts), then take the mean of all the mins (across patients)
  res[[ mean_name ]] <-  mean(x = apply(X = mean_dat,
                                        MARGIN = 1L,
                                        FUN = optfunc))
  if (endPoint == "CR"){
    # print('hi6')
    res[[ AU_name ]] <-  mean(x = apply(X = area_dat,
                                        MARGIN = 1L,
                                        FUN = optfunc))
  } else{ # N/A for RE
    # print('hi7')
    res[[ AU_name ]] <- NULL
  }

  return( res )
}

# generic defined in class_Optimal.R
setMethod(f = ".OptimalAsList",
          signature = c(object = "ITRSurvStep"),
          definition = function(object, ...) {
                return( .OptimalAsList(object = object@optimal) )
              })

#-------------------------------------------------------------------------------
# method to retrieve predicted values
#-------------------------------------------------------------------------------
# if findOptimal is TRUE, method stops with error
# if findOptimal is FALSE, method returns a Value object
#-------------------------------------------------------------------------------
setMethod(f = ".Predict",
          signature = c(object = "ITRSurvStep",
                        newdata = NULL),
          definition = function(object,
                                newdata,
                                Phase,
                                ...) {
            print(".Predict from class_ITRSurvStep.R when newdata = NULL: LINE 165")

              return( .Predict(object = object@survRF,
                               newdata = NULL,
                               Phase,
                               ...) )

            })

#-------------------------------------------------------------------------------
# method to predict value for new data
#-------------------------------------------------------------------------------
# if findOptimal is TRUE, method returns a list containing a Value object and
#   an Optimal object
# if findOptimal is FALSE, method returns a Value object
#-------------------------------------------------------------------------------
#' @include class_Optimal.R
#' @importFrom stats model.frame
#'
setMethod(f = ".Predict",
          signature = c(object = "ITRSurvStep",
                        newdata = "data.frame"),
          definition = function(object,
                                newdata,
                                Phase,
                                epName1,
                                endPoint,
                                ...,
                                eps0,
                                params,
                                findOptimal) {
            print("class_ITRSurvStep.R: LINE 195")
            # print(".Predict from class_ITRSurvStep.R when when newdata = dataframe")
            # print(sprintf("Phase: %s", Phase))

              # update model to remove response
              mod <- update(object@model, NULL ~ .)

              # ensure data contains all model covariates
              x <- tryCatch(expr = stats::model.frame(formula = mod,
                                                      data = newdata),
                            error = function(e) {
                                      stop("variables in the training data missing in newdata",
                                           call. = FALSE)
                                     })

              # remove response from x
              # this should no longer happen, but keeping anyway
              if (attr(x = terms(x = mod), which = "response") == 1L) {
                x <- x[,-1L,drop = FALSE]
              }

              if (findOptimal) {
                print("finding optimal. returning .PredictAll()")
                # if optimal is requested, make predictions for all possible
                # treatment options
                # print(3)
                # print(eps0)
                if (is.language(eps0)){
                  eps0 = c(eps0[[2]],eps0[[3]])
                }
                # View(eps0)
                print("resV: class_ITRSurvStep.R: line 268")
                message("epName1:", epName1)
                resV <- .PredictAll(object = object@survRF,
                                    Phase = Phase,
                                    epName1 = epName1,
                                    endPoint = endPoint,
                                    eps0 = eps0,
                                    newdata = newdata,
                                    params = params,
                                    model = mod,
                                    txName = object@txName,
                                    txLevels = object@txLevels)
                # print(4)
                if (Phase == 1){
                  resV_surv <<- resV
                  # message("resV_surv saved")
                } else if (Phase == 2){
                  # message("resV_ep saved")
                  resV_ep <<- resV
                }
                message("class_ITRSurvStep.R: LINE 237: resV_surv and resV_ep saved.")
                return(list("value" = resV$predicted, "optimal" = resV$optimal) )

              } else {
                print('not finding optimal. returning .Predict() instead of .PredictAll()')
                return( .Predict(object = object@survRF,
                                 newdata = x,
                                 params = params, ...) )
              }
          })


# Internal function to create random forest
#
# Function is not exported
#
# @param model A survival formula object for Phase 1, the rhs of which specifies the
#   covariates to be considered in the splitting algorithm

# @param model A survival formula object for Phase 2, the rhs of which specifies the
#   covariates to be considered in the splitting algorithm
#
# @params data A list for Phase 1 and Phase 2 data.frame objects containing covariate and treatment histories
#   For CR endpoint: they are the same.
#   For RE endpoint: they are different (one is for failure and one is for recurrent events)
#
# @params params A Parameters object.
#
# @params txName A character object or a character vector object. The names
#   of the treatment variables in data
#
# @params mTry An integer object or NULL. If integer, the maximum number of
#    covariates to consider for each split. If NULL, mTry is sqrt(np)
#
# @params sampleSize An integer object. The number of samples to draw for each
#    tree
#
#' @importFrom stats na.pass
#' @importFrom stats update
#' @importFrom stats terms
#' @importFrom stats complete.cases
#' @importFrom stats model.frame
#' @importFrom stats model.response
#' @include shiftMat.R survRF.R
.itrSurvStep <- function(...,
                         Phase,
                         eps0 = NULL,
                         epName,
                         idName,
                         model_surv,
                         model_ep,
                         # model_cause,
                         endPoint,
                         data,
                         ord_causeind,
                         ord_response,
                         params,
                         txName,
                         mTry,
                         sampleSize,
                         person_indicator) {
  # print("---STARTING ITRSURVSTEP---")

  if (Phase == "RE"){
    phase1_params <<- params@survivalparam
    params = params@endpointparam
  }
  ###########################################################################
  ############################ PHASE 1: SURVIVAL ############################
  ###########################################################################

  if (Phase != toString(endPoint)){

    dataset = data[[1]]
    # this is Phase 1: Survival

    mod <- model_surv$models
    # identify order 1 terms in formula
    order1 <- attr(x = stats::terms(x = mod), which = "order") == 1L
    if (any(order1)) {
      stageCov <- attr(x = stats::terms(x = mod), which = "term.labels")[order1]
    } else {
      stop("problem in identifying covariates, verify formula\n", call. = FALSE)
    }

    # warn about order > 1
    orderHigh <- attr(x = stats::terms(x = mod), which = "order") > 1L
    if (any(orderHigh)) message("interaction terms are ignored")

    # extract model frame
    x <- stats::model.frame(formula = mod,
                            data = dataset, # dataset for survival
                            na.action = na.pass)
    x_ep = x
    message("model ", appendLF = FALSE)
    tm <- as.character(mod)
    message(tm[2], " ~ ", tm[3])

    # identify individuals with complete data
    elig <- stats::complete.cases(x)

    # extract response and delta from model frame
    response0 <- stats::model.response(data = x)
    delta <- response0[,2L] # survival delta
    response_surv <- response0[,1L]

    if (endPoint == "CR"){
      # print('setting delta_cause = delta but we wont use it for survival part Phase1')
      delta_endpoint = delta # CR delta which is delta in Phase 1
    }
    if (endPoint == "RE"){
      # make sure this works with existing code...
      delta_endpoint = delta # dont need this variable for Phase 1 for RE
    }
    response = response_surv

    ###########################################################################
    ############################ PHASE 2: ENDPOINT ############################
    ###########################################################################
  } else{ # endpoint Phase 2
    dataset <<- data[[2]]
    testing_new_data <<- data[[1]]
    # print(0)

    mod <<- model_ep$models # for endpoint model
    mod_surv <<- model_surv$models
    # print(1)

    # identify order 1 terms in formula
    order1 <- attr(x = stats::terms(x = mod), which = "order") == 1L
    if (any(order1)) {
      stageCov <- attr(x = stats::terms(x = mod), which = "term.labels")[order1]
    } else {
      stop("problem in identifying covariates, verify formula\n", call. = FALSE)
    }
    # print(3)

    # warn about order > 1
    orderHigh <- attr(x = stats::terms(x = mod), which = "order") > 1L
    if (any(orderHigh)) message("interaction terms are ignored")

    # extract model frame
    x_endpoint <<- stats::model.frame(formula = mod,
                                     data = dataset,
                                     na.action = na.pass)
    message("model_endpoint ", appendLF = FALSE)
    tm <- as.character(mod)
    message(tm[2], " ~ ", tm[3])

    # identify individuals with complete data
    elig <<- stats::complete.cases(x_endpoint)
    # print(2)
    if (Phase == "RE"){
      long_x <<- cbind(dataset %>%
              dplyr::select(!!sym(idName)), x_endpoint)
      # elig <<- long_x %>%
      #   mutate(eligible = complete.cases(.)) %>%  # Determine eligibility for each row
      #   group_by(!!sym(idName)) %>%                          # Group by individual (ID)
      #   summarise(elig1 = all(eligible)) %>%       # Check if all rows are eligible for the person
      #   pull(elig1)                             # Extract the eligibility vector for each ID
      elig0 <<- long_x %>%
        mutate(eligible = complete.cases(.)) %>%    # Determine eligibility for each row
        group_by(!!sym(idName)) %>%                  # Group by individual (ID)
        summarise(elig00 = all(eligible)) %>%         # Check if all rows are eligible for the person
        dplyr::select(!!sym(idName), elig00)
      elig <<- merge(long_x, elig0, by = idName) %>%
        pull(elig00)
    }

    # extract response and delta from model frame
    response_tmp <<- stats::model.response(data = x_endpoint)
    if (endPoint == "CR"){
      delta_endpoint <<- response_tmp[,2L] # priority cause delta
      response_endpoint <<- response_tmp[,1L]
      response = response_endpoint # response of endpoint dataset
      x_surv <<- stats::model.frame(formula = mod_surv,
                                    data = dataset,
                                    na.action = na.pass)
      response_surv <<- stats::model.response(data = x_surv)
      delta <<- response_surv[,2L] # survival delta
      x = x_endpoint # covariates of endpoint dataset (this is same for CR; diff for RE)
      x_ep = x
    }
    if (endPoint == "RE"){
      delta_endpoint <<- response_tmp[,3L] # RE delta
      response_endpoint_start <<- response_tmp[,1L]
      response_endpoint_stop <<- response_tmp[,2L]
      response = response_endpoint_stop
      # Create d3 by joining data[[2]] with data[[1]] based on "id"
      data1 <<- data[[1]]
      data2 <<- data[[2]]
      data3 <<- data[[3]]
      var_names <- all.vars(mod_surv[[3]])
      # Extract the parts from `models_RE`
      survival_status_name <<- all.vars(model_surv[[1]])[2]
      tstop_name <<- all.vars(model_ep[[1]])[2]

      d3 <<- data2 %>%
        dplyr::select(!!sym(idName), !!sym(survival_status_name),
                      all_of(var_names)) %>%
        left_join(data1 %>%
                    dplyr::select(!!sym(idName),
                                  !!sym(tstop_name)),
                  by = idName)

      x_surv <<- stats::model.frame(formula = mod_surv,
                                    # we want to Phase 2 Dataset because multiple records per person for RE; and for CR they are the same dataset
                                    data = d3,
                                    na.action = na.pass)
      response_surv <<- stats::model.response(data = x_surv)
      delta <<- response_surv[,2L] # survival delta
      x = stats::model.frame(formula = mod_surv,
                             # covariates of survival dataset b/c 1 row per person
                             data = data1,
                             na.action = na.pass)
      x_ep = x_endpoint
    }
  }

  # View(data)
  # print(3)
  # print(data[[3]])
  # set id_vec which is needed for Phase2RE, in fortran, to calculate mff stuff with pr2 to get at risk for death in RE setting
  id_vec <<- data[[3]] %>% unlist()
  row_vec <<- c(1:length(id_vec))
  # idvec_test <<- id_vec
  # print(id_vec)
  # stop(" testing id_vec ")

 # old stuff located in scratch: old itrsurvstep code.R

  # remove response from x TO JSUT GET COVARIATES
  # print(head(x,10))
  if (attr(x = terms(x = mod), which = "response") == 1L) {
    x <- x[,-1L,drop = FALSE]
  }
  if (attr(x = terms(x = mod), which = "response") == 1L) {
    x_ep <- x_ep[,-1L,drop = FALSE]
  }
  # print(head(x,10))
  # if (Phase == "RE"){stop('tesitng x')}

  # we aren't doing multistage.
  # # responses that are zero indicate censored at a previous stage
  # zeroed <- abs(x = response) < 1e-8
  # elig <- elig & !zeroed

  if (sum(elig) == 0L) stop("no cases have complete data", call. = FALSE)
  # message("cases in stage: ", sum(elig))

  # maximum number of covariates to try
  if (is.null(x = mTry)) {
    mTry <- as.integer(x = ceiling(x = sqrt(x = ncol(x = x))))
    # message("maximum # of covariates considered for splitting set to ", mTry)
  } else if (mTry > ncol(x = x)) {
    # message("mTry reset as it is larger than the # of available covariates")
    mTry <- as.integer(x = ceiling(x = sqrt(x = ncol(x = x))))
    # message("maximum # of covariates considered for splitting set to ", mTry)
  } else {
    mTry <- as.integer(x = mTry)
  }
  # message('mTry: ', mTry)

  ### Transform the time variable to a probability mass vector ###
  # identify time points <= response
  # resp.tmp <<- response[elig]
  tp.tmp <<- .TimePoints(object = params)
  if (Phase == "RE"){
    tp.tmp.phase1<<- .TimePoints(object = phase1_params)
  }

  # print(4)
  # we use tSurv to get pr
  # tSurv shows at risk to not at risk.
  tSurv <- sapply(X = response[elig],
                  FUN = function(s, tp) { as.integer(x = {s < tp}) }, # 0 means at risk; 1 means NOT at risk
                  tp = .TimePoints(object = params)) # tp = sort(unique(timePointsEndpoint1)))
  tSurv.tmp<<-tSurv
  # time point nearest the status change (response (death), censoring, recurrent event, CR,etc) without going over
  # for RE: pr based on TStop times; pr2 based on the interval [TStart, TStop)

  # pr is the first time a person changes from being at risk (0) to not being at risk (1)
  # for pr: 1 indicates there was some sort of status change BEFORE the next tp
  # doesn't matter if it's censoring or death or a recurrent event
  # pr gives the indicator for the time of at-risk status change (either failure or censoring)
  # the time point at which they become not at risk anymore
  # in FORTRAN, we further modify to calculate at-risk or # events by multiplying by indicators of events of interest
  # and use pr2 for RE to get subjects at-risk
  # We also transpose before sending to Fortran

  # pr: {nTimes x nElig}
  pr <- {rbind(tSurv[-1L,],1)-tSurv}
  pr.tmp <- pr
  # colnames(pr.tmp) = tp.tmp
  # rownames(pr.tmp) = response_surv[,1]
  pr.tmp <<- pr.tmp
  # print(5)

  if (any(is.na(x = pr))) stop("NA not permitted in pr -- contact maintainer",
                               call. = FALSE)
  if (any(pr > 1.0) || any(pr < 0.0)) {
    stop("pr must obey 0 <= pr <= 1 -- contact maintainer", call. = FALSE)
  }

  if (Phase == "RE" & endPoint == "RE"){
    # print("pr_surv")
    # we use tSurv_surv to get pr_surv
    # changed response to be response_surv
    response_survival <<- response_surv[,1L]
    tSurv_surv <- sapply(X = response_survival[elig],
                         FUN = function(s, tp) { as.integer(x = {s < tp}) }, # 0 means at risk for death; 1 means NOT at risk
                         tp = .TimePoints(object = phase1_params)) # this is phase 1 TIME POINTS!!!
    tSurv_surv.tmp<<-tSurv_surv
    ttp <<- .TimePoints(object = phase1_params)
    # rownames(tSurv_surv.tmp) = ttp
    # colnames(tSurv_surv.tmp) = a1$id
    # time point nearest the status change (response (death) or censoring ONLY) without going over
    # prsurv: {nTimes_SURVIVAL x nElig (still records not people)}
    pr_surv <- {rbind(tSurv_surv[-1L,],1)-tSurv_surv}
    pr_surv.tmp <<- pr_surv
    # print(ttp)
    # print(rownames(pr_surv.tmp))
    # rownames(pr_surv.tmp) = ttp
    # colnames(pr_surv.tmp) = a1$id
    # print(rownames(pr_surv.tmp))
    # print(pr_surv.tmp)
    # View(pr_surv.tmp)
    # stop("stopping")

    if (any(is.na(x = pr_surv))) stop("NA not permitted in pr_surv -- contact maintainer",
                                      call. = FALSE)
    if (any(pr_surv > 1.0) || any(pr_surv < 0.0)) {
      stop("pr_surv must obey 0 <= pr_surv <= 1 -- contact maintainer", call. = FALSE)
    }

    # FOR RECURRENT EVENTS ONLY:
    # obtain pr2 to obtain number of people at risk for recurrent event set-up
    # pr2 shows 'at risk' to 'not at risk' for records in recurrent event set-up
    response_re <<- cbind(response_endpoint_start, response_endpoint_stop)
    elig<<-elig
    # print(elig)
    # print(response_re[elig,])
    # Apply function to each row of response_re
    ttpp <<- .TimePoints(object = params)
    pr2 <- apply(response_re[elig,], 1, function(row) {
      # start-stop interval is OPEN-CLOSED: (start, stop] when it's NOT a "first" record
      # start-stop interval is CLOSED-CLOSED: [start, stop] when it IS a "first" record
        # because we include 0, we need make it so that left hand interval is closed when it starts with 0
        # but when it's not starting with 0 because its a 2nd+ record, then we want it to be OPEN.
      # THIS IS TO IDENTIFY # AT RISK FOR RECORDS
      # row[1] is the start, row[2] is the stop

      sapply(.TimePoints(object = params), function(tp) {
        # 1 means at risk; 0 means NOT at risk (opposite of tSurv)
        #this is old, didn't account for 0: as.integer(row[1] < tp & row[2] >= tp)
        as.integer(((row[1] < tp & row[1] > 0) | (row[1] <= tp & row[1] == 0)) & row[2] >= tp)
      })
    })
    print(dim(pr2))
    pr2.tmp<<-pr2
    # August 2024: number of ppl at risk for recurrent event and death are the same because
    # if they are at-risk for RE then they are also at risk for death and vise versa
    # so, we can use pr2 for both in fortran when calculating at risk, which applies to both
    # (i.e., for KM survival and dRhat for calculating mff mu) # Ghosh and Lin, 2000
    # {nTimes x nElig} # should be same dimensions as pr
    if (any(is.na(x = pr2))) stop("NA not permitted in pr2 -- contact maintainer",
                                  call. = FALSE)
    if (any(pr2 > 1.0) || any(pr2 < 0.0)) {
      stop("pr2 must obey 0 <= pr2 <= 1 -- contact maintainer", call. = FALSE)
    }

    # pr2_surv: September 2024: not most efficient but easier to just plug into fortran b/c tp_surv is diff from tp_end and surv should reflect correct tp

    # print("pr2_surv")
    # FOR RECURRENT EVENTS ONLY:
    # obtain pr2_surv to obtain number of people at risk for death during recurrent event Phase
    # pr2_surv shows 'at risk' to 'not at risk' for people in recurrent event set-up
    # Apply function to each row of response_re
    pr2_surv <- apply(response_re[elig,], 1, function(row) {
      # start-stop interval is OPEN-CLOSED: (start, stop]
      # THIS IS TO IDENTIFY # AT RISK FOR RECORDS
      # row[1] is the start, row[2] is the stop
      # USING SURVIVAL TIME POINTS (includes 0 and tau)
      sapply(.TimePoints(object = phase1_params), function(tp) {
        # 1 means at risk; 0 means NOT at risk (opposite of tSurv)
        # old: doens't account for 0 : as.integer(row[1] < tp & row[2] >= tp)
        as.integer(((row[1] < tp & row[1] > 0) | (row[1] <= tp & row[1] == 0)) & row[2] >= tp)
      })
    })
    pr2_surv.tmp<<-as.matrix(pr2_surv)

    # {nTimes_survival x nElig}
    if (any(is.na(x = pr2_surv))) stop("NA not permitted in pr2_surv -- contact maintainer",
                                       call. = FALSE)
    if (any(pr2_surv > 1.0) || any(pr2_surv < 0.0)) {
      stop("pr2_surv must obey 0 <= pr2_surv <= 1 -- contact maintainer", call. = FALSE)
    }

    # stop("testing")

  } else{
    message("we set pr2=pr2_surv=pr_surv=pr since we don't use pr2,pr2_surv,pr_surv in this setting")
    pr2 = pr
    pr2_surv = pr
    pr_surv = pr
  }

  # identify tx levels in limited data
  if (is.factor(x = dataset[,txName])) {
    txLevels <- levels(x = factor(x = dataset[elig,txName]))
  } else {
    txLevels <- sort(x = unique(x = dataset[elig, txName]))
  }

  if (length(x = txLevels) == 1L) {
    message("***only one treatment level in dataset***")
  }


  # print('delta[elig]')
  # print(length(delta[elig]))
  # print(delta[elig])
  # print("ord_causeind[elig]")
  # print(length(ord_causeind[elig]))
  # print(ord_causeind[elig])
  #
  # print('delta')
  # print(length(delta))
  # print(delta)
  # print("ord_causeind")
  # print(length(ord_causeind))
  # print(ord_causeind)


  if (grepl("surv", Phase, ignore.case = TRUE) | Phase == 1){
    # print(Phase)
    # print(id_vec)
    id_vec = id_vec %>% unique()
    row_vec = row_vec
    # print(id_vec)
    # stop()
  }

  if (.Pooled(object = params)) {
      message("pooled analysis; treatments ", paste(txLevels,collapse=" "))
      # this will be a SurvRF object
      result <- .survRF(endPoint = endPoint,
                        Phase = Phase,
                        eps0 = eps0,
                        x = x[elig,,drop=FALSE],
                        x_ep = x_ep[elig,,drop=FALSE],
                        y = response[elig],
                        idvec = id_vec[elig],
                        rowvec = row_vec[elig],
                        pr = pr,
                        pr2 = pr2,
                        pr2_surv = pr2_surv,
                        pr_surv = pr_surv,
                        ord_causeind = ord_causeind[elig],
                        ord_response = ord_response[elig],
                        delta = delta[elig],
                        delta_endpoint = delta_endpoint[elig],
                        person_indicator = person_indicator[elig],
                        params = params,
                        mTry = mTry,
                        txLevels = txLevels,
                        model = mod,
                        sampleSize = sampleSize)

  } else {
    # message("stratified analysis")
    # result will be a list of SurvRF objects
    result <- list()
    # message("number of txLevels:", length(txLevels))
    # print(txLevels)
    # print(pr_surv)
    for (i in 1L:length(x = txLevels)) {
      # message("-----------")
      message("  treatment level ", txLevels[i])
      nms <<- as.character(x = txLevels[i])
      # Creates a logical vector di that checks if the treatment level in the dataset dataframe matches the current treatment level.
      # This is used for subsetting the dataset.
      di <<- {dataset[elig,txName] == txLevels[i]}
      # Creates another logical vector use that combines eligibility (elig) with the condition that the treatment variable (txName)
      # matches the current treatment level.

      # This is used to subset the dataset for the current treatment level.
      use <<- elig & {dataset[,txName] == txLevels[i]}
      # print(sprintf("starting .SurvRF for treatment level: %s", txLevels[i]))

      result[[ nms ]] <- .survRF(endPoint = endPoint,
                                 Phase = Phase,
                                 eps0 = eps0,
                                 x = x[use,,drop=FALSE], # subset of covariates for the current treatment level
                                 x_ep = x_ep[use,,drop=FALSE],
                                 y = response[use],
                                 idvec = id_vec[use], # only matters for Phase2RE
                                 rowvec = row_vec[use],
                                 pr = pr[,di],
                                 pr2 = pr2[,di],
                                 pr2_surv = pr2_surv[,di],
                                 pr_surv = pr_surv[,di],
                                 ord_causeind = ord_causeind[use], # only matters for Phase2CR
                                 ord_response = ord_response[use], # only matters for Phase2CR
                                 delta = delta[use],
                                 delta_endpoint = delta_endpoint[use],
                                 person_indicator = person_indicator[use],
                                 params = params,
                                 mTry = mTry,
                                 txLevels = txLevels[i],
                                 model = mod,
                                 sampleSize = sampleSize)
      # print(Phase)
      if (Phase == "Survival"){
        if (nms == 1){
          res_1_1 <<- result[[nms]]
          # View(res_1_1)
        } else{
          res_1_2 <<- result[[nms]]
          # View(res_1_2)
        }
      } else{
        if (nms == 1){
          res_2_1 <<- result[[nms]]
          # View(res_2_1)
        } else{
          res_2_2 <<- result[[nms]]
          # View(res_2_2)
        }
      }
    }
    # print(0)
    result <- new(Class = "SurvRFStratified", "strat" = result)
  }

  # print("====================================== end of forest =================================")

  # calculate the estimated values for all treatment levels
  # .PredictAll() is a method; called here for objects of class SurvRF, which
  # is defined in file class_SurvRF.R
  # print("WHAT2")
  # print("epName")
  # print(epName)

  # stop("testing")
  resV <- .PredictAll(Phase = Phase,
                      eps0 = eps0,
                      epName1 = epName,
                      endPoint = endPoint,
                      object = result,
                      newdata = dataset[elig,],
                      params = params,
                      model = mod,
                      txName = txName,
                      txLevels = txLevels)
  # if (Phase == "Survival"){stop("testing mean")}
  # print("WHAT3")
  # View(dataset[elig,])
  # message("View(resV)")
    # View(resV)
  if (Phase == "RE"){
    # View(resV[["predicted"]][["Func"]][[1]])
    # View(resV[["predicted"]][["Func"]][[2]])
    # View(resV[["predicted"]][["mean"]])
    # stop("testing resV")
  }

  result <- new(Class = "ITRSurvStep",
                "txName" = txName,
                "txLevels" = txLevels,
                "model" = mod,
                "survRF" = result,
                "eligibility" = elig,
                "valueAllTx" = resV$predicted,
                "optimal" = resV$optimal)
  # print("WHAT4")
  return( result )

}
