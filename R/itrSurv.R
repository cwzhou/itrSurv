#' Individualized Treatment Regime for Survival Analysis
#'
#' Provides methods for estimating single stage optimal individualized treatment
#'   regimes for survival outcomes with multiple endpoints. Currently only applicable to single stage diseases.
#'   1) Competing risk with priority cause.
#'   2) Recurrent events.
#'
#' @param ... Ignored. Present only to require named inputs.
#'
#' @param data A data.frame object.
#'   For CR endpoint: The full dataset including treatments received, all covariates,
#'   observed times, overall failure indicator, and cause 1 failure indicator.
#'   MUST BE SORTED BY OBSERVED FAILURE TIME IN ASCENDING ORDER for CR endpoint!
#'   For RE endpoint: The first column must be Individual ID. The full dataset including
#'   ID, treatment, covariates, start times, stop times,
#'   recurrent event indicator, death indicator.
#'   Can be provided as a matrix object if column headers are included.
#'   Can contain missing data coded as NA, but cannot contain NaN.
#'
#' @param endPoint A character object. Must be one of
#'   \{"CR", "RE"\}. The label for the ultimate end point type.
#'   For "CR": competing risks data to examine cumulative incidence function;
#'   For "RE": recurrent events data to examine mean frequency function
#'
#' @param yName A character object. The observed failure time variable name. Only required for endpoint CR.
#'
#' @param idName A character object. The ID variable name. Only required for endpoint RE.
#'
#' @param txName A character object. The treatment variable name.
#'
#' @param epName A character object. The endpoint variable name of interest.
#'   For "RE" endpoint: Refers to indicator if row is recurrent event (1) or not (0).
#'   Later used to subset Phase 1 survival dataset (one row per subject).
#'   For "CR" endpoint: Refers to status indicator for causes: 0 = censored, 1 = priority cause, 2 = any other cause.
#'   Should be inputted as one row per person for CR data.
#'   If doing csh_cr test for "CR" endpoint, then epName is not needed.
#'
#' @param models A list containing formulas defining the response as a Surv() object
#'   and the covariate structure of the model.
#'   Note that this model should not include any terms of order > 1.
#'    For CR endpoint: this should be a list where the first is formula for
#'    overall survival and the second is a formula for priority cause.
#'    For RE endpoint: the first list should be a formula for terminal (i.e death) events
#'    and the second should be a formula for recurrent events.
#'    Example: CR:
#'       Surv(obs_time, D.0) ~ Z1
#'       Surv(obs_time, D.1) ~ Z1
#'    Example: RE:
#'       Surv(TStop, D.0) ~ Z1
#'       Surv(TStart, TStop, D.1) ~ Z1
#'    For RE, TStart is start time, TStop is stop time, and D.0 represents death indicator for survival dataset
#'    D.1 represents recurrent event indicator for RE dataset (full dataset)
#'    The inputted dataset must reflect these variable names.
#'
#' @param timePointsSurvival A numeric vector object of the time points to be used.
#'   This should be the unique observed failure times.
#'   If 0 is not the first value, it will be concatenated by the software.
#'   If endPoint is "CR" then timePointsEndpoint is the same as timePointsSurvival.
#'
#' @param timePointsEndpoint A numeric vector object of the time points to be used.
#'   This should be the unique observed endpoint times.
#'   If endPoint is "RE" then timePointsEndpoint is a vector of the unique recurrent event times.
#'   If endPoint is "CR" then timePointsEndpoint can be the same as timePointsSurvival (since it is subset).
#'   If 0 is not the first value, it will be concatenated by the software.
#'
#' @param timePoints We recommend using the timePointsSurvival and timePointsEndpoint parameters.
#'   If both those and these are used, then those will overrule timePoints and nTimes parameters.
#'   A character object or a numeric vector object. If a character
#'   object, must be one of \{"quad", "uni", "exp"\} indicating the distribution
#'   from which the time points are to be calculated. For character input,
#'   input 'nTimes' must also be provided. If a numeric vector, the
#'   time points to be used. If 0 is not the first value, it will be
#'   concatenated by the software.
#'   We strongly recommend using the same number of decimal points between
#'   your failure times and the time points derived.
#'   For example, if T = 2.1, 2.2, 2.4, 5.1, then the time points here should
#'   also be to 1 decimal point.
#'   The largest time point should be the maximum of the observed failure time
#'   of your dataset.
#'   If you are using years and you have 2.1 total years, then you should use
#'   0.01 increments until you reach the maximum years (2.1).
#'   Default is NULL because we suggest inputting 'timePointsSurvival' and 'timePointsEndpoint'.
#'
#' @param nTimes An integer object. The total number of time points to be
#'   generated and considered. Used in conjunction with input 'timePoints'
#'   when 'timePoints' is a character; ignored otherwise.
#'   If inputting 'timePointsSurvival' and 'timePointsEndpoint' then leave 'nTimes' blank.
#'   Default is NULL because we recommend inputting 'timePointsSurvival' and 'timePointsEndpoint'.
#'
#' @param tau A numeric object or NULL. The study length. If NULL, the
#'   maximum timePoint is used. For RE, this is the maximum of all the stop times in
#'   the recurrent event dataset (multiple rows per person)
#'
#' @param criticalValue1 A character object. Must be one of
#'   \{"mean", "prob", "mean.prob.combo"\}. The estimator for the value
#'   of a treatment rule for Phase 1. For "mean": the mean survival/cif time; for
#'   "prob": the mean survival/cif probability at time 'evalTime';
#'   for "mean.prob.combo": first the mean survival/cif probability is used, if ties
#'   exist across treatments, the mean survival/cif time is used to identify the
#'   optimal. Typically, the same criticalValue should be used for Phase 1 and Phase 2.
#'   For "RE" endpoint, only "mean" is valid as the others are not yet coded.
#'
#' @param criticalValue2 A character object. Must be one of
#'   \{"mean", "prob", "mean.prob.combo"\}. The estimator for the value
#'   of a treatment rule for Phase 2. For "mean": the mean survival/cif time; for
#'   "prob": the mean survival/cif probability at time 'evalTime';
#'   for "mean.prob.combo": first the mean survival/cif probability is used, if ties
#'   exist across treatments, the mean survival/cif time is used to identify the
#'   optimal.Typically, the same criticalValue should be used for Phase 1 and Phase 2.
#'
#' @param tol1 Must be a vector where the length depends on criticalValue1. If criticalValue1
#'    is one of \{"mean", "prob", "area"\}, then length of tol1 is 2, where the first element
#'    is the tolerance specified for criticalValue under the curve for the treatment and
#'    counterfactual treatment for overall survival (Phase 1). The second element is for the
#'    specific scenario where one predicted criticalValue is 0 and the other is not.
#'    Then we take the absolute difference of the two and compare to the second tolerance criteria.
#'    When criticalValue1 is "mean.prob.combo", then we have a length of 4 elements, where the first
#'    two are for the ratio for mean and prob respectively, and the second two are for the difference
#'    for mean and prob, respectively. The defaults are 0.1 for mean ratio, 0.2 for prob ratio,
#'    0 for mean difference, and 0 for prob difference, respectively.
#'
#' @param evalTime A numeric object or NULL. This is T0. If numeric, the time at which
#'   the survival probability and/or cumulative incidence function is
#'   to be estimated to determine the optimal treatment rule;
#'   'criticalValue' must be one of  \{"mean.prob.combo", "prob"\}.
#'   If NULL, 'criticalValue' must be \{"mean"\}.
#'
#' @param splitRule1 A character object. splitRule for Phase 1
#'   Must be one of \{'logrank_surv', 'mean_surv'\} indicating the test used to determine an optimal split.
#' @param splitRule2 A character object. splitRule for Phase 2
#'   Must be one of \{'gray_cr', 'csh_cr', 'gray_re'\} indicating the test used to determine an optimal split.
#'
#'   If NULL and Phase = 'Survival' and 'criticalValue' = 'mean', then takes value 'mean'.
#'   If NULL and Phase = 'Survival' and 'criticalValue' = 'prob' or 'mean.prob.combo', then takes value 'logrank'.
#'   If NULL and Phase = 'CR' and endPoint = 'CR', takes value 'graycr'. # default over cshcr
#'   If NULL and Phase = 'RE' and endPoint = 'RE', takes value 'grayre'.
#'   ## survival endpoint ##
#'   logrank_surv = 1
#'   mean_surv = 2
#'   ## other endpoints (CR, RE) ##
#'   gray_cr = 3
#'   csh_cr = 4
#'   gray_re = 5
#'   if splitRule <=2 then Phase 1 survival
#'   if splitRule 3-4 then Phase 2 endpoint (CR)
#'   if splitRule >4 then Phase 2 endpoint (RE)
#'
#' @param ERT A logical object. If TRUE, the Extremely Randomized Trees
#'   algorithm is used to select the candidate variable.
#'
#' @param sampleSize A numeric object or NULL.
#'   The fraction (0 < sampleSize <= 1) of the data to be used for each
#'   tree in the forest. If NULL and 'ERT' is TRUE,
#'   sampleSize defaults to 1.0. If NULL and 'ERT'
#'   is FALSE, sampleSize defaults to 0.632.
#'
#' @param uniformSplit A logical object. If 'ERT' and 'uniformSplit' are TRUE,
#'   the random cutoff is sampled from a uniform distribution over the range
#'   of available covariate values. If 'ERT' is TRUE and 'uniformSplit' is
#'   FALSE, a case is randomly selected and the cutoff is taken to be the mean
#'   cutoff between it and the next largest covariate value. If 'ERT' is FALSE,
#'   input is ignored.
#'
#' @param replace A logical object or NULL. If TRUE, the sample drawn for each
#'   of the nTree trees may have duplicate records. If FALSE, no individual is
#'   present in the sample for than once. If NULL, 'replace' = !'ERT'.
#'
#' @param randomSplit A numeric object. The probability that a random split
#'   will occur. Must be 0 < randomSplit < 1.
#'
#' @param tieMethod A character object. Must be one of
#'   \{"first", "random"\}. If multiple splits lead to the same
#'   value, the method by which the tie is broken.
#'
#' @param minEvent An integer object. The minimum number of events that must be
#'   present in a node.
#'
#' @param nodeSize An integer object. The minimum number of individuals that
#'   must be present in a node.
#'
#' @param nTree An integer object. The number of trees to grow.
#'
#' @param mTry An integer or integer vector object. The maximum number of
#'   covariates to sample for each split. If a vector, each element
#'   corresponds to its respective decision point.
#'
#' @param pooled A logical object. If TRUE, data are pooled for the analysis.
#'    If FALSE, data is separated into groups based on treatment
#'    received and a tree is grown for each treatment group.
#'
#' @param stratifiedSplit A numeric object. The stratified random split
#'    coefficient. Covariates for which the number of splits (s_i) is less
#'    than s*stratifiedSplit/d are explored preferentially
#     (s is the total number of splits, d is the
#'    total number of covariates under consideration).
#'
#' @param stageLabel A character object. If using a common formula, the
#'    character used to separate the covariate from the decision point label.
#'    See details. IGNORE THIS FOR SINGLE STAGE.
#'
#' @references Zhou, C.W. and Kosorok, M.R.
#'   Estimating optimal individualized treatment regimes for survival data with competing risks.
#'   In prepraration.
#'
#'   Zhou, C.W. and Kosorok, M.R.
#'   Estimating optimal individualized treatment regimes for recurrent events and terminal event in survival data.
#'   In preparation.
#'
#' @include VerifyData.R VerifyTxName.R VerifyModels.R VerifySampleSize.R
#' @include VerifyEndPoint.R VerifyUsePrevTime.R class_Parameters.R
#' @include class_ITRSurvStep.R class_ITRSurv.R
#' @import methods
#' @export
#' @useDynLib itrSurv
#' @import survival
#'
#' @returns An S4 object of class ITRSurv containing the key results and
#'   input parameters of the analysis. The information contained therein
#'   should be accessed through convenience functions stage(), show(), print(),
#'   and predict().
#'
#' @seealso \code{\link{predict}} for retrieving the optimal treatment
#'    and/or the optimal survival curves. \code{\link{stage}} for retrieving stage
#'    results as a list. \code{\link{show}} for presenting the analysis results.
#'
#' @examples
#'
#'
#' dt <- data.frame("Y.1" = sample(1:100,100,TRUE),
#'                  "D.1" = rbinom(100, 1, 0.9),
#'                  "A.1" = rbinom(100, 1, 0.5),
#'                  "X.1" = rnorm(100))
#'
#' itrSurv(data = dt,
#'         txName = c("A.1"),
#'         models = list(Surv(Y.1,D.1)~X.1+A.1)
#'
#' # common formula
#' itrSurv(data = dt,
#'         txName = c("A.1"),
#'         models = Surv(Y,D)~X+A,
#'         usePrevTime = TRUE,
#'         stageLabel = ".")
#'
#' # common formula and pooled analysis
#' itrSurv(data = dt,
#'         txName = c("A.1"),
#'         models = Surv(Y,D)~X+A,
#'         stageLabel = ".",
#'         pooled = TRUE)
#'
#' dt <- data.frame("Y.1" = sample(1:100,100,TRUE), "Y.2" = sample(1:100,100,TRUE),
#'                  "D.1" = rbinom(100, 1, 0.9), "D.2" = rbinom(100,1,0.9),
#'                  "A.1" = rbinom(100, 1, 0.5), "A.2" = rbinom(100,1,0.5),
#'                  "X1" = rnorm(100), "X2" = rnorm(100))
#'
#' # common formula with only baseline covariates
#' itrSurv(data = dt,
#'         txName = c("A.1", "A.2"),
#'         models = Surv(Y,D)~X1+X2+A)
#'
#' # common formula with only baseline covariates
#' # cutoff selected from indices
#' itrSurv(data = dt,
#'         txName = c("A.1", "A.2"),
#'         models = Surv(Y,D)~X1+X2+A,
#'         ERT = TRUE, uniformSplit = FALSE)
#'
#' # common formula with only baseline covariates
#' # not extremely random trees
#' itrSurv(data = dt,
#'         txName = c("A.1", "A.2"),
#'         models = Surv(Y,D)~X1+X2+A,
#'         ERT = FALSE)
#'
#' # common formula with only baseline covariates
#' # survival probability
#' itrSurv(data = dt,
#'         txName = c("A.1", "A.2"),
#'         models = Surv(Y,D)~X1+X2+A,
#'         criticalValue = 'mean.prob.combo')
#'

#'
itrSurv <- function(data,
                    endPoint,
                    yName,
                    idName,
                    txName,
                    epName,
                    models,
                    ...,
                    timePointsSurvival,
                    timePointsEndpoint, # for CR, this is the same as timePointsSurvival
                    timePoints = NULL, #"quad",
                    nTimes = NULL, # 100L,
                    tau = NULL,
                    criticalValue1 = "mean", # delete later
                    criticalValue2 = "mean", # delete later
                    # criticalValue = "mean",
                    tol1 = c(0.1,0.1,0,0), # default is eps0_ratio = 0.1 and eps0_diff = 0;
                    evalTime = NULL,
                    splitRule1 = NULL,
                    splitRule2 = NULL,
                    ERT = TRUE,
                    uniformSplit = NULL,
                    sampleSize = NULL,
                    replace = NULL,
                    randomSplit = 0.2,
                    tieMethod = "random",
                    minEvent = 3L,
                    nodeSize = 6L,
                    nTree = 10L,
                    mTry = NULL,
                    pooled = FALSE,
                    stratifiedSplit = NULL,
                    stageLabel = ".") {

  #######################################################################################################
  #######################################################################################################
  #######################################################################################################
  #######################################################################################################
  #######################################################################################################
  # NEXTSTEPS: modify verifydata so that it checks for it being the correct format
  # if endPoint == "CR" then we need a delta AND a causeJ variable!!
  #######################################################################################################
  #######################################################################################################
  #######################################################################################################
  #######################################################################################################
  #######################################################################################################

  # ensure endPoint is one of {'CR', 'RE'}.
  # Methods return the original character possibly modified to be upper case.
  endPoint <- .VerifyEndPoint(endPoint = endPoint)

  # ensure that 'data' is provided as a data.frame or a matrix and does not
  # contain NaN values. If 'data' is appropriate, method returns a data.frame
  # object
  # Verify epName and Verify idName are ran for RE endpoint
  data_list <- .VerifyData(data = data,
                           yName = yName,
                           epName = epName,
                           endPoint = endPoint,
                           idName = idName)

  data_surv = data_list[[1]]
  data_ep = data_list[[2]]

  # ensure that 'txName' is provided as a character or character vector and
  # that the provided names are present in 'data'. If 'txName' is appropriate,
  # the object returned is the original input without modification.
  txName <- .VerifyTxName(txName = txName, data = data)

  # # total number of individuals in dataset
  nSamples <- nrow(data_surv)
  message("Dataset sample size: N = ", nSamples)

  # ignore nDP - leftover from multi-stage part that we DON'T do.
  # nDP = length(x = txName) # number of decision points in the analysis
  if (endPoint == "CR"){
    # need to check but in CR, it's used maybe to determine which model to use? (phase1 vs phase2)
    nDP <- 1
  } else{
    nDP <- 1
  }
  if (endPoint == "CR"){
    # print(models)
    message("Setting nCauses = 2 because there is one priority cause and everything else is lumped into Cause 2*")
    nCauses = 2 #there are always 2 causes because one is the priorty cause, and everything else is lumped into cause 2.
  } else{
    nCauses = 1
  }


  #######################################################################################################
  #######################################################################################################
  #######################################################################################################
  #######################################################################################################
  #######################################################################################################
  # NEXTSTEPS: modify verifymodels to include output for endpoint
  # i.e., within the script, set it so if endPoint == "CR" then response2, del2 ,and models2 = CR stuff needed for gray's test and CIF
  #######################################################################################################
  #######################################################################################################
  #######################################################################################################
  #######################################################################################################
  #######################################################################################################

  # ensure that 'models' is provided as a formula or a list of formula and
  # that the provided models can be generated by the data. If the input is
  # appropriate, the object returned is list containing
  #   "models" - the original input.
  models_surv <- .VerifyModels(models = models[[1]], #Surv(obs_time, D.0) ~ Z1
                           endPoint = endPoint,
                           nDP = nDP,
                           nCauses = nCauses,
                           data = data_surv,
                           txName = txName,
                           epName = epName,
                           stageLabel = stageLabel)

  if (endPoint == "CR"){
    #this line is redundant since earlier we specified
    # data=data_surv=data_ep for models_ep, but just to be safe :D
    models_ep <- .VerifyModels(models = models[[2]], #Surv(obs_time, D.1) ~ Z1
                               endPoint = endPoint,
                               nDP = nDP,
                               nCauses = nCauses,
                               data = data_surv, # same dataset, just different event indicator
                               txName = txName,
                               epName = epName,
                               stageLabel = stageLabel)
  } else if (endPoint == "RE"){
    models_ep <- .VerifyModels(models = models[[2]], # Surv(TStart, TStop, epName) ~ Z1
                               endPoint = endPoint,
                               nDP = nDP,
                               nCauses = nCauses,
                               data = data_ep,
                               txName = txName,
                               epName = epName,
                               stageLabel = stageLabel)
  } else{
    stop("itrSurv.R Line 419: endPoint is not specified. no models_ep created.")
  }

  response <- models_surv$response
  del <- models_surv$delta
  models <- models_surv$models

  response2 <- models_ep$response # 2-column vector for RE
  del1 <- models_ep$delta
  # for CR: del1 is priority event indicator
  # for RE: del1 is recurrent event indicator
  models1 <- models_ep$models


  # below is OLD code before adding in RE.
  # models0 <- .VerifyModels(models = models,
  #                     endPoint = endPoint,
  #                     nDP = nDP,
  #                     nCauses = nCauses,
  #                     data1 = data_surv,
  #                     data2 = data_ep,
  #                     txName = txName,
  #                     epName = epName,
  #                     stageLabel = stageLabel)
  # if (endPoint == "CR"){
  #   response <- models0$response
  #   del <- models0$delta[[1]] #overall failure
  #   models <- models0$models[[1]] #overall failure survival model
  #   # print(del)
  #   # print(models)
  #   if (nCauses > 0){
  #     message("nCauses = ", nCauses)
  #     del1 = models0$delta[[2]] #cause 1 failure
  #     models1 = models0$models[[2]] #cause 1 failure CR model
  #     # for (cause in 1:length(models0)){
  #     #   message("Cause ", cause)
  #     #   index = cause + 1
  #     #   d1 <- models0$delta[[index]]
  #     #   m1 <- models0$models[[index]]
  #     #   d1name = sprintf("del%s", cause)
  #     #   m1name = sprintf("models%s", cause)
  #     #   assign(d1name, d1)
  #     #   assign(m1name, m1)
  #     # }
  #   }
  # } else if (endPoint == "RE"){
  #
  # } else{
  #   message("ERROR: endPoint is NEITHER CR NOR RE --- ")
  #   response <- models0$response
  #   del <- models0$delta
  #   models <- models0$models
  # }

  # combine all inputs that regulate tree and specify analysis preferences
  # function returns a Parameters object
  # print(criticalValue1)

  if (endPoint == "CR"){
    if (!all(timePointsEndpoint == timePointsSurvival)){
      message("itrSurv.R Line 477: For CR endpoint, timePointsEndpoint should be the same as timePointsSurvival, NOT the subset for priority cause")
      print("-- changing timePointsEndpoint to be equal to inputted timePointsSurvival --")
      timePointsEndpoint = timePointsSurvival
    }
  }

  params <<- .parameters(endPoint = endPoint,
                        timePointsSurvival = timePointsSurvival,
                        timePointsEndpoint = timePointsEndpoint,
                        timePoints = timePoints,
                        tau = tau,
                        nTimes = nTimes,
                        response = response,
                        response_endpoint = response2,
                        nTree = nTree,
                        ERT = ERT,
                        uniformSplit = uniformSplit,
                        randomSplit = randomSplit,
                        splitRule1 = splitRule1,
                        splitRule2 = splitRule2,
                        replace = replace,
                        nodeSize = nodeSize,
                        minEvent = minEvent,
                        tieMethod = tieMethod,
                        criticalValue1 = criticalValue1,
                        criticalValue2 = criticalValue2,
                        survivalTime = evalTime,
                        endpointTime = evalTime, # need to change endpointTime to endpointTime in all scripts.
                        nSamples = nSamples,
                        pooled = pooled,
                        stratifiedSplit = stratifiedSplit)

  if (endPoint == "CR" & params@endpointparam@splitRule == "gray_cr"){
    # print("TESLFJKLSFJDKLSFJKLS")
    # ensure that 'epName' is provided as a character or character vector and
    # that the provided names are present in 'data'. For RE: this input defines the
    # dataset for the endpoint Phase analysis. For CR: this is needed for 'gray_cr' test.
    # If 'epName' is appropriate, the object returned is the original input without modification.
    epName <- .VerifyEpName(epName = epName, data = data, endPoint = endPoint)
    d0 = data
    sym(epName)
    cencode = 0 # this is non-negotiable
    d <- d0[order(response),]
    # Extract the column specified by epName
    ord_cause_dat <- d %>% dplyr::pull(!!sym(epName))
    censind <- ifelse(ord_cause_dat == cencode, 0, 1)
    # we do it this way JUST IN CASE there are more than 2 causes specified in data
    # priority cause is always first cause.
    if (is.factor(ord_cause_dat)) {
      uc <- table(ord_cause_dat[censind==1])
      uclab <- names(uc)[uc>0]
    } else {
      uclab <- sort(unique(ord_cause_dat[censind==1])) # as.numeric(names(uc)[uc>0])
    }
    causeind0 <- ifelse(ord_cause_dat==uclab[1],1,0)
    ord_causeind <<- 2*censind-causeind0
    ord_response <<- sort(response)
  } else{
    ord_causeind <- rep(0, nrow(data))
    ord_response <- rep(0, nrow(data))
  }

  print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
  print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
  print("Phase 1: Survival")
  print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
  print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
  # message("----------------Survival Parameters: Storing on Fortran Side----------------")
  # store basic information on fortran side

  ### Phase 1: SURVIVAL
  params1 = params@survivalparam
  # View(params1)
  crit1 <- .CriticalValueCriterion(params1, subcrit = criticalValue1)
  # View(crit1)
  tol1 = .VerifyTol1(tol1 = tol1, criticalValue = crit1)

  # retrieve index and fraction if survival type value estimator
  # set as 0's if mean estimator
  if (crit1 == "mean" | crit1 == "area") {
    ind1 = 0L
    frac1 = 0.0
  } else if (crit1 == "mean.prob.combo" || crit1 == "prob" || crit1 == "surv.prob") {
    # print(params1@sIndex)
    ind1 = params1@sIndex
    frac1 = params1@sFraction
  }
  # print('test2')

  # set basic parameter values in Fortran
  splitR_1 = ifelse(params1@splitRule == 'logrank_surv', 1, 2) # 1 for logrank; 2 for truncated mean
  res1 = .Fortran("setUpBasics",
                 t_nt = as.integer(x = .NTimes(object = params1)),
                 t_dt = as.double(x = .TimeDiff(object = params1)),
                 t_rs = as.double(x = params1@randomSplit),
                 t_ERT = as.integer(x = params1@ERT),
                 t_uniformSplit = as.integer(x = params1@uniformSplit),
                 t_nodeSize = as.integer(x = .NodeSize(object = params1)),
                 t_minEvent = as.integer(x = .MinEvent(object = params1)),
                 t_rule = as.integer(x = splitR_1),
                 t_sIndex = as.integer(x = ind1),
                 t_sFraction = as.double(x = frac1),
                 t_stratifiedSplit = as.double(x = params1@stratifiedSplit),
                 t_replace = as.integer(params1@replace),
                 PACKAGE = "itrSurv")
  # print('test9')
  # ensure that if given, sampleSize is 0 < sampleSize <= 1 and that
  # a value is provided for each decision point. If only 1 value is given,
  # it is assumed to be used for all decision points
  sampleSize1 <- .VerifySampleSize(sampleSize = sampleSize,
                                  ERT = params1@ERT,
                                  nDP = nDP
                                  )
  # message("sample size1: ", sampleSize1)

  # ensure that mTry is provided as a vector. At this point, there is
  # no verification of an appropriate value
  if (length(x = mTry) == 1L) {
    mTry <- rep(x = mTry, times = nDP)
  } else if (!is.null(x = mTry) && {length(x = mTry) != nDP}) {
    stop("if provided as vector, mTry must be provided for each dp",
         call. = FALSE)
  }

  # message("----------------End of Survival Storing on Fortran Side----------------")


  # Fitting generalized RSF
  # analysis (aka one stage analysis)
  # message("#######################################################################")
  message("RSF Starting: Phase 1: look at survival")
  # message("#######################################################################")
  # message("#######################################################################")
  ####################################################################################
  ####################################################################################
  phaseResults <- list()
  # print(params1@timeDiff)
  # message("...start of .itrSurvStep...")
  # print(models)
  # View(data)
  # View(params1)
  # print(txName[nDP]); print(mTry[nDP]); print(sampleSize1[nDP])
  Phase1Results <<- .itrSurvStep(Phase = "Survival",
                                eps0 = tol1,
                                model_surv = models_surv,
                                model_ep = models_ep,
                                # model_cause = NULL, #discontinued July 2024 after adding in RE endpoint
                                endPoint = endPoint,
                                data = data_list, # list of data_surv and data_ep
                                ord_causeind = ord_causeind, # only needed for Phase2 CR when using grays test
                                ord_response = ord_response, # only needed for Phase2 CR when using grays test
                                priorStep = NULL,
                                params = params1,
                                txName = txName[nDP],
                                mTry = mTry[nDP],
                                sampleSize = sampleSize1[nDP])
  # message("...end of .itrSurvStep...")
  # message("Phase1Results")
  # View(Phase1Results)
  assign("Phase1Results_survival", Phase1Results, envir = .GlobalEnv)
  phaseResults[[1]] <- Phase1Results

  print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
  print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
  print(sprintf("Phase 2: Endpoint: %s", endPoint))
  print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
  print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
  params2 = params@endpointparam
#   #   # we still want to run this for everyone- we just won't use the results for anyone who had ratio_stopping_ind = 1.
#   #   endPoint = "CR" #delete later

  # message("----------------EndPoint Parameters: Storing on Fortran Side----------------")
  ### Phase 2: ENDPOINT
  crit2 <- .CriticalValueCriterion(params2, subcrit = criticalValue2)
  # message("crit2:", crit2)
  if (endPoint == "CR" | endPoint == "RE"){
    if (crit2 == "mean" | crit2 == "area") {
      ind2 = 0L
      frac2 = 0.0
    } else if (crit2 == "mean.prob.combo" || crit2 == "prob" || crit2 == "cif.prob") {
      ind2 = params2@eIndex
      frac2 = params2@eFraction
    }
  }

  splitR_2 = ifelse(params2@splitRule == 'gray_cr', 3,
                    ifelse(params2@splitRule == 'csh_cr', 4,
                           ifelse(params2@splitRule == 'gray_re', 5, NA)))
  if (is.na(splitR_2)){stop("SPLIT RULE IS WRONG!!?")}
  res2 = .Fortran("setUpBasics",
                  t_nt = as.integer(x = .NTimes(object = params2)),
                  t_dt = as.double(x = .TimeDiff(object = params2)),
                  t_rs = as.double(x = params2@randomSplit),
                  t_ERT = as.integer(x = params2@ERT),
                  t_uniformSplit = as.integer(x = params2@uniformSplit),
                  t_nodeSize = as.integer(x = .NodeSize(object = params2)),
                  t_minEvent = as.integer(x = .MinEvent(object = params2)),
                  t_rule = as.integer(splitR_2),
                  t_sIndex = as.integer(x = ind2),
                  t_sFraction = as.double(x = frac2),
                  t_stratifiedSplit = as.double(x = params2@stratifiedSplit),
                  t_replace = as.integer(params2@replace),
                  PACKAGE = "itrSurv")

  # ensure that if given, sampleSize is 0 < sampleSize <= 1 and that
  # a value is provided for each decision point. If only 1 value is given,
  # it is assumed to be used for all decision points
  sampleSize2 <- .VerifySampleSize(sampleSize = sampleSize,
                                   ERT = params2@ERT,
                                   nDP = nDP)
  # message("sample size2: ", sampleSize2)

  # ensure that mTry is provided as a vector. At this point, there is
  # no verification of an appropriate value
  if (length(x = mTry) == 1L) {
    mTry <- rep(x = mTry, times = nDP)
  } else if (!is.null(x = mTry) && {length(x = mTry) != nDP}) {
    stop("if provided as vector, mTry must be provided for each dp",
         call. = FALSE)
  }
  # message("----------------End of EndPoint Storing on Fortran Side----------------")
  # message("#######################################################################")
  message(sprintf("RSF Starting: Phase 2: look at %s", endPoint))
  # message("#######################################################################")
#   if (endPoint == "CR"){
#     cause = 1 #TO DO: CHANGE THIS TO BECOME PRIORITY CAUSE LATER
#   # for (cause in 1:nCauses){
#     name = sprintf("Phase2Results_cause%s",cause);
#     # print(name)
#     # print("!!!!!!!!!!!!! models")
#     # print(models)
#     Phase2Results <<- .itrSurvStep(Phase = toString(endPoint),
#                                    eps0 = tol1, # don't really even need b/c dont use in Phase 2
#                                    model = models,
#                                    model_cause = get(sprintf("models%s", cause)),
#                                    endPoint = endPoint,
#                                    data = data_ep, # use Phase 2 data
#                                    params = params2,
#                                    txName = txName[nDP],
#                                    mTry = mTry[nDP],
#                                    sampleSize = sampleSize2[nDP])
#     # print("Phase2Results")
#     # print(name)
#     # View(Phase2Results)
#     assign(name, Phase2Results, envir = .GlobalEnv)
#     # print("assigned.")
#   # }
# }

  Phase2Results <<- .itrSurvStep(Phase = toString(endPoint),
                                 eps0 = tol1, # don't really even need b/c dont use in Phase 2
                                 model_surv = models_surv,
                                 model_ep = models_ep,
                                 # model_cause = get(sprintf("models%s", cause)), #discontinued July 2024 after adding in RE endpoint
                                 endPoint = endPoint,
                                 data = data_list,  # list of data_surv and data_ep
                                 ord_causeind = ord_causeind, # only needed for CR when using grays test
                                 ord_response = ord_response, # only needed for CR when using grays test
                                 params = params2,
                                 txName = txName[nDP],
                                 mTry = mTry[nDP],
                                 sampleSize = sampleSize2[nDP])
    assign("Phase2Results_endpoint", Phase2Results, envir = .GlobalEnv)

  phaseResults[[2]] <- Phase2Results
  p2 <<- Phase2Results
  p1 <<- Phase1Results
  # print(head(cbind(stop = p1@optimal@Ratio_Stopping_Ind,
  #                  p1_opttx = p1@optimal@optimalTx,
  #                  p2_opttx = p2@optimal@optimalTx,
  #                  opttx = (p1@optimal@Ratio_Stopping_Ind==1)*p1@optimal@optimalTx +
  #                    (p1@optimal@Ratio_Stopping_Ind==0)*p2@optimal@optimalTx),20))
  itrSurv_OptTx = (Phase1Results@optimal@Ratio_Stopping_Ind==1)*Phase1Results@optimal@optimalTx +
    (Phase1Results@optimal@Ratio_Stopping_Ind==0)*Phase2Results@optimal@optimalTx
  # print(head(itrSurv_OptTx,20))

  phaseResults[[3]] <- itrSurv_OptTx

  names(phaseResults) = c("SurvivalPhase1Results", "EndPointPhase2Results", "FinalOptimalTx_Recc")
  assign("phaseResults", phaseResults, envir = .GlobalEnv)

  # #
  # #   #######################################################################################################
  # #   #######################################################################################################
  # #   #######################################################################################################
  # #   #######################################################################################################
  # #   #######################################################################################################
  # #   # NEXTSTEPS: GO TO MEETING NOTES AND ADD IN VALUE FUNCTIONS AND SUMMARIES!!!
  # #   #######################################################################################################
  # #   #######################################################################################################
  # #   #######################################################################################################
  # #   #######################################################################################################
  # #   #######################################################################################################
  # #   # value obtained from the first stage analysis
  # this is in class_ITRSurvSTEP = need to make sure it uses MIN for CIF and max for CI or something

# print("itrSurv.R: Value Functions (line 607)")
print("")
print("")
message("--- itrSurv: Value Functions ---")
value1Train <- .meanValue(object = phaseResults,
                          Phase = "Survival")
value2Train <- mean(phaseResults[[1]]@optimal@Ratio_Stopping_Ind == 0)

if (endPoint == "CR" | endPoint == "RE"){
  # print(endPoint)
  # amongst those in value2Train (aka Phase1Results@optimal@Ratio_Stopping_Ind == 0), find mean CIF curve
  value3Train <- .meanValue(object = phaseResults,
                            Phase = toString(endPoint),
                            endPoint = endPoint)
} else{
  message("NULL value3Train b/c no endPoint")
  value3Train = NULL
}
valueTrain = list(value1Train, list(value2Train), value3Train)
names(valueTrain[[2]]) = "PropPhase2"
message("Estimated Value:", appendLF = FALSE)
for (i in 1L:length(valueTrain)) {
  v = valueTrain[[i]]
  for (v1 in 1L:length(v)){
    message(" ", names(v)[v1], ": ", v[[ v1]], appendLF = FALSE)
  }
}
# store values in call structure for returned object
cl <- match.call()
cl[[ 1L ]] <- as.name("itrSurv")

  itr_results <- new(Class = "ITRSurv",
                    "phaseResults" = phaseResults,
                    "value" = valueTrain,
                    "call" = cl,
                    "params" = params)
  # View(itr_results)

  names(itr_results@value) = c("V1", "V2", "V3")

  return( itr_results )
  # return(phaseResults)
}


#' Hidden methods
#'
#' @name itrSurv-internal-api
#' @keywords internal
#' @import methods
NULL
