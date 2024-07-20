#' Individualized Treatment Regime for Survival Analysis
#'
#' Provides methods for estimating single stage optimal individualized treatment
#'   regimes for survival outcomes with multiple endpoints. Currently only applicable to single stage diseases.
#'   1) Competing risk with priority cause.
#'   2) Recurrent events.
#'
#' @param ... Ignored. Present only to require named inputs.
#'
#' @param data A data.frame object. The full dataset including treatments
#'   received, all stage covariates, observed times, and censoring
#'   indicators.
#'   Can be provided as a matrix object if column headers are included.
#'   Can contain missing data coded as NA, but cannot contain NaN.
#'
#' @param txName A character vector object. The treatment variable name for
#'   each decision point. Each element corresponds to the respective decision
#'   point (element 1 = 1st decision; element 2 = 2nd decision, etc.).
#'
#' @param models A single formula defining the response as a Surv() object
#'   and the covariate structure of the model.
#'   Note that this model should not include any terms of order > 1.
#'
#' @param endPoint A character object. Must be one of
#'   \{"CR", "RE"\}. The label for the ultimate end point type.
#'   For "CR": competing risks data to examine cumulative incidence function;
#'   For "RE": recurrent events data to examine mean frequency function
#'
#' @param timePoints A character object or a numeric vector object. If a character
#'   object, must be one of \{"quad", "uni", "exp"\} indicating the distribution
#'   from which the time points are to be calculated. For character input,
#'   input 'nTimes' must also be provided. If a numeric vector, the
#'   time points to be used. If 0 is not the first value, it will be
#'   concatenated by the software.
#'
#' @param nTimes An integer object. The total number of time points to be
#'   generated and considered. Used in conjunction with input 'timePoints'
#'   when 'timePoints' is a character; ignored otherwise.
#'
#' @param tau A numeric object or NULL. The study length. If NULL, the
#'   maximum timePoint is used.
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
#' @param splitRule A character object OR NULL.
#'   Must be one of \{"logrank", "mean", "logrankcr", "meancr"\}
#'   indicating the test used to determine an optimal split. If NULL and
#'   'criticalValue' = 'mean', takes value 'mean'. If NULL and
#'   'criticalValue' = 'prob' or 'mean.prob.combo', takes value 'logrank'.
#'   If NULL and endPoint = 'CR', takes value 'logrankCR'.
#'   logrank = 1
#'   mean = 2
#'   logrankcr = 3
#'   meancr = 4 but this is same as mean (=2) which we change later in fortran
#'   if splitRule <=2 then Phase 1 survival
#'   if splitRule >=3 then Phase 2 endpoint (CR)
#'
#' @param ERT A logical object. If TRUE, the Extremely Randomized Trees
#'   algorithm is used to select the candidate variable.
#'
#' @param sampleSize A numeric object, numeric vector object, or NULL.
#'   The fraction (0 < sampleSize <= 1) of the data to be used for each
#'   tree in the forest. If only
#'   one value is given, it is assumed to be the fraction for all decision
#'   points. If a vector is given, the length must be equal to the total
#'   number of decision points and each element corresponds to its respective
#'   decision point. If NULL and 'ERT' is TRUE,
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
#'    See details.
#'
#' @references Zhou, C.W. and Kosorok, M.R.
#'   Estimating optimal individualized treatment regimes for survival data with competing risks.
#'   In prepraration.
#'
#'   Zhou, C.W. and Kosorok, M.R.
#'   Estimating optimal individualized treatment regimes for recurrent events and terminal event in survival data.
#'   In preparation.
#'
#'
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
                    txName,
                    models,
                    endPoint = "CR",
                    ...,
                    timePoints = "quad",
                    nTimes = 100L,
                    tau = NULL,
                    criticalValue1 = "mean", # delete later
                    criticalValue2 = "mean", # delete later
                    # criticalValue = "mean",
                    tol1 = c(0.1,0.1,0,0), # default is eps0_ratio = 0.1 and eps0_diff = 0;
                    evalTime = NULL,
                    splitRule = NULL,
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
  # ensure that 'data' is provided as a data.frame or a matrix and does not
  # contain NaN values. If 'data' is appropriate, method returns a data.frame
  # object
  data <- .VerifyData(data = data, endPoint = endPoint)

  # total number of individuals in dataset
  nSamples <- nrow(x = data)

  # ensure that 'txName' is provided as a character or character vector and
  # that the provided names are present in 'data'. This input defines the
  # number of decision points for the analysis. If 'txName' is appropriate,
  # the object returned is the original input without modification.
  txName <- .VerifyTxName(txName = txName, data = data)

  # # # number of decision points in the analysis
  nDP <- length(x = txName)
  if (endPoint == "CR"){
    # print(models)
    nCauses = 2 #TO DO: make this general later!!
    # print("nCauses")
    # print(nCauses)
  } else{
    nCauses = 1
  }

  # ensure endPoint is one of {'CR', 'RE', 'MC'}.
  # Methods return the original character possibly modified to be upper case.
  endPoint <- .VerifyEndPoint(endPoint = endPoint)

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
  models0 <- .VerifyModels(models = models,
                          endPoint = endPoint,
                          nDP = nDP,
                          nCauses = nCauses,
                          data = data,
                          txName = txName,
                          stageLabel = stageLabel)

  if (endPoint == "CR"){
    response <- models0$response
    del <- models0$delta[[1]] #overall failure
    models <- models0$models[[1]] #overall failure survival model
    # print(del)
    # print(models)
    if (nCauses > 0){
      message("nCauses = ", nCauses)
      del1 = models0$delta[[2]] #cause 1 failure
      models1 = models0$models[[2]] #cause 1 failure CR model
      # for (cause in 1:length(models0)){
      #   message("Cause ", cause)
      #   index = cause + 1
      #   d1 <- models0$delta[[index]]
      #   m1 <- models0$models[[index]]
      #   d1name = sprintf("del%s", cause)
      #   m1name = sprintf("models%s", cause)
      #   assign(d1name, d1)
      #   assign(m1name, m1)
      # }
    }
  } else{
    response <- models0$response
    del <- models0$delta
    models <- models0$models
  }

  # combine all inputs that regulate tree and specify analysis preferences
  # function returns a Parameters object
  print(criticalValue1)
  params <- .parameters(endPoint = endPoint,
                        timePoints = timePoints,
                        tau = tau,
                        nTimes = nTimes,
                        response = response,
                        nTree = nTree,
                        ERT = ERT,
                        uniformSplit = uniformSplit,
                        randomSplit = randomSplit,
                        splitRule = splitRule,
                        replace = replace,
                        nodeSize = nodeSize,
                        minEvent = minEvent,
                        tieMethod = tieMethod,
                        criticalValue1 = criticalValue1,
                        criticalValue2 = criticalValue2,
                        survivalTime = evalTime,
                        CIFTime = evalTime,
                        nSamples = nSamples,
                        pooled = pooled,
                        stratifiedSplit = stratifiedSplit)
  # message("End of params.")

  # print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
  # print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
  print("Phase 1: Survival")
  # print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
  # print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
  # message("----------------Survival Parameters: Storing on Fortran Side----------------")
  # store basic information on fortran side

  ### Phase 1: SURVIVAL
  # retrieve index and fraction if survival type value estimator
  # set as 0's if mean estimator
  params1 = params@survivalparam
  # View(params1)
  crit1 <- .CriticalValueCriterion(params1, subcrit = criticalValue1)
  # View(crit1)
  tol1 = .VerifyTol1(tol1 = tol1, criticalValue = crit1)

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
  if (params1@splitRule == 'logrank'){
    #logrank for survival
    splitR = as.integer(params1@splitRule == 'logrank')
  } else{
    # truncated mean test for survival
    splitR = as.integer(2)
  }
  # message("FJDSD;JF: as.integer(x = .NTimes(object = params1)):", as.integer(.NTimes(object = params1)))
  # print(.NTimes(object = params1))
  # View(params1)
  res1 = .Fortran("setUpBasics",
                 t_nt = as.integer(x = .NTimes(object = params1)),
                 t_dt = as.double(x = .TimeDiff(object = params1)),
                 t_rs = as.double(x = params1@randomSplit),
                 t_ERT = as.integer(x = params1@ERT),
                 t_uniformSplit = as.integer(x = params1@uniformSplit),
                 t_nodeSize = as.integer(x = .NodeSize(object = params1)),
                 t_minEvent = as.integer(x = .MinEvent(object = params1)),
                 t_rule = splitR, # as.integer(x = params1@splitRule == 'logrank')
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
                                  nDP = nDP)
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
                                model = models,
                                model_cause = NULL,
                                endPoint = endPoint,
                                data = data,
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


  # print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
  # print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
  print(sprintf("Phase 2: Endpoint: %s", endPoint))
  # print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
  # print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
  params2 = params@endpointparam
#   #   # we still want to run this for everyone- we just won't use the results for anyone who had ratio_stopping_ind = 1.
#   #   endPoint = "CR" #delete later

  # message("----------------EndPoint Parameters: Storing on Fortran Side----------------")
  ### Phase 2: ENDPOINT
  crit2 <- .CriticalValueCriterion(params2, subcrit = criticalValue2)
  # message("crit2:", crit2)
  if (endPoint == "CR"){
    if (crit2 == "mean" | crit2 == "area") {
      ind2 = 0L
      frac2 = 0.0
    } else if (crit2 == "mean.prob.combo" || crit2 == "prob" || crit2 == "cif.prob") {
      ind2 = params2@cIndex
      frac2 = params2@cFraction
    }
  }

  # set basic parameter values in Fortran
  if (params2@splitRule == 'logrankcr'){
    splitR_endpoint = as.integer(3)
  } else if (params2@splitRule == 'meancr' | params2@splitRule == 'mean'){
    splitR_endpoint = as.integer(4)
  } else{
    message("splitRule is unclear and we assigned splitR_endpoint as 99")
    splitR_endpoint = as.integer(99)
  }

  # print(splitR_endpoint)
  res2 = .Fortran("setUpBasics",
                  t_nt = as.integer(x = .NTimes(object = params2)),
                  t_dt = as.double(x = .TimeDiff(object = params2)),
                  t_rs = as.double(x = params2@randomSplit),
                  t_ERT = as.integer(x = params2@ERT),
                  t_uniformSplit = as.integer(x = params2@uniformSplit),
                  t_nodeSize = as.integer(x = .NodeSize(object = params2)),
                  t_minEvent = as.integer(x = .MinEvent(object = params2)),
                  t_rule = splitR_endpoint, # as.integer(x = params2@splitRule == 'logrank')
                  t_sIndex = as.integer(x = ind2),
                  t_sFraction = as.double(x = frac2),
                  t_stratifiedSplit = as.double(x = params2@stratifiedSplit),
                  t_replace = as.integer(params2@replace),
                  PACKAGE = "itrSurv")
# #
  # message("res2:")
  # View(res2)

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
  if (endPoint == "CR"){
    cause = 1 #TO DO: CHANGE THIS TO BECOME PRIORITY CAUSE LATER
  # for (cause in 1:nCauses){
    name = sprintf("Phase2Results_cause%s",cause);
    # print(name)
    # print("!!!!!!!!!!!!! models")
    # print(models)
    Phase2Results <<- .itrSurvStep(Phase = toString(endPoint),
                                   eps0 = tol1, # don't really even need b/c dont use in Phase 2
                                   model = models,
                                   model_cause = get(sprintf("models%s", cause)),
                                   endPoint = endPoint,
                                   data = data,
                                   params = params2,
                                   txName = txName[nDP],
                                   mTry = mTry[nDP],
                                   sampleSize = sampleSize2[nDP])
    # print("Phase2Results")
    # print(name)
    # View(Phase2Results)
    assign(name, Phase2Results, envir = .GlobalEnv)
    # print("assigned.")
  # }
}

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
message("--- itrSurv: Value Functions ---")
value1Train <- .meanValue(object = phaseResults,
                          Phase = "Survival")
value2Train <- mean(phaseResults[[1]]@optimal@Ratio_Stopping_Ind == 0)

if (endPoint == "CR"){
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

