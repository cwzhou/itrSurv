# Internal function to generate forest
#
# Function is not exported
#
# @param x A data.frame object {n x p}. The model.frame for covariates to be
#   considered in splitting
#
# @param delta An integer vector object {n}. Vector of indicators for censoring
#   (1 = not censored; 0 = censored)
#
# @param pr A matrix object {nt x n}. Probability mass vector of survival
#   function
#
# @param pr2 A matrix object {nt x n}. Used to calculate at-risk subjects for RE
#   which has multiple records per subject
#
# @param pr2_surv A matrix object {nt_surv x n}. Used to calculate at-risk subjects for death during RE
#   which has one record per subject
#
# @param ord_causeind A response-ordered vector of cause status (for CR). Needed for
#   Gray's Test for node-splitting. 0 = censored, 1 = priority cause, 2 = any other causes.
#   Vector of 0s for RE or CR where test is not Gray's test.
#
# @param ord_response An ordered vector of response (for CR). Needed for
#   Gray's Test for node-splitting. Vector of 0s for RE or CR where test is not Gray's test.
#
# @param params A Parameters object. All information that regulates tree and
#   specifies analysis preferences.
#
# @param mTry An integer object. The maximum number of covariates to use
#   for splitting algorithm.
#
# @params sampleSize An integer object. The number of samples to draw for each
#    tree
#
#' @include class_Parameters.R
#' @include class_CriticalValue.R
#' @include class_CriticalValueArea.R class_CriticalValueMean.R
#' @include class_CriticalValueSurv.R
#' @include class_CriticalValueEndpoint.R
#' @include class_SurvRF.R
#' @import parallel
#'
# result <- .survRF(x = x[elig,,drop=FALSE],
# y = response[elig],
# pr = pr,
# delta = delta[elig],
# delta_endpoint = delta_endpoint[elig],
# params = params,
# mTry = mTry,
# txLevels = txLevels,
# model = mod,
# sampleSize = sampleSize)

.survRF <- function(..., Phase, eps0, x, idvec,
                    pr, pr2, pr2_surv = NULL, pr_surv = NULL,
                    ord_causeind, ord_response,
                    delta, delta_endpoint,
                    person_indicator,
                    params, mTry, sampleSize) {
  # message("---------- starting .survRF function from survRF.R ----------------")

  sampleSize_frac = sampleSize
  message("sampleSize_frac", sampleSize_frac)

  # if x_i is an unordered factor, nCat_i is the number of levels
  # if x_i is not a factor, nCat is 0
  # if x_i is an ordered factor, nCat is 1
  xLevels <- lapply(X = x, FUN = levels)
  # View(xLevels)

  nCat <- sapply(X = x, FUN = nlevels)
  # View(nCat)

  nCat <- ifelse(test = sapply(X = x, FUN = is.ordered), yes = 1L, no = nCat)
  # View(nCat)

  # FIGURE OUT HOW TO DO WITHOUT CRASHING
  # number of individuals in training data
  nSamples <- nrow(x = x) # RE: nrow(x) is for RE for Phase 2 data (# records)
  print(nSamples)
  if (Phase == "RE"){
    nSamples_surv = length(unique(idvec))
  } else{
    nSamples_surv = nSamples
  }
  # message('number of records in training data: ', nSamples)
  # message('number of individuals in training data: ', nSamples_surv)

  # tpSurv =
  # View(pr)
  # View(pr2)
  # message("dim pr:", dim(pr))
  # if (Phase == "RE"){
  #   message("dim pr2:", dim(pr2))
  # }

  # number of time points
  nTimes <- nrow(x = pr)
  if (is.null(nTimes)){
    nTimes = length(x = pr)
  }
  # message('number of time points: ', nTimes)

  # total number of trees to be grown in the forest
  # .NTree() is a getter method defined for Parameters objects
  nTree <- .NTree(object = params)
  # message("total number of trees to be grown in forest: ", nTree)

  # determine the number of samples to include in each tree
  sampleSize <- ceiling(x = sampleSize_frac * nSamples)
  # message("number of samples to include in each tree: ", sampleSize)
  # for RE: this is # records

  # maximum number of nodes in a tree
  maxNodes <- 2L * sampleSize + 1L
  # message("maximum nodes in a tree based on samples: ", maxNodes)

  if (Phase == "RE"){
    sampleSize_surv <- ceiling(x = sampleSize_frac * nSamples_surv)
    maxNodes_surv <- 2L * sampleSize_surv + 1L
    # message("number of individuals to include in each tree: ", sampleSize_surv)
    # message("maximum nodes in a tree based on people: ", maxNodes_surv)
  } else{
    sampleSize_surv <- sampleSize
    maxNodes_surv = maxNodes
  }

  # convert factors to integers
  x = data.matrix(frame = x)
  nr = nrow(x = x) # number of individuals based on covariate length
  # message("number of individuals, nr: ", nr)

  # message("setUpInners: Send info to Fortran")
  # send step specific x, pr, delta, mTry, nCat to Fortran

  dd <<- delta
  oo <<- ord_causeind
  rr <<- ord_response
  ii <<- idvec


  # if (Phase == "RE"){
  #   # View(tSurv)
  #   # View(tSurv_surv)
  #   # View(t(pr))
  #   # View(t(pr_surv))
  #   stop("Testing RE Phase.")
  # }

    res = .Fortran("setUpInners",
                 t_n = as.integer(x = nSamples), # number of subjects for Phase1/2CR, number of records for Phase2RE
                 t_n_surv = as.integer(x = nSamples_surv), # number of subjects
                 t_idvec = as.integer(x = idvec), # id labels (1 row per person for Phase1/Phase2CR, multiple rows per person for Phase2RE to later obtain pr2 subset for at risk for death in mff in Fortran)
                 t_person_ind = as.integer(x = person_indicator), # needed for 2RE
                 t_np = as.integer(x = ncol(x = x)), # number of covariates
                 t_x = as.double(x = x), # covariates
                 t_pr = as.double(x = t(x = pr)), # transpose(pr): dim: n x nt #used to get number of events
                 t_pr2 = as.double(x = t(x = pr2)), # pr2 to get at-risk for RE during isPhase2RE
                 t_pr2surv = as.double(x = t(x = pr2_surv)), #to get at-risk for death during isPhase2RE
                 t_prsurv = as.double(x = t(x = pr_surv)), # transpose(prsurv): dim: n_records x nt_death #used to get number of events
                 t_ord_causeind = as.integer(x = ord_causeind), # CR: response-ordered cause status #RE: vec of 0s
                 t_ord_response = as.double(x = ord_response), # CR: ordered response #RE: vec of 0s
                 t_delta = as.integer(x = delta), # delta failure
                 t_delta_m = as.integer(x = delta_endpoint), # endpoint delta
                 t_mTry = as.integer(x = mTry),
                 t_nCat = as.integer(x = nCat),
                 t_sampleSize = as.integer(x = sampleSize),
                 t_sampleSize_surv = as.integer(x = sampleSize_surv),
                 t_nTree = as.integer(x = params@nTree),
                 t_nrNodes = as.integer(x = maxNodes), # for 2nd endpoint
                 t_nrNodes_surv = as.integer(x = maxNodes_surv), # for survival
                 PACKAGE = "itrSurv")

  #message(" dfgdfgd ================= Phase: ", Phase)
  if (grepl("surv", Phase, ignore.case = TRUE) | Phase == 1){
    res_pooled0_surv <<- res
    # print("survTree: survTree in Fortran")
    # nTimes = maximum number of time points
    # nr = number of rows
    # message("number of rows/people in dataset: nr = ", nr)
    # message("maximum number of time points: nTimes = ", nTimes)
    Tree <- .Fortran("survTree",
                       forestSurvFunc = as.double(numeric(nTimes*nr)),#survival function averaged over forest
                       forestMean = as.double(numeric(nr)),#mean survival time averaged over forest
                       forestSurvProb = as.double(numeric(nr)),#survival probability at t0 (evalTime/survivalTime) averaged over forest
                       PACKAGE = "itrSurv")

    # print(nr)
    Tree_Surv <<- Tree

  } else{ #if (grepl("CR", Phase, ignore.case = TRUE)){
    res_pooled0_cif <<- res
    Tree <- .Fortran("cifTree",
                         forestSurvFunc = as.double(numeric(nTimes*nr)),
                         forestMean = as.double(numeric(nr)),
                         forestSurvProb = as.double(numeric(nr)),
                         PACKAGE = "itrSurv")
    # print(nr)
    Tree_Cif <<- Tree
  }
  message(sprintf("%s Tree for Phase %s", Phase, Phase))

  # retrieve trees from Fortran
  # message("survRF.R: Retrieving trees from Fortran")
  trees <- list()
  # print(nr)
  for (i in 1L:nTree) {
    # print(sprintf("survRF.R: Tree %s", i))

    trees[[ i ]] <- list()
    # if (i == 2 & nSamples == 339){
    #   tree2list_339 <<- trees
    #   # View(tree2list_339[[1]])
    # }
    # if (i == 2 & nSamples == 94){
    #   tree2list_94 <<- trees
    #   # View(tree2list_94[[1]])
    # }

    treeDims <- .Fortran("treeDim",
                         iTree = as.integer(x = i),
                         nr = as.integer(x = 1L),
                         nc = as.integer(x = 1L),
                         PACKAGE = "itrSurv")

    # print("treeDims")
    # # print(treeDims)
    # print(treeDims$nr)

    temp <- .Fortran("getTree",
                     iTree = as.integer(x = i),
                     nr = as.integer(x = treeDims$nr),
                     nc = as.integer(x = treeDims$nc),
                     nodes = as.double(x = numeric(length = treeDims$nc*treeDims$nr)),
                     Func = as.double(x = numeric(length = treeDims$nr*nTimes)),
                     mean = as.double(x = numeric(length = treeDims$nr)),
                     Prob = as.double(x = numeric(length = treeDims$nr)),
                     PACKAGE = "itrSurv")

    gettree <<-temp
    trees[[ i ]]$nodes <- matrix(data = temp$nodes,
                                 nrow = treeDims$nr,
                                 ncol = treeDims$nc)
    trees[[ i ]]$Func <- matrix(data = temp$Func,
                                    nrow = nTimes,
                                    ncol = treeDims$nr)
    trees[[ i ]]$mean <- temp$mean
    trees[[ i ]]$Prob <- temp$Prob
    if (i == 1){
      tree1list <<- trees
      # View(tree1list[[i]])
    }
  }
  forest <- list()
  forest[[ "Func" ]] <- matrix(data = Tree$forestSurvFunc,
                                   nrow = nTimes, ncol = nr)
  forest[[ "mean" ]] <- Tree$forestMean

  # print(Phase)
  crit <- .CriticalValueCriterion(params)
  if (crit %in% c("prob", "mean.prob.combo")) {
    forest[[ "Prob" ]] <- Tree$forestSurvProb
  }

  return( new(Class = "SurvRF",
              "trees" = trees,
              "forest" = forest,
              "variables" = colnames(x = x),
              "mTry" = mTry,
              "nCat" = nCat,
              "xLevels" = xLevels) )

  message("Done with .survRF")

}
