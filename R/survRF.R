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
#' @include class_CriticalValue.R class_CriticalValueMean.R
#' @include class_CriticalValueSurvival.R class_CriticalValueSurvivalatT0.R
#' #' @include class_CriticalValueCR.R
#' @include class_SurvRF.R
#' @import parallel
#'
# result <- .survRF(x = x[elig,,drop=FALSE],
# y = response[elig],
# pr = pr,
# delta = delta[elig],
# delta_cause = delta_cause[elig],
# params = params,
# mTry = mTry,
# txLevels = txLevels,
# model = mod,
# sampleSize = sampleSize)

.survRF <- function(..., Phase, eps0, x, delta, delta_cause, pr, params, mTry, sampleSize) {
  # message("starting .survRF from survRF.R")

  # if x_i is an unordered factor, nCat_i is the number of levels
  # if x_i is not a factor, nCat is 0
  # if x_i is an ordered factor, nCat is 1
  xLevels <- lapply(X = x, FUN = levels)
  # View(xLevels)

  nCat <- sapply(X = x, FUN = nlevels)
  # View(nCat)

  nCat <- ifelse(test = sapply(X = x, FUN = is.ordered), yes = 1L, no = nCat)
  # View(nCat)

  # number of individuals in training data
  nSamples <- nrow(x = x)
  # message('number of individuals in training data: ', nSamples)

  # number of time points
  nTimes <- nrow(x = pr)
  # message('number of time points: ', nTimes)

  # total number of trees to be grown in the forest
  # .NTree() is a getter method defined for Parameters objects
  nTree <- .NTree(object = params)
  # message("total number of trees to be grown in forest: ", nTree)

  # determine the number of samples to include in each tree
  sampleSize <- ceiling(x = sampleSize * nSamples)
  # message("number of samples to include in each tree: ", sampleSize)

  # maximum number of nodes in a tree
  maxNodes <- 2L * sampleSize + 1L
  # message("maximum nodes in a tree: ", maxNodes)

  # convert factors to integers
  x = data.matrix(frame = x)
  nr = nrow(x = x)
  # message("nr: ", nr)

  # message("setUpInners: Send info to Fortran")
  # send step specific x, pr, delta, mTry, nCat to Fortran
  res = .Fortran("setUpInners",
                 t_n = as.integer(x = nrow(x = x)),
                 t_np = as.integer(x = ncol(x = x)),
                 t_x = as.double(x = x),
                 t_pr = as.double(x = t(x = pr)),
                 t_delta = as.integer(x = delta),
                 t_delta_m = as.integer(x = delta_cause),
                 t_mTry = as.integer(x = mTry),
                 t_nCat = as.integer(x = nCat),
                 t_sampleSize = as.integer(x = sampleSize),
                 t_nTree = as.integer(x = params@nTree),
                 t_nrNodes = as.integer(x = maxNodes),
                 PACKAGE = "itrSurv")

  # if (sampleSize == 339){
  #   res339 <<- res
  #   trees2_339 <<- trees
  #   # View(res339)
  # }
  # if (sampleSize == 94){
  #   res94 <<- res
  #   trees2_94 <<- trees
  #   # View(res94)
  # }

  # message("================= Phase: ", Phase)
  if (grepl("surv", Phase, ignore.case = TRUE)){
    res_pooled0_surv <<- res
    # message("survTree: survTree in Fortran")
    Tree <- .Fortran("survTree",
                       forestSurvFunc = as.double(numeric(nTimes*nr)),
                       forestMean = as.double(numeric(nr)),
                       forestSurvProb = as.double(numeric(nr)),
                       PACKAGE = "itrSurv")

    # print(nr)
    Tree_Surv <<- Tree

  } else if (grepl("CR", Phase, ignore.case = TRUE)){
    res_pooled0_cif <<- res
    Tree <- .Fortran("cifTree",
                         forestSurvFunc = as.double(numeric(nTimes*nr)), #survival function averaged over forest
                         forestMean = as.double(numeric(nr)), #mean survival time averaged over forest
                         forestSurvProb = as.double(numeric(nr)), #survival probability at t0 averaged over forest
                         PACKAGE = "itrSurv")
    # print(nr)
    Tree_Cif <<- Tree
  }
  # message(sprintf("%s Tree for Phase %s", Phase, Phase))

  # retrieve trees from Fortran
  # message("survRF.R: Retrieving trees from Fortran")
  trees <- list()
  # print(nr)
  for (i in 1L:nTree) {

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

  # message("Done with .survRF")

}
