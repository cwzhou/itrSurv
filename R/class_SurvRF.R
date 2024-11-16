# Virtual class to denote objects are arising from survRF step
#
# Methods
#   .Predict(object, newdata, ...) {new; not allowed}
setClass(Class = "SurvRFObject", # Defines a virtual class named "SurvRFObject," which indicates that it's a class for which methods will be defined in subclasses.
         contains = c("VIRTUAL"))

#-------------------------------------------------------------------------------
# method to make predictions for new data or to retrieve fitted values
#-------------------------------------------------------------------------------
setGeneric(name = ".Predict", # Define a generic function named ".Predict" for making predictions on new data or retrieving fitted values.
           def = function(object, newdata, ...) {
             standardGeneric(".Predict")
             })

setMethod(f = ".Predict", # method cannot be used on object of class "SurvRFObject". Specific subclasses should implement their own ".Predict" methods.
          signature = c(object = "ANY",
                        newdata = "ANY"),
          definition = function(object, newdata, ...) {
            View(object)
            View(newdata)
            stop("class_SurvRF Line 18: not allowed")
            })

#-------------------------------------------------------------------------------
# method to make predictions for new data at each tx level
#-------------------------------------------------------------------------------
setGeneric(name = ".PredictAll", #Define a generic function named ".PredictAll" for making predictions on new data at each treatment level.
           def = function(object, ...) {
                   standardGeneric(".PredictAll")
               # message("class_SurvRF.R: .PredictAll from LINE 30")
                 })

setMethod(f = ".PredictAll", # method is defined as "not allowed" for objects of class "SurvRFObject." Specific subclasses should provide their own implementations for this method.
          signature = c(object = "ANY"),
          definition = function(object, ...) {
            stop("not allowed")
            })

# Class for storing survRF results for pooled analysis
#
# Class is not exported and is for internal convenience only
#
#  @slot trees A list object. The results of the tree building algorithm for
#    each tree in the forest
#
#  @slot forest A list object. The values averaged across all trees in the
#    forest
#
#  @slot variable A character vector. The variables considered in the
#    analysis
#
#  @slot mTry An integer. The maximum number of covariates considered for
#    splitting
#
#  @slot nCat An integer vector. The number of categories for each covariate
#    considered. >=2 unordered factor, 1 ordered factor, 0 continuous
#
#  @slot xLevels A list object. The categories in each covariate considered.
#
# Methods
#   .Predict(object, newdata, ...) {defined}
#
setClass(Class = "SurvRF",
         slots = c("trees" = "list",
                   "forest" = "list",
                   "variables" = "character",
                   "mTry" = "integer",
                   "nCat" = "integer",
                   "xLevels" = "list"),
         contains = c("SurvRFObject"))

# # Extracting and organizing tree-related results from a random forest survival analysis, considering both stratified and non-stratified scenarios.
# # It provides a structured way to access the tree and forest information generated during the analysis.
# .stageSurvRF <- function(object) {
#   message('.stageSurvRF')
#   res <- list() # The function creates an empty list res to hold the extracted information.
#   # Tree Extraction for Stratified Analysis:
#   if (is(object = object, class2 = "SurvRFStratified")) { # The function checks whether the input object is of class SurvRFStratified using the is() function. If it is, it suggests that the input object contains stratified survival random forest results. Stratified analysis is often used to analyze different subgroups or strata of data separately.
#     for (i in 1L:length(x = object@strat)) { # If the input object is of class SurvRFStratified, the function iterates through the strata (groups) and extracts tree-related information for each stratum.
#       res[[ i ]] <- list() # For each stratum, it creates a list that contains:
#       res[[ i ]][[ "trees" ]] <- object@strat[[ i ]]@trees # A component containing the tree objects related to that stratum.
#       res[[ i ]][[ "forest" ]] <- object@strat[[ i ]]@forest #A component containing the forest of trees related to that stratum.
#     }
#   # Tree Extraction for Non-Stratified Analysis (directly extracts tree-related information):
#   } else {
#     res[[ "trees" ]] <- object@trees # A component containing the tree objects.
#     res[[ "forest" ]] <- object@forest # A component containing the forest of trees.
#   }
# }

#-------------------------------------------------------------------------------
# method to make predictions for new data or to retrieve fitted values
#-------------------------------------------------------------------------------
# method with NULL retrieves fitted values from SurvRF object
#-------------------------------------------------------------------------------
# method returns a list object containing Func, mean, and? Prob
#-------------------------------------------------------------------------------
setMethod(f = ".Predict",
          signature = c(object = "SurvRF",
                        newdata = NULL),
          definition = function(object, newdata, ...) {
              # message("class_SurvRF.R: line 102")
              return( object@forest )
            })


#-------------------------------------------------------------------------------
# method to make predictions for new data or to retrieve fitted values
#-------------------------------------------------------------------------------
# method with data.frame calculates value for new data from SurvRF object
#-------------------------------------------------------------------------------
# method returns a list object containing Func, mean, and? Prob
#-------------------------------------------------------------------------------
#' @include predictSurvTree.R
#' #Internal function to predict value of tree
setMethod(f = ".Predict",
          signature = c(object = "SurvRF",
                        newdata = "data.frame"),
          definition = function(object,
                                newdata,
                                ...,
                                params) {
            # message("class_SurvRF.R: LINE 123: method .Predict for SurvRF")

              # verify that there are not new levels in the the data
              # this assumes that newdata has been passed in with
              # covariates in the order used in the analysis. This is
              # guaranteed if the data.frame is created from the model.
              xLevels <- lapply(X = newdata, FUN = levels)
              # print('test95')

              for (i in length(x = xLevels)) {
                if (is.null(x = xLevels[[ i ]]) &&
                    is.null(x = object@xLevels[[ i ]])) next
                if (any(! {xLevels[[ i ]] %in% object@xLevels[[ i ]]})) {
                  stop("new factor levels not present in the training data",
                       call. = FALSE)
                }
              }
              # print('test96')

              # verify that type of data is the same as the training data
              # type means numeric (nCat = 0), ordered factor (nCat = 1), or
              # factor (nCat = length(levels(x)))
              nCat <- sapply(X = xLevels, FUN = length)

              nCat <- ifelse(test = sapply(X = newdata, FUN = is.ordered),
                             yes = 1L,
                             no = nCat)

              if (any(unlist(x = object@nCat) != unlist(x = nCat))) {
                stop("type of predictors in newdata do not match the training data",
                       call. = FALSE)
              }
              # print('test97')
              nTree <- length(x = object@trees)

              # predict for first tree
              # .predictSurvTree() is an internal function defined in predictSurvTree.R

              newResult <- .predictSurvTree(x = newdata,
                                            params = params,
                                            nCat = nCat,
                                            nodes = object@trees[[ 1L ]])
              # print('test98')
              iTree <- 2L
              while (iTree <= nTree) {

                # predict for tree iTree; sum result
                # .predictSurvTree() is an internal function defined in predictSurvTree.R

                tmp <- .predictSurvTree(x = newdata,
                                        params = params,
                                        nCat = nCat,
                                        nodes = object@trees[[ iTree ]])
                # print('tmp')
                # print(tmp)

                newResult$Func <- newResult$Func + tmp$Func
                newResult$mean <- newResult$mean + tmp$mean
                if (!is.null(x = newResult$Prob)) {
                  newResult$Prob <- newResult$Prob + tmp$Prob
                }

                iTree <- iTree + 1L
              }
              # print("test99")

              newResult$Func <- newResult$Func / nTree
              newResult$mean <- newResult$mean / nTree
              if (!is.null(x = newResult$Prob)) {
                newResult$Prob <- newResult$Prob / nTree
              }
              # print("test100")

              return( newResult )

            })

#-------------------------------------------------------------------------------
# method to make predictions for new data or to retrieve fitted values
#-------------------------------------------------------------------------------
# method with data.frame calculates value for new data from SurvRF object
#-------------------------------------------------------------------------------
# method returns a list object containing Func, mean, and? Prob
#-------------------------------------------------------------------------------
setMethod(f = ".PredictAll",
          signature = c(object = "SurvRF"),
          definition = function(Phase, eps0, object, ..., newdata, model, txLevels, txName, params) {
            # message("method .PredictAll for SurvRF: LINE 211")
            # message("comes from class_ITRSurv.R line 355")

            # set all cases to receive first treatment
            newdata[[ txName ]] <- txLevels[1L]
            # print(2)
            # extract new model frame
            x <- stats::model.frame(formula = model, data = newdata)

            # remove response from x
            if (attr(x = terms(x = model), which = "response") == 1L) {
              x <- x[,-1L,drop = FALSE]
            }
            # print(3)
            # calculate the estimated values for this treatment level
            # .Predict() is a method; called here for objects of class SurvRF, which
            # is defined in this file
            res <- .Predict(object = object,
                            newdata = x,
                            params = params, ...)
            # CZ adding in AUS
            # message('class_SurvRF.R: Pooled Line 230: calculating AUS')
            Func_res = as.matrix(res$Func)
            new_tp0 = list()
            # id_test = 1
            # plot(Func_res[,id_test],
            #      main = sprintf("Plot Func_res for subject %s for Phase %s",
            #                     id_test,Phase))

            # Calculating time points so that we only keep the time point after last event
            # area_tp0 = list() # this is old as of August 2024
            for (subj in 1:ncol(Func_res)){
              # message("for subject: ", subj)
              probabilities0 = Func_res[,subj]
              # Find the index right after the last change in probability
              last_change_index0 <- 1
              for (i0 in 2:length(probabilities0)) {
                if (probabilities0[i0] != probabilities0[i0 - 1]) {
                  last_change_index0 <- i0
                }
              }
              # Time point (index) right after the last change in probability
              time_point_after_last_change0 <- last_change_index0 + 1
              new_tp0[[subj]] = params@timePoints[1:time_point_after_last_change0]
              # area_tp0[[subj]] = new_tp0 # this is old as of august 2024
            }

            res$Func <- list(res$Func)
            res$mean <- list(res$mean)
            res$Prob <- list(res$Prob)
            res$AUS_cut_tp = list(new_tp0)
            # res$AUS = list(area_trt0) # old as of August 2024

            # other treatments
            i <- 2L
            # repeat this process for each tx level
            while (i <= length(x = txLevels)) {
              # set all cases to receive the ith treatment
              newdata[[ txName ]][] <- txLevels[i]
              # extract new model frame
              x <- stats::model.frame(formula = model, data = newdata)
              # remove response from x
              if (attr(x = terms(x = model), which = "response") == 1L) {
                x <- x[,-1L,drop=FALSE]
              }
              # calculate the estimated values for this treatment level
              # .Predict() is a method; called here for objects of class SurvRF, which
              # is defined in file class_SurvRF.R
              tt <-  .Predict(object = object,
                              newdata = x,
                              params = params, ...)

              Func_tt = as.matrix(tt$Func)
              # area_tp_tt = list() # old as of August 2024
              new_tp_tt = list()
              # plot(Func_tt[,id_test],
              #      main = sprintf("Plot Func_tt for subject %s for Phase %s",
              #                     id_test,Phase))

              for (subj in 1:ncol(Func_tt)){
                probabilities_tt = Func_tt[,subj]
                # Find the index right after the last change in probability
                last_change_index_tt <- 1
                for (itt in 2:length(probabilities_tt)) {
                  if (probabilities_tt[itt] != probabilities_tt[itt - 1]) {
                    last_change_index_tt <- itt
                  }
                }
                # Time point (index) right after the last change in probability
                time_point_after_last_change_tt <- last_change_index_tt + 1
                # print(time_point_after_last_change_tt)
                new_tp_tt[[subj]] = params@timePoints[1:time_point_after_last_change_tt]
                # area_tp_tt[[subj]] = new_tp_tt # old as of August 2024

              }

              res[[ "Func" ]][[ i ]] <- tt$Func
              res[[ "mean" ]][[ i ]] <- tt$mean
              res[[ "Prob" ]][[ i ]] <- tt$Prob
              res[[ "AUS_cut_tp" ]][[i]] <- new_tp_tt # area_tp_tt is old as of August 2024

              i <- i + 1L
            }

            # message("%%%%%%%%%% AUS start")
            # Calculate AUS using all timepoints from params@timepoints

            # Initialize a list to store the areas under the curve for each treatment
            area_trt_list <- vector("list", length = length(x = txLevels))
            # Iterate over each treatment level
            treatment_index = 1L # start w/ treatment one then go through each treatment
            while (treatment_index <= length(x = txLevels)) {
              # message("treatment index: ", treatment_index)
              # Initialize variable to store the areas under the curve for the current treatment
              area_trt <- numeric(length = ncol(res[["Func"]][[treatment_index]]))

              # Iterate over each subject
              for (subject in seq_len(ncol(res[["Func"]][[treatment_index]]))) {
                tp_data = params@timePoints
                # Calculate the length of treatment data for the current subject and treatment
                length_tp <- length(tp_data)
                # Calculate the area under the curve for the current treatment level
                area_trt[subject] <- .calculate_area_under_curve(time_points = tp_data,
                                                                 survival_probabilities = res[["Func"]][[treatment_index]][, subject])
              }
              # Store the areas under the curve for the current treatment level
              area_trt_list[[treatment_index]] <- area_trt
              treatment_index <- treatment_index + 1L
            }

            # Assign the calculated areas under the curve to the result object
            i = 1L # start w/ treatment one then go through each treatment
            while (i <= length(x = txLevels)) {
              res[["AUS"]][[i]] <- area_trt_list[[i]]
              i = i + 1L
            }
            # message("%%%%%%%%%% AUS end")


            # message("%%%%%%%%%% start AUS_cut")
            # Initialize lists to store the maximum time points and corresponding treatment index for each subject
            max_tp_list <- vector("list", length = ncol(res[["Func"]][[1]]))
            max_tp_trtindex_list <- vector("list", length = ncol(res[["Func"]][[1]]))

            # Iterate over each subject
            for (subject in seq_len(ncol(res[["Func"]][[1]]))) {
              # Initialize max_tp and corresponding trt index for the current subject
              max_tp_subject <- 0
              max_tp_index_subject <- 0

              # Iterate over each treatment index
              for (treatment_index in seq_along(txLevels)) {
                # Get the length of treatment data for the current subject and treatment
                t_tp <- length(res[["AUS_cut_tp"]][[treatment_index]][[subject]])

                # Update max_tp_subject and max_tp_index_subject if necessary
                if (t_tp > max_tp_subject) {
                  max_tp_subject <- t_tp
                  max_tp_index_subject <- treatment_index
                }
              }
              # Store the maximum time point and corresponding index for the current subject
              max_tp_list[[subject]] <- max_tp_subject
              max_tp_trtindex_list[[subject]] <- max_tp_index_subject
            }

            # Initialize a list to store the areas under the curve for each treatment
            area_trt_list <- vector("list", length = length(x = txLevels))

            # Iterate over each treatment level
            treatment_index = 1L # start w/ treatment one then go through each treatment
            while (treatment_index <= length(x = txLevels)) {
              message("treatment index: ", treatment_index)
              # Initialize variable to store the areas under the curve for the current treatment
              area_trt <- numeric(length = ncol(res[["Func"]][[treatment_index]]))

              # Iterate over each subject
              for (subject in seq_len(ncol(res[["Func"]][[treatment_index]]))) {
                # Identify maximum time points per subject over all treatments
                tp_ind = max_tp_list[[subject]]
                tp_trtindex_ind = max_tp_trtindex_list[[subject]]
                # message("for subject ", subject, "the max time point is :", tp_ind, " and this is from treatment:", tp_trtindex_ind)
                tp_data = res[["AUS_cut_tp"]][[tp_trtindex_ind]][[subject]]
                # Calculate the length of treatment data for the current subject and treatment
                length_tp <- length(tp_data)
                # Calculate the area under the curve for the current treatment level
                area_trt[subject] <- .calculate_area_under_curve(time_points = tp_data,
                                                                 survival_probabilities = res[["Func"]][[treatment_index]][, subject])
              }
              # Store the areas under the curve for the current treatment level
              area_trt_list[[treatment_index]] <- area_trt
              treatment_index <- treatment_index + 1L
            }

            # Assign the calculated areas under the curve to the result object
            i = 1L # start w/ treatment one then go through each treatment
            while (i <= length(x = txLevels)) {
              res[["AUS_cut"]][[i]] <- area_trt_list[[i]]
              i = i + 1L
            }
            # message("%%%%%%%%%% end AUS_cut")

            ############################# below is old stuff for only 2 trts. current code is for 2 or more. #############################
            # area_trt0 = numeric()
            # area_trt1 = numeric()
            #
            # ## CURRENTLY THIS IS ONLY CODED FOR TWO TREATMENTS
            # for (subject in 1:ncol(res[["Func"]][[1]])){
            #   t0_tp = length(res[["AUS_cut_tp"]][[1]][[subject]])
            #   t1_tp = length(res[["AUS_cut_tp"]][[2]][[subject]])
            #   # WE WANT TO CALCULATE THE AREA UNDER THE CURVE FROM 0 TO THE MAXIMUM TIMEPOINT OF THE TWO
            #   tp_ind = max(t0_tp, t1_tp)
            #   if (tp_ind == t0_tp){
            #     tp = res[["AUS_cut_tp"]][[1]][[subject]]
            #   } else if (tp_ind == t1_tp){
            #     tp = res[["AUS_cut_tp"]][[2]][[subject]]
            #   } else{
            #     stop("class_SurvRF.R LINE 583")
            #   }
            #   print('calculating area')
            #   area_trt0[subject] =
            #     .calculate_area_under_curve(time_points = tp,
            #                                 survival_probabilities = res[["Func"]][[1]][,subject])
            #   #time_points = params@timePoints,
            #   # survival_probabilities = Func_res[,subj])
            #   area_trt1[subject] =
            #     .calculate_area_under_curve(time_points = tp,
            #                                 survival_probabilities = res[["Func"]][[2]][,subject])
            #   # area_trt1[subj] =
            #   #   .calculate_area_under_curve(time_points = params@timePoints,
            #   #                               survival_probabilities = Func_tt[,subj])
            # }
            #
            # res[["AUS"]][[1]] = area_trt0
            # res[["AUS"]][[2]] = area_trt1
            ##############################################################################################################################

            # message("TESITNGSLDFJSDKLF;J")

            # if (Phase == 1 | Phase == "Survival"){
            #   res_tmp1 <<- res
            # }
            #
            # # for now, just two treatments for MEAN ONLY - need to code for PROB
            # res$Ratio_Trt0 <- res[["AUS"]][[1]]/res[["AUS"]][[2]] #this is trt1/trt2. I want NONOPT/OPT
            # res$Ratio_Trt1 <- res[["AUS"]][[2]]/res[["AUS"]][[1]] #this is trt2/trt1. I want NONOPT/OPT.
            # res[["Ratio"]] = cbind(res$Ratio_Trt0, res$Ratio_Trt1)
            #
            # # if 0/X = 0 or X/0 = Inf then set ratio = difference
            # inf0indices = which(apply(res[["Ratio"]], 1, function(row) sum(row == 0) == 1))
            # res[["Ratio"]][inf0indices,] = cbind(res[["AUS"]][[1]][inf0indices]-res[["AUS"]][[2]][inf0indices],
            #                                      res[["AUS"]][[2]][inf0indices]-res[["AUS"]][[1]][inf0indices])
            #
            # if (Phase == 1 | Phase == "Survival"){
            #   res_tmp2 <<- res
            # }


            # so, we take optTx below in .optimal function and if optTx = 1 then we keep ratio;
            # but if optTx = 2 then we take the inverse of Ratio. Compare this to 1-eps0=1-0.1=0.9.
            #   mutate(Ratio = NonOpt_AUS/Opt_AUS) %>%
            #   #if ratio >= 0.9 then ratio_stop_ind = 1 (aka we stop b/c we know which treatment to pick; otherwise ind = 0 aka we keep going)
            #   mutate(Ratio_Stop_Ind = ifelse(Ratio >= (1-eps0), 1, 0)) #%>% filter(Ratio_Stop_Ind == 1)

            trtlevel <<- txLevels
            opt <- .optimal(Phase = Phase,
                            eps0 = eps0,
                            params = params,
                            predicted = res,
                            txLevels = txLevels)
            return( list("predicted" = res, "optimal" = opt) )
          })


# Class for storing survRF results for stratified analysis
#
# Class is not exported and is for internal convenience only
#
#  @slot strat A list object. The results of the tree building algorithm for
#   each tx group. List will be a list of SurvRF objects
#
# Methods
#   .Predict(object, newdata, ...) {defined}
#
setClass(Class = "SurvRFStratified",
         slots = c("strat" = "list"),
         contains = c("SurvRFObject"))

#-------------------------------------------------------------------------------
# method to make predictions for new data or to retrieve fitted values
#-------------------------------------------------------------------------------
# method with NULL retrieves fitted values from SurvRFStratified object
#-------------------------------------------------------------------------------
# method returns a list object each element a list containing Func, mean,
#   and? Prob
#-------------------------------------------------------------------------------
setMethod(f = ".Predict",
          signature = c(object = "SurvRFStratified",
                        newdata = NULL),
          definition = function(object, newdata, ...) {

            # print("method .Predict for SurvRFStratified for NULL newdata: LINE 470")
              res <- list()
              for (i in 1L:length(x = object@strat)) {
                # print(i)
                res[[ i ]] <- .Predict(object@strat[[ i ]], newdata = NULL, ...)
              }
              return( res )
            })

#-------------------------------------------------------------------------------
# method makes predictions for new data or retrieves fitted values
#-------------------------------------------------------------------------------
# method with data.frame calculated value for new data from SurvRF object
#-------------------------------------------------------------------------------
# method returns a list object each element a list containing Func, mean,
#   and? Prob
#-------------------------------------------------------------------------------
setMethod(f = ".Predict",
          signature = c(object = "SurvRFStratified",
                        newdata = "data.frame"),
          definition = function(object, newdata, ..., params) {

            # print("method .Predict for SurvRFStratified for data.frame newdata")

              # trt 0
              res <- list()
              res[[ 1L ]] <- .Predict(object = object@strat[[ 1L ]],
                                      newdata = newdata,
                                      params = params, ...)

              # trt 1
              i <- 2L
              while (i <= length(x = object@strat)) {
                res[[ i ]] <- .Predict(object = object@strat[[ i ]],
                                       newdata = newdata,
                                       params = params, ...)
                i <- i + 1L
              }

              return( res )

            })


#-------------------------------------------------------------------------------
# method makes predictions for new data or retrieves fitted values
#-------------------------------------------------------------------------------
# method with data.frame calculated value for new data from SurvRF object
#-------------------------------------------------------------------------------
# method returns a list object each element a list containing Func, mean,
#   and? Prob
#-------------------------------------------------------------------------------
setMethod(f = ".PredictAll",
          signature = c(object = "SurvRFStratified"),
          definition = function(Phase, eps0, epName = epName, object, ..., newdata, model, params, txLevels) {
            # message("class_SurvRF.R: LINE 392")
            # message(sprintf("method .PredictAll for SurvRFStratified for Phase: %s\n", Phase))

            if (Phase == "RE"){
              newdata_cov = newdata %>%
                     filter(!!sym(epName) == 0)
              # extract new model frame
              x <- stats::model.frame(formula = model, data = newdata_cov)
            } else{
              # extract new model frame
              x <- stats::model.frame(formula = model, data = newdata)
            }
            # remove response from x
            if (attr(x = terms(x = model), which = "response") == 1L) {
              x <- x[,-1L,drop=FALSE]
            }

            # calculate the estimated values for this treatment level
            # .Predict() is a method; called here for objects of class SurvRFStratified, which
            # is defined in this file

            # object <<- object
            # x <<- x
            # params <<- params

            # define object and params so it's phase-specific
            # predicting for first trt
            # message("Below uses .Predict from class_SurvRF.R: .Predict from LINE 14")
            res <- .Predict(object = object@strat[[ 1L ]], # trt0
                            newdata = x,
                            params = params,
                            ...)

            testing_paramstp <<- params@timePoints
            # message("end of .Predict from class_SurvRF.R: .Predict from LINE 14")
            res0 <<- res

            # below is updated for this section for stratified analysis
            # need to update code following this for pooled analysis (have not done yet as of Nov 2024 for RE.)
            if (Phase == "RE"){
              res$Func <- list(res$Func)
              res$mean <- list(res$mean)
              res$Prob <- list(res$Prob)
              i <- 2L
              while (i <= length(x = object@strat)) {
                tt_RE <<- .Predict(object = object@strat[[ i ]],
                               newdata = x,
                               params = params, ...)

                res[[ "Func" ]][[ i ]] <- tt_RE$Func
                res[[ "mean" ]][[ i ]] <- tt_RE$mean #cbind(res$mean, tt_RE$mean)
                res[[ "Prob" ]][[ i ]] <- tt_RE$Prob #cbind(res$Prob, tt_RE$Prob)

                i <- i + 1L
              }
            }

            if (Phase == "Survival" | Phase == 1 | Phase == "CR"){
              # CZ adding in AUS
              # message('class_SurvRF.R: Line 419: calculating AUS')
              Func_res = as.matrix(res$Func)
              new_tp0 = list()

              # below. we obtain the time points to calculate AUS_cut which is for AUS under timepoints from 0 to tiempoint after last event for each subject (varies by subject)
              for (subj in 1:ncol(Func_res)){
                # message("for subject: ", subj)
                probabilities0 = Func_res[,subj]
                # Find the index right after the last change in probability
                last_change_index0 <- 1
                for (i0 in 2:length(probabilities0)) {
                  if (probabilities0[i0] != probabilities0[i0 - 1]) {
                    last_change_index0 <- i0
                  }
                }
                # Time point (index) right after the last change in probability
                time_point_after_last_change0 <- last_change_index0 + 1
                new_tp0[[subj]] = params@timePoints[1:time_point_after_last_change0] #testing_paramstp
              }

              res$Func <- list(res$Func)
              res$mean <- list(res$mean)
              res$Prob <- list(res$Prob)
              res$AUS_cut_tp <- list(new_tp0)

              # other treatments
              i <- 2L
              while (i <= length(x = object@strat)) {
                # message(i)
                tt <- .Predict(object = object@strat[[ i ]],
                               newdata = x,
                               params = params, ...)
                Func_tt = as.matrix(tt$Func)

                new_tp_tt = list()
                for (subj in 1:ncol(Func_tt)){
                  # message("for subject: ", subj)

                  probabilities_tt = Func_tt[,subj]
                  # Find the index right after the last change in probability
                  last_change_index_tt <- 1
                  for (itt in 2:length(probabilities_tt)) {
                    if (probabilities_tt[itt] != probabilities_tt[itt - 1]) {
                      last_change_index_tt <- itt
                    }
                  }
                  # Time point (index) right after the last change in probability
                  time_point_after_last_change_tt <- last_change_index_tt + 1
                  # print(time_point_after_last_change_tt)
                  new_tp_tt[[subj]] = params@timePoints[1:time_point_after_last_change_tt] #testing_paramstp
                }

                res[[ "Func" ]][[ i ]] <- tt$Func
                res[[ "mean" ]][[ i ]] <- tt$mean
                res[[ "Prob" ]][[ i ]] <- tt$Prob
                res[[ "AUS_cut_tp"]][[i]] <- new_tp_tt

                i <- i + 1L
              }
              stratobject <<- object@strat
              current_res <<- res

              # message("%%%%%%%%%%")
              #######################################################################################################################################################################
              #######################################################################################################################################################################
              #######################################################################################################################################################################
              #######################################################################################################################################################################
              # OLD: only for 2 trts - current code is generalized to more than 2 treatments
              # area_trt0 = numeric()
              # area_trt1 = numeric()
              #
              # print("CURRENTLY THIS IS ONLY CODED FOR TWO TREATMENTS")
              # ## CURRENTLY THIS IS ONLY CODED FOR TWO TREATMENTS
              # for (subject in 1:ncol(res[["Func"]][[1]])){
              #   t0_tp = length(res[["AUS_cut_tp"]][[1]][[subject]])
              #   t1_tp = length(res[["AUS_cut_tp"]][[2]][[subject]])
              #   # WE WANT TO CALCULATE THE AREA UNDER THE CURVE FROM 0 TO THE MAXIMUM TIMEPOINT OF THE TWO
              #   tp_ind = max(t0_tp, t1_tp)
              #   if (tp_ind == t0_tp){
              #     tp = res[["AUS_cut_tp"]][[1]][[subject]]
              #   } else if (tp_ind == t1_tp){
              #     tp = res[["AUS_cut_tp"]][[2]][[subject]]
              #   } else{
              #     stop("class_SurvRF.R LINE 583")
              #   }
              #   area_trt0[subject] =
              #     .calculate_area_under_curve(time_points = tp,
              #                                 survival_probabilities = res[["Func"]][[1]][,subject])
              #   #time_points = params@timePoints,
              #   # survival_probabilities = Func_res[,subj])
              #   area_trt1[subject] =
              #     .calculate_area_under_curve(time_points = tp,
              #                                 survival_probabilities = res[["Func"]][[2]][,subject])
              #   # area_trt1[subj] =
              #   #   .calculate_area_under_curve(time_points = params@timePoints,
              #   #                               survival_probabilities = Func_tt[,subj])
              # }
              #
              # res[["AUS"]][[1]] = area_trt0
            }

            #######################################################################################################################################################################
            #######################################################################################################################################################################
            #######################################################################################################################################################################
            #######################################################################################################################################################################

            if (Phase == "Survival" | Phase == "CR"){
              # message("%%%%%%%%%% AUS start")
              # Calculate AUS using all timepoints from params@timepoints

              # Initialize a list to store the areas under the curve for each treatment
              area_trt_list <- vector("list", length = length(x = stratobject))
              # Iterate over each treatment level
              treatment_index = 1L # start w/ treatment one then go through each treatment
              while (treatment_index <= length(x = stratobject)) {
                # message("treatment index: ", treatment_index)
                # Initialize variable to store the areas under the curve for the current treatment
                area_trt <- numeric(length = ncol(res[["Func"]][[treatment_index]]))

                # Iterate over each subject
                for (subject in seq_len(ncol(res[["Func"]][[treatment_index]]))) {
                  tp_data = params@timePoints
                  # Calculate the length of treatment data for the current subject and treatment
                  length_tp <- length(tp_data)
                  # Calculate the area under the curve for the current treatment level
                  area_trt[subject] <- .calculate_area_under_curve(time_points = tp_data,
                                                                   survival_probabilities = res[["Func"]][[treatment_index]][, subject])
                }
                # Store the areas under the curve for the current treatment level
                area_trt_list[[treatment_index]] <- area_trt
                treatment_index <- treatment_index + 1L
              }

              # Assign the calculated areas under the curve to the result object
              i = 1L # start w/ treatment one then go through each treatment
              while (i <= length(x = stratobject)) {
                res[["AUS"]][[i]] <- area_trt_list[[i]]
                i = i + 1L
              }
              # message("%%%%%%%%%% AUS end")

              # message("%%%%%%%%%% AUS_cut start")
              # NOW WE CALCULATE AUS FOR THE CUT OFF TPS , AKA VARIES PER PERSON, USUALLY FOR WHEN TIMEPOINTS OR NTIMES WAS ENTERED - NOT WHEN TIMEPOITNSSURVIVAL/ENDPOINT WERE ENTERED
              # WE DO THIS TO AVOID WHEN THERE ARE LONG TAILS IN THE CURVES
              # Initialize lists to store the maximum time points and corresponding treatment index for each subject
              max_tp_list <- vector("list", length = ncol(res[["Func"]][[1]]))
              max_tp_trtindex_list <- vector("list", length = ncol(res[["Func"]][[1]]))
              # Iterate over each subject
              for (subject in seq_len(ncol(res[["Func"]][[1]]))) {
                # Initialize max_tp and corresponding trt index for the current subject
                max_tp_subject <- 0
                max_tp_index_subject <- 0

                # Iterate over each treatment index
                for (treatment_index in seq_along(stratobject)) {
                  # Get the length of treatment data for the current subject and treatment
                  t_tp <- length(res[["AUS_cut_tp"]][[treatment_index]][[subject]])

                  # Update max_tp_subject and max_tp_index_subject if necessary
                  if (t_tp > max_tp_subject) {
                    max_tp_subject <- t_tp
                    max_tp_index_subject <- treatment_index
                  }
                }
                # Store the maximum time point and corresponding index for the current subject
                max_tp_list[[subject]] <- max_tp_subject
                max_tp_trtindex_list[[subject]] <- max_tp_index_subject
              }
              # Initialize a list to store the areas under the curve for each treatment
              area_trt_list <- vector("list", length = length(x = stratobject))
              # Iterate over each treatment level
              treatment_index = 1L # start w/ treatment one then go through each treatment
              while (treatment_index <= length(x = stratobject)) {
                # message("treatment index: ", treatment_index)
                # Initialize variable to store the areas under the curve for the current treatment
                area_trt <- numeric(length = ncol(res[["Func"]][[treatment_index]]))

                # Iterate over each subject
                for (subject in seq_len(ncol(res[["Func"]][[treatment_index]]))) {
                  # Identify maximum time points per subject over all treatments
                  tp_ind = max_tp_list[[subject]]
                  tp_trtindex_ind = max_tp_trtindex_list[[subject]]
                  # message("for subject ", subject, "the max time point is :", tp_ind, " and this is from treatment:", tp_trtindex_ind)
                  tp_data = res[["AUS_cut_tp"]][[tp_trtindex_ind]][[subject]]
                  # Calculate the length of treatment data for the current subject and treatment
                  length_tp <- length(tp_data)
                  # Calculate the area under the curve for the current treatment level
                  area_trt[subject] <- .calculate_area_under_curve(time_points = tp_data,
                                                                   survival_probabilities = res[["Func"]][[treatment_index]][, subject])
                }
                # Store the areas under the curve for the current treatment level
                area_trt_list[[treatment_index]] <- area_trt
                treatment_index <- treatment_index + 1L
              }

              # Assign the calculated areas under the curve to the result object
              i = 1L # start w/ treatment one then go through each treatment
              while (i <= length(x = stratobject)) {
                res[["AUS_cut"]][[i]] <- area_trt_list[[i]]
                i = i + 1L
              }
              # message("%%%%%%%%%% AUS_cut end")
            } # end of AUS for Phase 1 or Phase 2CR

            # below is old stuff. I think I can delete.
            # if (Phase == 1 | Phase == "Survival"){
            #   res_tmp1 <<- res
            # }
            #
            # # for now, just two treatments for MEAN ONLY - need to code for PROB
            # res$Ratio_Trt0 <- res[["AUS"]][[1]]/res[["AUS"]][[2]] #this is trt1/trt2. I want NONOPT/OPT
            # res$Ratio_Trt1 <- res[["AUS"]][[2]]/res[["AUS"]][[1]] #this is trt2/trt1. I want NONOPT/OPT.
            # res[["Ratio"]] = cbind(res$Ratio_Trt0, res$Ratio_Trt1)
            #
            # # if 0/X = 0 or X/0 = Inf then set ratio = difference
            # inf0indices = which(apply(res[["Ratio"]], 1, function(row) sum(row == 0) == 1))
            # res[["Ratio"]][inf0indices,] = cbind(res[["AUS"]][[1]][inf0indices]-res[["AUS"]][[2]][inf0indices],
            #                                      res[["AUS"]][[2]][inf0indices]-res[["AUS"]][[1]][inf0indices])
            #
            # if (Phase == 1 | Phase == "Survival"){
            #   res_tmp2 <<- res
            # }

            # so, we take optTx below in .optimal function and if optTx = 1 then we keep ratio;
            # but if optTx = 2 then we take the inverse of Ratio. Compare this to 1-eps0=1-0.1=0.9.
            #   mutate(Ratio = NonOpt_AUS/Opt_AUS) %>%
            #   #if ratio >= 0.9 then ratio_stop_ind = 1 (aka we stop b/c we know which treatment to pick; otherwise ind = 0 aka we keep going)
            #   mutate(Ratio_Stop_Ind = ifelse(Ratio >= (1-eps0), 1, 0)) #%>% filter(Ratio_Stop_Ind == 1)
            # if (Phase == "Survival"){
            #   View(res)
            # }

            # print("WE ARE NOW RUNNING .OPTIMAL FROM STRATIFIED CLASS_SURVRF.R LINE 760")
            opt <- .optimal(Phase = Phase,
                            eps0 = eps0,
                            params = params,
                            predicted = res,
                            txLevels = txLevels)
            # print('end of .optimal')
            if (Phase == "Survival" | Phase == 1){
              opt_p1 <<- opt
            }
            if (Phase == "CR" | Phase == "RE" | Phase == 2){
              opt_p2 <<- opt
            }
            result_list =  list("predicted" = res,
                                "optimal" = opt)
            return(result_list)
          })

###############################################################################################
###############################################################################################
###################################### .optimal function ######################################
###############################################################################################
###############################################################################################
# Phase = 1; eps0 = mean_tol1; params = pa; predicted = re; txLevels = trtlevel
.optimal <- function(Phase, eps0, params, predicted, txLevels) {

  # print("%%%%% beginning .optimal function in class_SurvRF.R %%%%%%%")
  # message("Phase:", Phase)
  # if (Phase == "RE"){
    # View(predicted)
  # }
  # print(length(optTx))
  # if (Phase == "RE"){stop("stpoing")}
  # crit can only be mean, prob, area, or mean.prob.combo
  # 'mean', 'area', 'prob', 'mean.prob.comb'
  crit <- .CriticalValueCriterion(params)
  # message(".optimal function with crit: ", crit, " for Phase ", Phase)

  # initialize empty matrices for mean_trts and area_trts based on trt1
  mean_trts = area_trts = prob_trts = matrix(nrow = ncol(predicted[["Func"]][[1]]), ncol = length(txLevels))

  # Get the function based on the Phase
  if (Phase == "Survival" | Phase == 1){
    findMax1 = TRUE # aka find max
    opp1 = FALSE # aka find min
    predd_surv <<- predicted
  } else if (Phase == 2 | Phase == "CR" | Phase == "RE"){ # for now Phase = 2 means CR. but need to change/edit later when we add RE..
    findMax1 = FALSE # aka find min
    opp1 = TRUE # aka find max
    predd_ep <<- predicted
  } else{
    stop("Phase is not coded yet:", Phase)
  }

  # View(predicted)
  for (trt in seq_along(txLevels)){
    mean_trts[, trt] <- unlist(predicted$mean[[trt]])
    if (Phase != "RE"){
      area_trts[, trt] <- unlist(predicted$AUS[[trt]])
    }
    if (crit == "mean.prob.combo" | crit == "prob"){
      ttmp <<-  unlist(predicted$Prob[[trt]])
      prob_trts[, trt] <- unlist(predicted$Prob[[trt]])
    }
  }
  # print("test1")
  # if (Phase == 1 | Phase == "Survival"){
  #   pp1_mean <<- mean_trts
  #   pp1_area <<- area_trts
  #   if (crit == "mean.prob.combo"){
  #     pp1_prob <<- prob_trts
  #   }
  # }
  # if (Phase == 2 | Phase == "CR"){
  #   pp2_mean <<- mean_trts
  #   pp2_area <<- area_trts
  #   if (crit == "mean.prob.combo"){
  #     pp2_prob <<- prob_trts
  #   }
  # }

  # starting based on what crit is
  if (crit == "mean" | crit == "area" | crit == "prob") {

    eps0_ratio = eps0[1] # we care about this one; later delete difference and edit in verifytol etc
    eps0_diff = eps0[2]

    if (crit == "mean"){
      # print("test2")
      dat_trts = mean_trts
      pred_trts = predicted[["mean"]]
    } else if (crit == "area"){
      # print("test2.2")
      dat_trts = area_trts
      pred_trts = predicted[["AUS"]]
    } else{
      # print("test2.3")
      dat_trts = prob_trts
      pred_trts = predicted[["Prob"]]
    }
    # print("test3")
    # View(dat_trts)
    # View(pred_trts)
    # identify which Trt contains the maximum expected survival time or min expected CI time
    optTx <- apply(X = dat_trts,
                   MARGIN = 1L, # across rows (each person)
                   FUN = .whichExtremum,
                   tieMethod = params@tieMethod,
                   findMax = findMax1)
    # print("test4")
    # extract the survival or cumulative incidence mean/probability at optimal tx
    V_star <- matrix(data = 0.0,
                     nrow = length(x = optTx),
                     ncol = 1)
    for (i in 1L:length(optTx)) { # i = subject
      # message("subject:", i)
      V_star[i] <- pred_trts[[optTx[i]]][i]
    }
    # print("test5")
    tol_V_star = (1-eps0_ratio)*V_star

    if (!(Phase %in% c(1, "Survival", 2, "CR", "RE"))) {
      stop("Phase must be either 1, 'Survival', 2, or 'CR', or 'RE'.")
    }
    if (Phase == 1 | Phase == "Survival"){
      # for each subject, compare all trt survival prob/means to V_star (aka argmax) and see if its within eps0_ratio

      V_mat = matrix(nrow = length(optTx), ncol = length(seq_along(pred_trts)))

      for (subj in 1L:length(optTx)){
        # message("subj:", subj)
        tol_V_star_subj = tol_V_star[subj]

        for (trt in seq_along(pred_trts)){
          # cat("\ntrt:",trt,"\n")
          V_trt = pred_trts[[trt]][subj]

          V_mat[subj,trt] = ifelse(V_trt >= tol_V_star_subj, 1, 0)
        } # end trts
      } # end subject
      # print("test6")
      # Calculate the row sums to get number of chosen trts
      row_sums <- rowSums(V_mat); all(row_sums > 0) # should always be true b/c one trt will always equal the opt b/c of how derived
      # Add the row sums as a new column to the matrix
      V_df <- cbind(V_mat, NumTrts = row_sums) %>% as.data.frame()
      # Create Stop_Ind column based on NumTrts. Stop_Ind = 1 if we stop at Phase 1; = 0 if we continue to Phase 2.
      V_df1 <- V_df %>%
        mutate(Stop_Ind = ifelse(NumTrts > 1, 0, 1))
      num_trts = V_df1$NumTrts
      stop_ind = V_df1$Stop_Ind
      predd_surv[["Stopping_Ind"]] <<- stop_ind
      # print(stop_ind)

    } # end of Phase = 1
  } # end of crit = mean | crit = area

  else if (crit == "mean.prob.combo") {
    message("mean.prob.combo")
    # View(area_trts)

    # identify which element contains
    # the maximum expected survival time
    # or min expected CI time
    optTx <- apply(X = area_trts,
                   MARGIN = 1L,
                   FUN = .whichExtremum,
                   tieMethod = params@tieMethod,
                   findMax = findMax1)


    # index for which the survival probability is max
    # or which index for CI prob is min
    tmp <- apply(X = prob_trts,
                 MARGIN = 1L,
                 FUN = .whichExtremum,
                 tieMethod = "NA",
                 findMax = findMax1)

    # for those that are not tied, replace mean survival/cif time with mean survival/cif probability
    isna <<- is.na(x = tmp)
    mean_index = which(isna)
    prob_index = which(!isna)
    optTx[!isna] <- tmp[!isna]

    # mean.prob.combo is actually area.prob.combo (aka we use the truncated area for determining if we continue to P2 or not but use mean for any actual calculations.)
    V_mat = tol_V_star_list = Num_Trts = Stop_Ind = list()
    V_star_AUS = V_star_Prob = numeric()
    V_mat[["AUS"]] <- matrix(NA, nrow = length(optTx), ncol = length(predicted[["AUS"]]))
    V_mat[["Prob"]] <- matrix(NA, nrow = length(optTx), ncol = length(predicted[["Prob"]]))
    eps0_ratio_mean = eps0[1]
    eps0_ratio_prob = eps0[2]
    eps0_diff_mean = eps0[3]
    eps0_diff_prob = eps0[4]

    if (!(Phase %in% c(1, "Survival", 2, "CR", "RE"))) {
      stop("Phase must be either 1, 'Survival', 2, or 'CR' or 'RE'.")
    }

    if (Phase == 1 | Phase == "Survival"){

      # extract the survival or cumulative incidence mean/probability at optimal tx
      V_star <- matrix(data = 0.0,
                       nrow = length(x = optTx),
                       ncol = 1)

      for (i in 1L:length(optTx)) { # i = subject
        # message("subject:", i)
        V_star_AUS[i] <-predicted[["AUS"]][[optTx[i]]][i]
        V_star_Prob[i] <-predicted[["Prob"]][[optTx[i]]][i]
      }

      tol_V_star_list[["AUS"]] = (1-eps0_ratio_mean)*V_star_AUS
      tol_V_star_list[["Prob"]] = (1-eps0_ratio_prob)*V_star_Prob

      for (scenario in c("AUS", "Prob")) {
        # print(scenario)
        pred_vals <- predicted[[scenario]]
        tol_V_star <- tol_V_star_list[[scenario]]

        for (subj in seq_along(optTx)) {
          # obtain (1-alpha)V* for subject for current scenario
          tol_V_star_subj <- tol_V_star[subj]

          for (trt in seq_along(pred_vals)) {
            V_trt <- pred_vals[[trt]][subj]
            # indicator to continue to P2 ( = 1 ) ==> this is a potential optimal trt candidate
            V_mat[[scenario]][subj,trt] <- ifelse(V_trt >= tol_V_star_subj, 1, 0)
          } # end trt
        } # end subj

        # Calculate the row sums to get number of chosen trts
        row_sums <- rowSums(V_mat[[scenario]]); if (!all(row_sums > 0)){stop("ERROR: not all rows_sums > 0")} # should always be true b/c one trt will always equal the opt b/c of how derived
        # Add the row sums as a new column to the matrix
        V_df <- cbind(V_mat[[scenario]], NumTrts = row_sums) %>% as.data.frame()
        # Create Stop_Ind column based on NumTrts. Stop_Ind = 1 if we stop at Phase 1; = 0 if we continue to Phase 2.
        V_df1 <- V_df %>%
          mutate(Stop_Ind = ifelse(NumTrts > 1, 0, 1))
        Num_Trts[[scenario]] = V_df1$NumTrts
        Stop_Ind[[scenario]] = V_df1$Stop_Ind

      } # end scenario

      Num_Trts0 = Num_Trts[["AUS"]]
      Num_Trts0[!isna] <- Num_Trts[["Prob"]][!isna]
      Stop_Ind0 = Stop_Ind[["AUS"]]
      Stop_Ind0[!isna] <- Stop_Ind[["Prob"]][!isna]
      num_trts = Num_Trts0
      stop_ind = Stop_Ind0
      predd_surv[["Stopping_Ind"]] <<- stop_ind
    } # end of Phase = 1
  } # end of mean.prob.combo
  else{
    stop("CRIT NOT DEFINED")
  }

  # print("test7")
  if (Phase == 2 | Phase == "CR" | Phase == "RE"){
    # message('phase 2')
    # dont care about going to another phase
    num_trts = rep(99, length(optTx))
    stop_ind = rep(99, length(optTx))
    predd_ep[["Stopping_Ind"]] <<- stop_ind
    } # end of CIF
  # print("test8")
  # extract the survival or cumulative incidence function at optimal tx
  optSv <- matrix(data = 0.0,
                  nrow = length(x = optTx),
                  ncol = .NTimes(params))

  for (i in 1L:length(optTx)) {
    optSv[i,] <- predicted$Func[[optTx[i]]][,i]
  }

  # print("test9")
  # BELOW IS OLD STUFF.
  # # initialize empty matrices for mean_trts and area_trts based on trt1
  # mean_trts = area_trts = matrix(nrow = ncol(predicted[["Func"]][[1]]), ncol = length(txLevels))
  #
  # # Get the function based on the Phase
  # if (Phase == "Survival" | Phase == 1){
  #   findMax1 = TRUE # aka find max
  #   opp1 = FALSE # aka find min
  #   predd_surv <<- predicted
  # } else if (Phase == 2 | Phase == "CR"){ # for now Phase = 2 means CR. but need to change/edit later when we add RE..
  #   findMax1 = FALSE # aka find min
  #   opp1 = TRUE # aka find max
  #   predd_ep <<- predicted
  # } else{
  #   stop("Phase is not coded yet:", Phase)
  # }
  #
  # # 'mean', 'prob', 'mean.prob.comb'
  # crit <- .CriticalValueCriterion(params)
  # message("\ncrit: ", crit)
  #
  # for (trt in seq_along(txLevels)){
  #   print(trt)
  #   mean_trts[, trt] <- unlist(predicted$mean[[trt]])
  #   area_trts[, trt] <- unlist(predicted$AUS[[trt]])
  # }
  # # mean_trts <<- cbind(predicted$mean[[1]],predicted$mean[[2]])
  # # area_trts <<- cbind(predicted$AUS[[1]],predicted$AUS[[2]])
  #
  # if (crit == "mean") {
  #   # identify which element contains
  #   # the maximum expected survival time
  #   # or min expected CI time
  #   if (Phase == 1 | Phase == "Survival"){
  #     pp1_mean <<- mean_trts
  #     pp1_area <<- area_trts
  #   }
  #   if (Phase == 2 | Phase == "CR"){
  #     pp2 <<- mean_trts
  #     pp2_area <<- area_trts
  #   }
  #   optTx <- apply(X = mean_trts,
  #                  MARGIN = 1L, # across rows (each person)
  #                  FUN = .whichExtremum,
  #                  tieMethod = params@tieMethod,
  #                  findMax = findMax1)
  #   nonoptTx <- apply(X = mean_trts,
  #                     MARGIN = 1L,
  #                     FUN = .whichExtremum,
  #                     tieMethod = params@tieMethod,
  #                     findMax = opp1)
  #
  #   if (Phase == 1 | Phase == "Survival"){
  #     # for now, just two treatments for MEAN ONLY - need to code for PROB
  #     predicted$Ratio_Trt0 <- predicted[["mean"]][[1]]/predicted[["mean"]][[2]] #this is trt1/trt2. I want NONOPT/OPT
  #     predicted$Ratio_Trt1 <- predicted[["mean"]][[2]]/predicted[["mean"]][[1]] #this is trt2/trt1. I want NONOPT/OPT.
  #     predicted[["Ratio"]] = cbind(predicted$Ratio_Trt0, predicted$Ratio_Trt1)
  #
  #     # if 0/X = 0 or X/0 = Inf then set ratio = difference
  #     inf0indices = which(apply(predicted[["Ratio"]], 1, function(row) sum(row == 0) == 1))
  #     predicted[["Ratio"]][inf0indices,] = cbind(predicted[["mean"]][[1]][inf0indices]-predicted[["mean"]][[2]][inf0indices],
  #                                                predicted[["mean"]][[2]][inf0indices]-predicted[["mean"]][[1]][inf0indices])
  #
  #     }
  #
  #   # print('test6')
  # } else if (crit == "mean.prob.combo") {
  #   # print(eps0)
  #   # if (Phase == 2 | Phase == "CR"){
  #   #   print(Phase)
  #   #   print(crit)
  #   #   print(opp1)
  #   # }
  #   prob_trts <<- cbind(predicted$Prob[[1]],predicted$Prob[[2]])
  #   # identify which element contains
  #   # the maximum expected survival time
  #   # or min expected CI time
  #   optTx <- apply(X = mean_trts,
  #                  MARGIN = 1L,
  #                  FUN = .whichExtremum,
  #                  tieMethod = params@tieMethod,
  #                  findMax = findMax1)
  #
  #   Ratio_Trt0_mean <- predicted[["mean"]][[1]]/predicted[["mean"]][[2]] #this is trt1/trt2. I want NONOPT/OPT
  #   Ratio_Trt1_mean <- predicted[["mean"]][[2]]/predicted[["mean"]][[1]] #this is trt2/trt1. I want NONOPT/OPT.
  #   Ratio_mean = cbind(Ratio_Trt0_mean, Ratio_Trt1_mean)
  #
  #   # index for which the survival probability is max
  #   # or which index for CI prob is min
  #   tmp <- apply(X = prob_trts,
  #                MARGIN = 1L,
  #                FUN = .whichExtremum,
  #                tieMethod = "NA",
  #                findMax = findMax1)
  #
  #   Ratio_Trt0_prob <- predicted[["Prob"]][[1]]/predicted[["Prob"]][[2]] #this is trt1/trt2. I want NONOPT/OPT
  #   Ratio_Trt1_prob <- predicted[["Prob"]][[2]]/predicted[["Prob"]][[1]] #this is trt2/trt1. I want NONOPT/OPT.
  #   Ratio_prob = cbind(Ratio_Trt0_prob, Ratio_Trt1_prob)
  #
  #   # for those that are not tied, replace mean survival/cif time with mean survival/cif probability
  #   isna <<- is.na(x = tmp)
  #   mean_index = which(isna)
  #   prob_index = which(!isna)
  #   optTx[!isna] <- tmp[!isna]
  #   ratio = Ratio_mean
  #   ratio[!isna] <- Ratio_prob[!isna]
  #   predicted[["Ratio"]] = ratio
  #
  #   nonoptTx <- apply(X = mean_trts,
  #                  MARGIN = 1L,
  #                  FUN = .whichExtremum,
  #                  tieMethod = params@tieMethod,
  #                  findMax = opp1)
  #   # index for which the survival probability is max
  #   # or which index for CI prob is min
  #   tmp_non <- apply(X = prob_trts,
  #                MARGIN = 1L,
  #                FUN = .whichExtremum,
  #                tieMethod = "NA",
  #                findMax = opp1)
  #   # for those that are not tied, replace mean survival/cif time
  #   # with mean survival/cif probability
  #   isna_non <- is.na(x = tmp_non)
  #   nonoptTx[!isna_non] <- tmp_non[!isna_non]
  #
  #
  #   # if 0/X = 0 or X/0 = Inf then set ratio = difference
  #   inf0indices <<- which(apply(predicted[["Ratio"]], 1, function(row) sum(row == 0) == 1))
  #   predicted[["Ratio"]][inf0indices,] = ifelse(inf0indices %in% which(!isna),
  #                                               cbind(
  #                                                 predicted[["Prob"]][[1]][inf0indices]-predicted[["Prob"]][[2]][inf0indices],
  #                                                 predicted[["Prob"]][[2]][inf0indices]-predicted[["Prob"]][[1]][inf0indices]
  #                                                 ),
  #                                               cbind(
  #                                                 predicted[["mean"]][[1]][inf0indices]-predicted[["mean"]][[2]][inf0indices],
  #                                                 predicted[["mean"]][[2]][inf0indices]-predicted[["mean"]][[1]][inf0indices]
  #                                               )
  #   )
  #   # print('test7')
  # } else if (crit == "area") {
  #   # identify which element contains the maximum expected survival probability
  #   optTx <- apply(X = area_trts,
  #                  MARGIN = 1L,
  #                  FUN = .whichExtremum,
  #                  tieMethod = params@tieMethod,
  #                  findMax = findMax1)
  #   nonoptTx <- apply(X = area_trts,
  #                  MARGIN = 1L,
  #                  FUN = .whichExtremum,
  #                  tieMethod = params@tieMethod,
  #                  findMax = opp1)
  #   # print('test9')
  #
  #   if (Phase == 1 | Phase == "Survival"){
  #     # for now, just two treatments for MEAN ONLY - need to code for PROB
  #     predicted$Ratio_Trt0 <- predicted[["AUS"]][[1]]/predicted[["AUS"]][[2]] #this is trt1/trt2. I want NONOPT/OPT
  #     predicted$Ratio_Trt1 <- predicted[["AUS"]][[2]]/predicted[["AUS"]][[1]] #this is trt2/trt1. I want NONOPT/OPT.
  #     predicted[["Ratio"]] = cbind(predicted$Ratio_Trt0, predicted$Ratio_Trt1)
  #
  #     # if 0/X = 0 or X/0 = Inf then set ratio = difference
  #     inf0indices = which(apply(predicted[["Ratio"]], 1, function(row) sum(row == 0) == 1))
  #     predicted[["Ratio"]][inf0indices,] = cbind(predicted[["AUS"]][[1]][inf0indices]-predicted[["AUS"]][[2]][inf0indices],
  #                                                predicted[["AUS"]][[2]][inf0indices]-predicted[["AUS"]][[1]][inf0indices])
  #   }
  #
  # } else{
  #   stop("CRIT NOT DEFINED")
  # }
  # # print('test10')
  # # extract the survival or cumulative indcidence function at optimal tx
  # optSv <- matrix(data = 0.0,
  #                 nrow = length(x = optTx),
  #                 ncol = .NTimes(params))
  #
  # for (i in 1L:length(optTx)) {
  #   # print(i)
  #   optSv[i,] <- predicted$Func[[optTx[i]]][,i]
  # }
  # # View(optSv)
  #
  # if (Phase == 1 | Phase == "Survival"){
  #   # extract the ratio at nonopt tx to select non-opt/opt ratio
  #   your_matrix0 <<- predicted$Ratio
  #   #always do non-opt over opt aka ratio is less than 1.
  #
  #   # if 0/X = 0 or X/0 = Inf then set ratio = difference
  #   inf0indices = which(apply(predicted[["Ratio"]], 1, function(row) any(row == 0 | row == Inf)))
  #   # this should be 0
  #   # print(inf0indices)
  #
  #   # edit predicted$Ratio so we avoid NaN and Inf/0
  #   for (row1 in 1:nrow(predicted$Ratio)){
  #     for (col1 in 1:2){
  #       # if 0/0 = Nan then set ratio =1
  #       if (is.nan(predicted$Ratio[row1,col1])){
  #         predicted$Ratio[row1,col1] = 1
  #       }
  #     }
  #   }
  #
  #   opt_tx_tmp <<- optTx
  #   p <<- predicted
  #   your_matrix1 <<- predicted$Ratio
  #
  #   min_ratio <<- apply(predicted$Ratio,
  #                       1, #across rows (each person)
  #                       min) # we know denom is always bigger or equal to num
  #   if (any(min_ratio>1)){stop("class_SurvRF.R LINE 732: min ratio cant be more than 1")}
  #
  #   if (crit == "mean.prob.combo"){
  #     eps0_ratio_mean = eps0[1]
  #     eps0_ratio_prob = eps0[2]
  #     eps0_diff_mean = eps0[3]
  #     eps0_diff_prob = eps0[4]
  #
  #     # if the difference is really small then we continue to P2, otherweise we stop
  #     stopping_ind <<- ifelse(min_ratio > 0,
  #                             # ratio
  #                             ifelse(seq_along(min_ratio) %in% prob_index,
  #                                    ifelse(min_ratio >= 1-eps0_ratio_prob, 0, 1), # 0 = go to P2; 1 = stop
  #                                    ifelse(min_ratio >= 1-eps0_ratio_mean, 0, 1)),
  #                             # difference
  #                             ifelse(seq_along(min_ratio) %in% prob_index,
  #                                    ifelse(abs(min_ratio) <= eps0_diff_prob, 0, 1),
  #                                    ifelse(abs(min_ratio) <= eps0_diff_mean, 0, 1)
  #                                    )
  #     )
  #   } else{
  #     eps0_ratio = eps0[1]
  #     eps0_diff = eps0[2]
  #     # message("\nTOLERANCES ARE: ", eps0_ratio, "\n", eps0_diff, "\n")
  #
  #     stopping_ind <<- ifelse(min_ratio > 0,
  #                             ifelse(min_ratio >= 1-eps0_ratio, 0, 1),
  #                             ifelse(abs(min_ratio)<=eps0_diff, 0, 1)) # if the difference is really small then we continue to P2, otherweise we stop
  #     # if bigger than 0.9 then they are similar
  #     # if bigger than 0.9 then curves are similar and we go to phase2.
  #     # if ratio >= 0.9 then ratio_stop_ind = 1 (aka we stop b/c we know which treatment to pick; otherwise ind = 0 aka we keep going)
  #   }
  #   stopped_indices <<- which(stopping_ind == 1)
  #   Ratio_Stopping_Ind = stopping_ind
  #   Ratio = min_ratio
  #   # View(stopping_ind)
  #   # View(min_ratio)
  #   # checking_stop_ind_dataset <<- cbind(subj = 1:length(stopping_ind),
  #   #                                     predd$AUS[[1]],
  #   #                                     predd$AUS[[1]]/predd$AUS[[2]],
  #   #                                     predd$Ratio[,1],
  #   #                                     predd$AUS[[2]],
  #   #                                     predd$AUS[[2]]/predd$AUS[[1]],
  #   #                                     predd$Ratio[,2],
  #   #                                     min_ratio,
  #   #                                     stopping_ind)
  #   # View(checking_stop_ind_dataset)
  #   predd_surv[["Ratio"]] <<- Ratio
  #   predd_surv[["Stopping_Ind"]] <<- Ratio_Stopping_Ind
  # } else{
  #   # print(sprintf("Phase is %s -- so we replace Ratio_Stopping_Ind and Ratio with 99", Phase))
  #   Ratio_Stopping_Ind = rep(99, length(optTx))
  #   Ratio = rep(99, length(optTx))
  # }

  # print("LAST TEST")
  return( new(Class = "Optimal",
              "optimalTx" = txLevels[optTx],
              "optimalY" = optSv,
              "type" = crit,
              "NumTrts" = num_trts, #NumTrts,
              # "NonOpt_Opt_Ratio" = Ratio,
              "Ratio_Stopping_Ind" = stop_ind)) #Ratio_Stopping_Ind) )
}

###############################################################################################
###############################################################################################
###############################################################################################

# Create a function to calculate the area under a stepwise function
.calculate_area_under_curve <- function(time_points, survival_probabilities) {
  # print(time_points)
  # print(survival_probabilities)
  tp <- length(time_points)
  # print(tp)
  area <- 0
  # print("starting")
  for (i in 1:(tp - 1)) {
    width <- time_points[i + 1] - time_points[i]
    # print("width")
    # print(width)
    height <- survival_probabilities[i]
    # print("height")
    # print(height)
    # print(area)
    area <- area + width * height
    # print(area)
    # message("area: ", area)
  }
  return(area)
}

#-------------------------------------------------------------------------------
# internal function to identify the maximum value of input vector x
#  @param x a vector of values
#  @param tieMethod a character indicating the method to be used to
#    breaks {first, random, NA}
#-------------------------------------------------------------------------------
# function returns a single numeric or NA
#-------------------------------------------------------------------------------
.whichMax <- function(x, tieMethod) {
  # print("class_SurvRF.R: .whichMax")
  ind <- which(x >= {max(x) - 1e-8}) # for ties
  if (length(x = ind) == 1L) return( ind )
  if (tieMethod == "first") return( ind[1L] )
  if (tieMethod == "random") return( resample(x = ind, size = 1L) )
  return( NA )
}

.whichMin <- function(x, tieMethod) {
  # print("class_SurvRF.R: .whichMin")
    ind <- which(x <= {min(x) + 1e-8})
    if (length(x = ind) == 1L) return(ind)
    if (tieMethod == "first") return(ind[1L])
    if (tieMethod == "random") return(resample(x = ind, size = 1L))
    return(NA)
}

.whichExtremum <- function(x, tieMethod, findMax = TRUE) {
  if (findMax) {
    # if (Phase == 2 | Phase == "CR"){
      # print("Finding max")
    # }
    extremum_val <- max(x) - 1e-8  # Find the maximum value
    ind <- which(x >= extremum_val)  # Find the indices where the maximum value occurs
  } else {
    # if (Phase == 2 | Phase == "CR"){
      # print("Finding min")
    # }
    extremum_val <- min(x) + 1e-8  # Find the minimum value
    ind <- which(x <= extremum_val)  # Find the indices where the minimum value occurs
  }
  if (length(ind) == 1L) return( ind )
  if (tieMethod == "first") return( ind[1L] )
  if (tieMethod == "random") return( resample(ind, size = 1L) )
  return( NA )
}

