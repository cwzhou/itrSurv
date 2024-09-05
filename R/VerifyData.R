# Verify input 'data'
#
# methods are not exported and are for internal convenience only
#
# ensures that 'data' is provided as a data.frame or a matrix with named
# columns and that the object does not contain NaN values.
#
# for CR datasets: these are one-row per subject with columns, sorted in ascending order by yName:
#                  id, delta, delta_j, covariates, treatment assignment, obs_time
# for RE datasets: these are in andersen-gill (more than one row per id) format with the columns:
#                  id, start, stop, delta_r, delta_d, covariates, treatment assignment
#
# successful methods return a data.frame object containing the data
#
setGeneric(name = ".VerifyData",
           def = function(data, ...) { standardGeneric(".VerifyData") })

# the default method generates an error
setMethod(f = ".VerifyData",
          signature = c(data = "ANY"),
          definition = function(data, ..., yName, epName, endPoint, idName) {
              stop("data must be a data.frame or a matrix with named columns",
                   call. = FALSE)
            })

setMethod(f = ".VerifyData",
          signature = c(data = "data.frame"),
          definition = function(data, ..., yName, epName, endPoint, idName) {

              if (any(sapply(X = data, FUN = is.nan))) {
                stop("data cannot include NaN values", call. = FALSE)
              }

            # 'yName' is the column name for observed time
            if (endPoint == "CR") {
              yName <- .VerifyYName(yName = yName, data = data)
              vector_y = data %>% dplyr::select(!!sym(yName)) %>% unlist()
              # sorted = sort(vector)
              if (is.unsorted(vector_y)) {
                stop("For CR endpoint, data is not sorted in ascending order by yName. Data must be sorted before added as an argument for itrSurv (CR).")
              }
            }

            if (endPoint == "RE"){
              message("Endpoint: RE")

              # ensure that 'idName' is provided as a character and
              # that the provided name is present in 'data'. If 'idName' is appropriate,
              # the object returned is the original input without modification.
              idName <- .VerifyIdName(idName = idName, data = data)

              # ensure that 'epName' is provided as a character or character vector and
              # that the provided names are present in 'data'. For RE: this input defines the
              # dataset for the endpoint Phase analysis. For CR: this is needed for 'gray_cr' test.
              # If 'epName' is appropriate, the object returned is the original input without modification.
              epName <- .VerifyEpName(epName = epName, data = data, endPoint = endPoint)

              # Determining Phase 1 and Phase 2 Datasets (in case they differ like in RE)
              # message("First, we identify failure dataset")
              # epName = "status" #equal to 1 if recurrent event.
              data_surv = data %>% filter(!!sym(epName) == 0) #filter to NON RE (one row per id)
              # message("Next, we identify recurrent event dataset")
              data_ep = data # entire dataset

              # sample size
              nsamp = data %>% distinct(!!sym(idName)) %>% nrow()

              # id vector for endpoint data
              id_ep = data %>% dplyr::select(!!sym(idName)) #do we need this unlist? %>% unlist()

              if ( nrow(data_surv) != nsamp) {
                View(data_surv)
                View(data_ep)
                View(id_ep)
                print(sprintf("Sample size: %s", nsamp))
                print(sprintf("Survival dataset sample size based on filtering to non-recurrent events (epName == 0): %s", nrow(data_surv)))
                message("Please fix your dataset before continuing. Failure dataset and Endpoint datasets are printed for your convenience, as well as sample size based on idName variable.")
                stop("The total number of individuals does not match the length of survival dataset. \n Make sure your dataset has a censored row using tau if there are only recurrent events for that individual \n (i.e., (last TStop time, tau] as (TStart, TStop] with Status == 0 and Status_D == 0 (censored)).")
              }

              } else{ # CR
                idName = NULL
                # for CR, the dataset is the same for both phase 1 and phase 2
                message("Phase 1 and Phase 2 datasets are the same.")
                data_surv = data_ep = data
                id_ep = c(1:nrow(data_ep))
              }

              return( list(data_surv, data_ep, id_ep) )
            })

setMethod(f = ".VerifyData",
          signature = c(data = "matrix"),
          definition = function(data, ..., yName, epName, endPoint, idName) {

              if (is.null(x = colnames(x = data))) {
                stop("if a matrix, data must include column headers",
                     call. = FALSE)
              }

              return( .VerifyData(data = as.data.frame(x = data),
                                  yName = yName, epName = epName,
                                  endPoint = endPoint, idName = idName, ...) )
            })
