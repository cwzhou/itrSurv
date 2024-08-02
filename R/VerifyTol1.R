### LATER: NEED TO CHANGE SO THERE IS NO DIFF AND UPDATE CORRESOPNDING PARTS IN OTHER SCRIPTS
# WE LOOK AT RATIO ONLY.

# Verify input 'tol1', which is the tolerance only for Phase1
#
# methods are not exported and are for internal convenience only
#
# ensures that 'tol1' is provided if criticalValue is one of
# {'mean', 'mean.prob.combo'}.
#
# successful methods return the numeric vector tol1 or NULL.
#
setGeneric(name = ".VerifyTol1",
           def = function(tol1, ...) {
             standardGeneric(".VerifyTol1")
           })

# the default method generates an error
setMethod(f = ".VerifyTol1",
          signature = c(tol1 = "ANY"),
          definition = function(tol1, ...) {
            stop("tol1 must be numeric or NULL",
                 call. = FALSE)
          })

setMethod(f = ".VerifyTol1",
          signature = c(tol1 = "numeric"),
          definition = function(tol1, ..., criticalValue) {

            if (length(x = tol1) < 2L | length(x = tol1) > 4L){
              stop("tol1 must have length 2, 3, or 4 depending on criticalValue for Phase 1. Minimally, one for ratio tolerance and one for difference tolerance.")
            }

            # 'mean', 'prob', 'area', 'mean.prob.combo'
            if (criticalValue %in% c("mean", "prob", "area")) {
              # message(criticalValue," only needs 1 ratio tolerance value + 1 difference tolerance value.")
              if (length(x = tol1) == 2L){
                message(criticalValue, " ratio tolerance: ", tol1[1],
                        "\n",criticalValue, " diff tolerance: ", tol1[2],
                        "\nIf this is wrong: enter the correct order of 'ratio', 'difference'.")
              } else{
                stop("tol1 has length ", length(x = tol1), " but it should only need 2 for ", criticalValue)
              }
            } else if (criticalValue %in% c("mean.prob.combo")) {
              if (length(x = tol1) == 4L){
                message(criticalValue, " mean ratio tolerance: ", tol1[1],
                        "\n", criticalValue, " prob ratio tolerance: ", tol1[2],
                        "\n", criticalValue, " mean diff tolerance: ", tol1[3],
                        "\n",criticalValue, " prob diff tolerance: ", tol1[4],
                        "\nMake sure the correct order is: 'mean ratio', 'prob ratio', 'mean diff', 'prob diff'.")
              } else{
                stop("tol1 has length ", length(x = tol1), " but it should have 4 elements for ", criticalValue)
              }
            } else{
              stop("criticalValue ", criticalValue, " not defined.")
            }

            return( tol1 )
          })

setMethod(f = ".VerifyTol1",
          signature = c(tol1 = "NULL"),
          definition = function(tol1, ..., tau) {
            stop("tol1 is null...")
          })
