# Verify input 'endPoint'
#
# methods are not exported and are only for internal convenience
#
# ensures that 'endPoint' is one of {'CR', 'RE', 'MC'}. 
#
# successful methods return the original character possibly modified to be 
# all upper case.
#
setGeneric(name = ".VerifyEndPoint",
           def = function(endPoint, ...) { 
             standardGeneric(".VerifyEndPoint") 
           })

# the default method generates an error
setMethod(f = ".VerifyEndPoint",
          signature = c(endPoint = "ANY"),
          definition = function(endPoint, ...) { 
            stop("endPoint must be one of ",
                 "{'CR', 'RE', 'MC'}", 
                 call. = FALSE)
          })

setMethod(f = ".VerifyEndPoint",
          signature = c(endPoint = "character"),
          definition = function(endPoint, ...) { 
            
            endPoint <- toupper(x = endPoint)
            
            if (endPoint %in% c('CR', 'RE', 'MC'))
              return( endPoint )
            
            stop("endPoint must be one of ",
                 "{'CR', 'RE', 'MC'}", 
                 call. = FALSE)
            
          })
