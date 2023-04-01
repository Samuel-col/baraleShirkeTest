library(Rcpp)
library(RcppArmadillo)

sourceCpp("baraleShirke.cpp")

baraleshirke.test <- function(X1,X2,
                              depth = "Mahalanobis",
                              NIter = 1000,
                              alpha = 0.05,
                              returnDepths = F,
                              returnSamples = F){
  
  # Validar objetos ---------------------------
  
  ## X1 y X2
  if(missing(X1)) stop("Missing X1.")
  if(missing(X2)) stop("Missing X2.")
  
  if(length(base::intersect(class(X1),c("matrix","array","data.frame"))) == 0)
    stop("X1 must be of class matrix, array or data.frame.")
  if(length(base::intersect(class(X2),c("matrix","array","data.frame"))) == 0)
    stop("X2 must be of class matrix, array or data.frame.")
  
  X1 <- as.matrix(X1)
  X2 <- as.matrix(X2)
  
  if(!is.numeric(X1)) stop("X1 must only contain numeric data.")
  if(!is.numeric(X2)) stop("X2 must only contain numeric data.")
  
  if(ncol(X1)!=ncol(X2)) stop("X1 and X2 must have the same number of columns.")
  
  ## depth
  if(!is.character(depth)) 
    stop("depth must be one of the following characters: 'betaSkeleton', 'halfspace', 'L2', 'Mahalanobis', 'MahalanobisMCD', 'projection', 'potential', 'qhpeeling', 'simplicial', 'simplicialVolume', 'spatial' and 'zonoid'.")
  if(!(depth %in% c('betaSkeleton', 'halfspace', 'L2', 'Mahalanobis', 'MahalanobisMCD', 'projection', 'potential', 'qhpeeling', 'simplicial', 'simplicialVolume', 'spatial', 'zonoid'))) 
    stop("depth must be one of the following characters: 'betaSkeleton', 'halfspace', 'L2', 'Mahalanobis', 'MahalanobisMCD', 'projection', 'potential', 'qhpeeling', 'simplicial', 'simplicialVolume', 'spatial' and 'zonoid'.")
  
  ## NIter
  if(!is.numeric(NIter)) stop("NIter must be numeric.")
  NIter <- as.integer(NIter)
  
  ## alpha
  if(!is.numeric(alpha)) stop("alpha must be numeric.")
  if(alpha < 0 || alpha > 1) stop("alpha must be greater than 0 and less than 1.")
  
  ## returnDepths
  if(!is.logical(returnDepths)) stop("returnDepths must be either TRUE or FALSE.")
  
  ## returnSamples
  if(!is.logical(returnSamples)) stop("returnSamples must be either TRUE or FALSE.")
  
  
  # Correr prueba ----------------------------------
  test.results <- baraleShirkeTest(X1,X2,depth,
                                   NIter,alpha,
                                   returnDepths,
                                   returnSamples,
                                   options("width")$width-5)
  
  class(test.results) <- "bstest"
  
  return(test.results)
}

print.bstest <- function(bs.test, ...){
  cat(bs.test$Message)
}

cat("Use baraleshirke.test function to perform the test.\n")
