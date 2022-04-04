#' choc
#' carry out a choc analysis on a multivariate time series
#'
#' @param mydata a data frame or matrix with one column per time series and one row per observation (event). Marginal distributions are assumed to follow a normal distribution
#' @param H either a function from library ks to estimate a bandwith matrix, or directly a bandwith matrix. Exemples of function \code{\link[ks]{Hpi}} or \code{\link[ks]{Hscv}}...
#' @param timevar a vector specifying the time step for each observation. Several observations per time step are required. Observations of a same time step are assumed to be replicates.
#' @param weights weights of the observations
#' @param resolution grid resolution on which the densities of probability will be computed. The resolution does not affect the result: high resolution increases the resolution of final diagrams but increases computation time and memory usage. The number of observations per time step should be roughly similar
#' @param ncores Number of cores used. The parallelization will take place only if OpenMP is supported.
#'
#' @return a chocR object, i.e. a list with 5 elements: \enumerate{
#' \item list_data the list of data per time step
#' \item grid a dafaframe with one row per time series and a column tau that gives the tau of Kendall trend test
#' \item cholH the cholesky decomposition of H
#' \item list_weights the weights of observation
#' \item root_i a transformation of inv(trimatu(cholH))
#' }
#' @importFrom pcaPP cor.fk
#' @importFrom grDevices chull
#' @importFrom pracma inpolygon
#' @import Rcpp
#' @examples
#' #generate artificial data set
#' #two time series measured on 40 time steps with 100 observations per time step.
#' #the two series follow a multinormal time series with a tend on means and a constant
#' #covariance matrix
#' if (require(MASS) & require(ks)){
#'   tvar <- rep(1:40,times=100) #times steps
#'   meansX <-tvar/40 #trend on 1st variable
#'   meansY <- -0.5*tvar/40 #trend on 2nd variable
#'   sigma <- matrix(c(1,.1,.1,1),2,2) #covariance matrix
#'   values <- t(apply(cbind(meansX,meansY),1,function(mu) mvrnorm(1,mu,sigma))) #generate the values
#'   H <- Hpi #choose the default bandwith
#'   res_choc <- choc(values,H,tvar)
#' }
#'
#' @export
choc <- function(mydata, H, timevar, weights= NULL, resolution = 100,ncores=1) {
  if (is.null(weights))
    weights <- rep (1/nrow(mydata), nrow(mydata))

  if (class(H) == "function") {
    H <- H(mydata)
  } else if (class(H) != "matrix") {
    stop ("H is not valid, should be a function or a matrix")
  }
  cholH <- chol(H)
  root_i=get_root_i(cholH)


  #create a list of data per timestep
  list_data <- lapply(unique(timevar), function(y) {
    sub_mydata <- subset(mydata, timevar == y)
    sub_mydata
  })
  names(list_data) <- unique(timevar)

  list_weights <- lapply(unique(timevar), function(y) {
    sub_weights <- subset(weights, timevar == y)
    sub_weights <- sub_weights / sum(sub_weights)
    sub_weights
  })
  names(list_weights) <- unique(timevar)

  #get the convex hull
  liste <- chull(cbind(mydata))
  hull <- mydata[liste, ]




  #build a grid on which densities of probability will be computed
  grid <-
    expand.grid(as.list(as.data.frame(apply(mydata, 2, function(x)
      seq(min(x), max(x), length.out = resolution)))))
  names(grid) <- colnames(mydata)
  grid <- grid[inpolygon(grid[, 1], grid[, 2], hull[, 1], hull[, 2]), ]

  #compute density for each point of the grid and each time step
  densprob <-
    sapply(seq_along(list_data), function(i){
      sdata <- list_data[[i]]
      weights <- list_weights[[i]]
      dKernel(grid = as.matrix(grid),
              obs = sdata,
              probs = weights,
              rooti = root_i)
      })

  #compute tau kendall
  years <- seq_len(length(list_data))
  grid$tau <- apply(densprob, 1, function(x)
    cor.fk(x, years))

  res <- list(list_data=list_data,grid=grid,cholH=cholH,
              list_weights=list_weights,
              root_i=root_i)
  class(res) <- "chocR"
  return(res)
}
