#' choc
#' carry out a choc analysis on a multivariate time series
#'
#' @param mydata a data frame or matrix with one column per time series and one row per observation (event). Marginal distributions are assumed to follow a normal distribution
#' @param H either a function from library ks to estimate a bandwith matrix, or directly a bandwith matrix. Exemples of function \code{\link[ks]{Hpi}} or \code{\link[ks]{Hscv}}...
#' @param timevar a vector specifying the time step for each observation. Several observations per time step are required. Observations of a same time step are assumed to be replicates.
#' @param resolution grid resolution on which the densities of probability will be computed. The resolution does not affect the result: high resolution increases the resolution of final diagrams but increases computation time and memory usage. The number of observations per time step should be roughly similar
#' @param ncores Number of cores used. The parallelization will take place only if OpenMP is supported.
#'
#' @return a list with 3 elements: \enumerate{
#' \item the list of data per time step
#' \item a dafaframe with one row per time series and a column tau that gives the tau of Kendall trend test
#' \item the cholesky decomposition of H
#' }
#' @importFrom mvnfast dmvn
#' @importFrom Kendall MannKendall
#' @examples
#' #generate artificial data set
#' #two time series measured on 40 time steps with 365 observations per time step.
#' #the two series follow a multinormal time series with a tend on means and a constant
#' #covariance matrix
#' library(MASS)
#' require(ks)
#' tvar <- rep(1:40,times=100) #times steps
#' meansX <-tvar/40 #trend on 1st variable
#' meansY <- -0.5*tvar/40 #trend on 2nd variable
#' sigma <- matrix(c(1,.1,.1,1),2,2) #covariance matrix
#' values <- t(apply(cbind(meansX,meansY),1,function(mu) mvrnorm(1,mu,sigma))) #generate the values
#' H <- Hpi #choose the default bandwith
#' res_choc <- choc(values,H,tvar)
#'
#' @export
choc <- function(mydata, H, timevar, resolution = 100,ncores=1) {
  if (class(H) == "function") {
    H <- H(mydata)
  } else if (class(H) != "matrix") {
    stop ("H is not valid, should be a function or a matrix")
  }
  cholH <- chol(H)
  #create a list of data per timestep
  list_data <- lapply(unique(timevar), function(y) {
    sub_mydata <- subset(mydata, timevar == y)
  })
  names(list_data) <- unique(timevar)

  #build a grid on which densities of probability will be computed
  grid <-
    expand.grid(as.list(as.data.frame(apply(mydata, 2, function(x)
      seq(min(x), max(x), length.out = resolution)))))
  names(grid) <- colnames(mydata)

  #compute density for each point of the grid and each time step
  densprob <-
    sapply(list_data, function(sdata){
      rowMeans(apply(sdata,1,function(points){
        dmvn(as.matrix(grid),points,sigma=cholH,isChol=TRUE,ncores=ncores)
      }))
    })

  #compute tau kendall
  grid$tau <- apply(densprob, 1, function(x)
    MannKendall(x)$tau)

  return(list(list_data=list_data,grid=grid,cholH=cholH))
}
