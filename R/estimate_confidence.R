#' estimate confidence intervals for choc analysis
#' Two methods are available: perm permutates the kernell per time step and estimates Kendall tau on permutations.
#' kern fits a kernell on the whole dataset (assuming that there is not time trend) and uses this overall kernell to
#' generate surrogate data sets on which kendall tau are estimated. Permutations is a good solution when there is seasonnality
#' within time step to preserve internal seasonality, however, it requires more time steps. kern is a good solution when there
#' is no seasonnality within time step or when the number of observations per time step is important enough.
#'
#' @param mychoc a list as returned by \link{choc}
#' @param method either "perm" (default) or "kern", see details
#' @param conf size of the confidence interval
#' @param nb_replicates number of replicates used to assess confidence intervals
#'
#' @return an updated version of mychoc with two columns added to mychoc$grid which corresponds to the bounds of the confidence interval
#' @importFrom ks kde
#' @importFrom ks rkde
#' @importFrom stats quantile
#' @importFrom Kendall MannKendall
#' @examples
#' #generate artificial data set
#' #two time series measured on 40 time steps with 365 observations per time step.
#' #the two series follow a multinormal time series with a tend on means and a constant
#' #covariance matrix
#' library(MASS)
#' library(ks)
#' tvar <- rep(1:40,times=100) #times steps
#' meansX <- tvar/40 #trend on 1st variable
#' meansY <- -0.5*tvar/40 #trend on 2nd variable
#' sigma <- matrix(c(1,.1,.1,1),2,2) #covariance matrix
#' values <- t(apply(cbind(meansX,meansY),1,function(mu) mvrnorm(1,mu,sigma))) #generate the values
#' H<-Hpi #choose the default bandwith
#' res <- choc(values,H,tvar)
#' #here we put a low number of replicates to limit computation time
#' res_confid <- estimate_confidence(res,"perm",0.95,50)
#'
#' @export


estimate_confidence <-
  function(mychoc,
           method = "perm",
           conf = 0.95,
           nb_replicates = 500) {
    thresholds <- NA
    grid_points <- mychoc$grid[, -ncol(mychoc$grid)]
    H <- mychoc$kernels[[1]]$H
    if (method == "kern") {
      overall_data <-
        do.call("rbind", lapply(mychoc$kernels, function(x)
          x$x))
      overall_kernel <- kde(overall_data, H = H)
      nb_time_step <- length(mychoc$kernels)
      nb_obs <- nrow(overall_data)
      thresholds <- apply(replicate(nb_replicates, {
        #sample nb_obs observations from overall_kernell
        mock_data <- data.frame(rkde(nb_obs, overall_kernel))
        mock_tvar <-
          as.vector(sapply(1:length(mychoc$kernels), function(i)
            rep(i, nrow(
              mychoc$kernels[[i]]$x
            ))))
        mock_kernels <- lapply(unique(mock_tvar), function(y) {
          sub_mock_data <- subset(mock_data, mock_tvar == y)
          kde(sub_mock_data, H = H)
        })
        mock_dens <-
          sapply(mock_kernels, function(kern)
            kde(
              x = kern$x,
              H = kern$H,
              eval.points = grid_points
            )$estimate)
        tau <- apply(mock_dens, 1, function(x)
          MannKendall(x)$tau)
      }),
      1,
      quantile,
      probs = c((1 - conf) / 2, 1 - (1 - conf) / 2))
    } else if (method == "perm") {
      thresholds <- apply(replicate(nb_replicates, {
        perm_kernels <-
          mychoc$kernels[sample.int(length(mychoc$kernels), replace = TRUE)]
        perm_dens <-
          sapply(perm_kernels, function(kern)
            kde(
              x = kern$x,
              H = H,
              eval.points = grid_points
            )$estimate)
        perm_tau <- apply(perm_dens, 1, function(x)
          MannKendall(x)$tau)
      }),
      1,
      quantile,
      probs = c((1 - conf) / 2, 1 - (1 - conf) / 2))
    } else{
      stop("wrong method")
    }
    mychoc$grid$binf <- thresholds[1, ]
    mychoc$grid$bsup <- thresholds[2, ]
    mychoc
  }
