#' estimate_confidence
#' estimate confidence intervals for choc analysis
#'
#' @param mychoc a list as returned by \link{choc}
#' @param method either "perm" (default) or "kern", see details
#' @param conf size of the confidence interval
#' @param nb_replicates number of replicates used to assess confidence intervals
#' @param ncores Number of cores used. The parallelization will take place only if OpenMP is supported (default 1)
#'
#' @section Details:
#' Two methods are available: perm permutates the kernell per time step and estimates Kendall tau on permutations.
#' kern fits a kernell on the whole dataset (assuming that there is not time trend) and uses this overall kernell to
#' generate surrogate data sets on which kendall tau are estimated. Permutations is a good solution when there is seasonnality
#' within time step to preserve internal seasonality, however, it requires more time steps. kern is a good solution when there
#' is no seasonnality within time step or when the number of observations per time step is important enough.
#'
#' @return an updated version of mychoc with two columns added to mychoc$grid which corresponds to the bounds of the confidence interval
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom stats quantile
#' @importFrom pcaPP cor.fk
#' @examples
#' #retrieve results of a choc function
#' data(res_choc)
#' #here we put a low number of replicates to limit computation time
#' #res_confid <- estimate_confidence(res_choc,"perm",0.95,50)
#'
#' @export


estimate_confidence <-
  function(mychoc,
           method = "perm",
           conf = 0.95,
           nb_replicates = 500,
           ncores=1) {
    thresholds <- NA
    years <- seq_len(length(mychoc$list_data))
    grid_points <- mychoc$grid[, -ncol(mychoc$grid)]
    cholH <- mychoc$cholH
    pb <- txtProgressBar(min = 0, max = nb_replicates, style = 3)
    if (method == "kern") {
      overall_data <-
        do.call("rbind", mychoc$list_data)
      overall_weight <-
        do.call("c", mychoc$list_weights)
      overall_weight <- overall_weight / sum(overall_weight)
      nb_time_step <- length(mychoc$list_data)
      nb_obs <- nrow(overall_data)
      center <- rep(0,ncol(overall_data))
      thresholds <- apply(sapply(1:nb_replicates, function(r) {
        mock_data <- rKernel(nb_obs,
                             overall_data,
                             mychoc$cholH,
                             overall_weight)
        mock_tvar <- do.call("c",
                             lapply(seq_along(mychoc$list_weights),
                                    function(i) rep(i,
                                                    length(mychoc$list_weights[[i]]))))
        mock_list_data <- lapply(unique(mock_tvar), function(y) {
          sub_mock_data <- subset(mock_data, mock_tvar == y)
        })

        mock_weights <- lapply(mock_list_data, function(sdata){
          rep(1/nrow(sdata),nrow(sdata))
        })
        mock_dens <-
          sapply(seq_along(mock_list_data), function(i){
            sdata <- mock_list_data[[i]]
            weights <- mock_weights[[i]]
            dKernel(grid = as.matrix(grid_points),
                    obs = sdata,
                    probs = weights,
                    rooti = mychoc$root_i)
        })

        tau <- apply(mock_dens, 1, function(x)
          cor.fk(x, years))
        setTxtProgressBar(pb,r)
        tau
      }),
      1,
      quantile,
      probs = c((1 - conf) / 2, 1 - (1 - conf) / 2))
    } else if (method == "perm") {
      thresholds <- apply(sapply(1:nb_replicates, function(r){
        iperm <- sample.int(length(mychoc$list_data), replace = TRUE)
        perm_list_data <-
          mychoc$list_data[iperm]
        perm_list_weight <-
          mychoc$list_weight[iperm]
        perm_dens <-
          sapply(seq_along(perm_list_data), function(i){
            sdata <- perm_list_data[[i]]
            weights <- perm_list_weight[[i]]
            dKernel(grid = as.matrix(grid_points),
                    obs = sdata,
                    probs = weights,
                    rooti = mychoc$root_i)
          })

        perm_tau <- apply(perm_dens, 1, function(x)
          cor.fk(x, years))
        setTxtProgressBar(pb,r)
        perm_tau
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
