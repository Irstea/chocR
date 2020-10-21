#' estimate_confidenceIvlev
#' estimate confidence intervals for chocIvlev analysis
#'
#' @param mychocIvlev an object returned by \link{chocIvlev}
#' @param conf size of the confidence interval
#' @param nb_replicates number of replicates used to assess confidence intervals
#' @param ncores Number of cores used. The parallelization will take place only if OpenMP is supported (default 1)
#'
#'
#' @return an updated version of mychoc with two columns added to mychoc$grid which corresponds to the bounds of the confidence interval
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom stats quantile
#' @importFrom pcaPP cor.fk
#'
#' @export


estimate_confidenceIvlev <-
  function(mychocIvlev,
           conf = 0.95,
           nb_replicates = 500,
           ncores=1) {
    thresholds <- NA
    pb <- txtProgressBar(min = 0, max = nb_replicates, style = 3)
    grid_points <- mychocIvlev$grid[, -ncol(mychocIvlev$grid)]
    years <- seq_len(length(mychocIvlev$list_ivlev))
    ivlevs <- sapply(mychocIvlev$list_ivlev, function (li) li$ivlev)
    thresholds <- apply(sapply(1:nb_replicates, function(r){
        iperm <- sample.int(length(mychocIvlev$list_ivlev), replace = TRUE)
        perm_list_ivlev <- ivlevs[, ivlevs]
        perm_tau <- apply(ivlevs, 1, function(x)
          cor.fk(x, years))
        setTxtProgressBar(pb,r)
        perm_tau
      }),
      1,
      quantile,
      probs = c((1 - conf) / 2, 1 - (1 - conf) / 2))

    mychocIvlev$grid$binf <- thresholds[1, ]
    mychocIvlev$grid$bsup <- thresholds[2, ]
    mychocIvlev
  }
