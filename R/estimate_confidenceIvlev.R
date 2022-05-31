#' estimate_confidenceIvlev
#' estimate confidence intervals for chocIvlev analysis
#'
#' @param mychocIvlev an object returned by \link{chocIvlev}
#' @param conf size of the confidence interval
#' @param nb_replicates number of replicates used to assess confidence intervals
#' @param ncores Number of cores used. The parallelization will take place only if OpenMP is supported (default 1)
#' @param progressbar (default TRUE) show progressbar (might be a bit slower)
#'
#'
#' @return an updated version of mychoc with two columns added to mychoc$grid which corresponds to the bounds of the confidence interval
#' @importFrom pbapply pbsapply
#' @importFrom stats quantile
#' @importFrom pcaPP cor.fk
#' @importFrom dplyr coalesce
#'
#' @export


estimate_confidenceIvlev <-
  function(mychocIvlev,
           conf = 0.95,
           nb_replicates = 500,
           ncores = 1,
           progressbar = TRUE) {
    if (! progressbar)
      getOption("pboptions")$type == "none"
    thresholds <- NA
    grid_points <- mychocIvlev$grid[, -ncol(mychocIvlev$grid)]
    years <- seq_len(length(mychocIvlev$list_ivlev))
    ivlevs <- sapply(mychocIvlev$list_ivlev, function (li) li$ivlev)

    parallel <- FALSE
    cl <- NULL #default no cluster
    if (requireNamespace("parallel", quietly = TRUE) & ncores > 1) {
      cl <- parallel::makeCluster(min(ncores,
                                      parallel::detectCores()-1))
      parallel <- TRUE
      parallel::clusterEvalQ(cl, {
        library(ks)
        library(chocR)
      })
      parallel::clusterExport(cl, list("mychocIvlev",
                                       "ivlevs",
                                       "replicatefunction"),
                              envir = environment())
    } else if (ncores > 1){
      print("package parallel should be installed to use several cores")
    }

    replicatefunction <- function (r){
      iperm <- sample.int(length(mychocIvlev$list_ivlev), replace = TRUE)
      perm_list_ivlev <- ivlevs[, ivlevs]
      perm_tau <- apply(ivlevs, 1, function(x){
        return(coalesce(cor.fk(x, years), 0))
      })
      perm_tau
    }



    thresholds <- apply(pbsapply(
      seq_len(nb_replicates),
      replicatefunction,
      cl = cl),
      1,
      quantile,
      probs = c((1 - conf) / 2, 1 - (1 - conf) / 2))




    mychocIvlev$grid$binf <- thresholds[1, ]
    mychocIvlev$grid$bsup <- thresholds[2, ]
    mychocIvlev
  }
