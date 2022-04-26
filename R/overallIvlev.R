#' overallIvlev
#' computes the average preferenum over the whole study period
#'
#' @param chocRealised a chocR object of realised niche \link{choc}
#' @param chocAvailable a chocR object of available niche consistent \link{choc}
#' with chocRealised
#'
#' @return a data.frame with the grid and an ivlev index
#' @importFrom ks kde
#'
#' @export
overallIvlev <- function(chocRealised, chocAvailable) {
  if (class(chocRealised) != "chocR")
    stop("chocRealised should be a chocR object")
  grid <- chocRealised$grid[, -ncol(chocRealised$grid)]
  if (class(chocAvailable) != "chocR")
    stop("chocAvailable should be a chocR object")
  if (!all(grid == chocAvailable$grid[, -ncol(chocAvailable$grid)]))
    stop("the choc objects grids should be similar")
  overallAvailabeData <- do.call(rbind.data.frame, chocAvailable$list_data)
  overallAvailabeWeights <- do.call(c, chocAvailable$list_weights)
  overallRealisedData <- do.call(rbind.data.frame, chocRealised$list_data)
  overallRealisedWeights <- do.call(c, chocRealised$list_weights)


  realised <-  kde(eval.points = as.matrix(grid),
                   x = as.matrix(overallRealisedData),
                   w = overallRealisedWeights / sum(overallRealisedWeights),
                   H = chocRealised$H,
                   binned = TRUE)$estimate

  avalaible <-  kde(eval.points = as.matrix(grid),
                   x = as.matrix(overallAvailabeData),
                   w = overallAvailabeWeights / sum(overallAvailabeWeights),
                   H = chocAvailable$H,
                   binned = TRUE)$estimate

    res <- cbind.data.frame(grid,
                     data.frame(ivlev = (realised - avalaible) /
                                  (realised + avalaible)))

  return(res)
}
