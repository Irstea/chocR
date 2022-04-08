#' overallNiche
#' computes the average niche over the whole study period
#'
#' @param mychoc a chocR object of realised niche \link{choc}
#' with chocRealised
#'
#' @return a data.frame with the grid and and the corresponding density
#'
#' @export
overallNiche <- function(mychoc) {
  if (class(mychoc) != "chocR")
    stop("chocRealised should be a chocR object")
  grid <- mychoc$grid[, -ncol(mychoc$grid)]

  overallData <- do.call(rbind.data.frame, mychoc$list_data)
  overallWeights <- do.call(c, mychoc$list_weights)
  niche <- dKernel(grid = as.matrix(grid),
                        obs = as.matrix(overallData),
                        probs = overallWeights,
                        rooti = mychoc$root_i)
    res <- cbind.data.frame(grid,
                     data.frame(dens = niche))

  return(res)
}
