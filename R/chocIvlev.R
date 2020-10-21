#' chocIvlev
#' carry out a choc analysis by comparing a chocR of realised niche and a chocR
#' of available niche
#'
#' @param chocRealised a chocR object of realised niche \link{choc}
#' @param chocAvailable a chocR object of available niche consistent \link{choc}
#' with chocRealised
#'
#' @return a list with 3 elements: \enumerate{
#' \item chocRealised
#' \item chocAvailable
#' \item list_ivlev the list of ivelev per time step computed on grid
#' \item grid the Kendall tau on the grid
#' }
#'
#' @export
chocIvlev <- function(chocRealised, chocAvailable) {
  if (class(chocRealised) != "chocR")
    stop("chocRealised should be a chocR object")
  grid <- chocRealised$grid[, -ncol(chocRealised$grid)]
  if (class(chocAvailable) != "chocR")
    stop("chocAvailable should be a chocR object")
  if (!all(grid == chocAvailable$grid[, -ncol(chocAvailable$grid)]))
    stop("the choc objects grids should be similar")
  list_ivlev <- lapply(seq_along(chocRealised$list_data) , function(i){
    cbind.data.frame(grid,
                     data.frame(ivlev = computeIvlev(chocRealised[[i]],
                                                     chocAvailable[[i]])
                                ))
  })

  years <- seq_len(length(list_ivlev))
  grid$tau <- apply(sapply(list_ivlev, function (li) li$ivlev),
                    1,
                    function(x) cor.fk(x, years))

  res <- list(chocRealised = chocRealised,
              chocAvailable = chocAvailable,
              list_ivlev = list_ivlev,
              grid = grid)

  class(res) <- "chocIvlev"

  return(res)
}
