#' plot_choc
#' Plot the result of choc analysis
#'
#' @param mychoc a list returned by \code{\link{estimate_confidence}}
#' @param ivar the two variables that should be plotted
#'
#' @section Details:
#' red means that there is a significant decrease of occurence of the association, pale red stands for non significant trend
#' pale green for non significant positive trend and green for significant positive trend
#'
#' @return a ggplot object
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 xlab
#' @importFrom ggplot2 ylab
#' @importFrom ggplot2 geom_tile
#' @importFrom ggplot2 scale_fill_manual
#' @examples
#' #example of results provided by estimate_confidence
#' data(res_confid)
#' g<-plot_choc(res_confid)
#' print(g)
#'
#' @export

plot_choc <- function(mychoc, ivar = c(1, 2)) {
  overall_data <-
    do.call("rbind", mychoc$list_data)
  mygrid <- mychoc$grid
  if (ncol(overall_data) > 2) {
    ###if we have more than 2 dimensions, we plot the diagram at the quantile of other dimensions
    lev <- NA
    apply(overall_data, 2, function(x) {
      lev <- unique(x)
      lev[round(length(lev) / 2)]
    })
    for (i in 1:ncol(mygrid)) {
      if (!i %in% ivar) {
        mygrid <- subset(mygrid, mygrid[, i] == lev[i])
      }
    }
  }

  mygrid$icolor <- NA
  mygrid$icolor <-
    ifelse(
      mygrid$tau > 0,
      ifelse(mygrid$tau > mygrid$bsup, 1, 2),
      ifelse(mygrid$tau < mygrid$binf, 4, 3)
    )
  mygrid$icolor <- as.factor(mygrid$icolor)
  ggplot(mygrid, aes(x = mygrid[, ivar[1]], y = mygrid[, ivar[2]])) +
    geom_tile(aes(fill = mygrid$icolor)) +
    scale_fill_manual(
      values = c( "green", "palegreen1", "pink1", "red"),
      guide = FALSE
    ) +
    xlab(names(overall_data[, ivar[1]])) +
    ylab(names(overall_data[, ivar[2]]))

}
