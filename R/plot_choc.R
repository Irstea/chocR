#' Plot the result of choc analysis
#' red means that there is a significant decrease of occurence of the association, pale red stands for non significant trend
#' pale green for non significant positive trend and green for significant positive trend
#'
#' @param mychoc a list returned by \code{\link{estimate_confidence}}
#' @param ivar the two variables that should be plotted
#'
#' @return a ggplot object
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom grDevices chull
#' @importFrom ggplot2 xlab
#' @importFrom ggplot2 ylab
#' @importFrom ggplot2 geom_raster
#' @importFrom sp point.in.polygon
#' @importFrom ggplot2 scale_fill_manual
#' @examples
#' #example of results provided by estimate_confidence
#' data(res_confid)
#' g<-plot_choc(res_confid)
#' print(g)
#'
#' @export

plot_choc<-function(mychoc,ivar=c(1,2)){
  overall_data <-
    do.call("rbind", lapply(mychoc$kernels, function(x)
      x$x))
  mygrid<-mychoc$grid
  if (ncol(overall_data)>2){ ###if we have more than 2 dimensions, we plot the diagram at the quantile of other dimensions
    lev <- NA
    apply(overall_data,2,function(x) {
      lev<-unique(x)
      lev[round(length(lev)/2)]
    })
    for (i in 1:ncol(mygrid)){
      if (!i%in%ivar){
        mygrid<-subset(mygrid,mygrid[,i]==lev[i])
      }
    }
    }

  mychull <- chull(overall_data[,ivar])
  mychull <- c(mychull,mychull[1])

  mygrid$inhull <- point.in.polygon(mygrid[,ivar[1]],mygrid[,ivar[2]],pol.x=overall_data[mychull,1],pol.y=overall_data[mychull,2])
  mygrid$icolor<-ifelse(mygrid$tau>0,ifelse(mygrid$tau>mygrid$bsup,1,2),ifelse(mygrid$tau<mygrid$binf,4,3))
  mygrid$icolor<-ifelse(mygrid$inhull,mygrid$icolor,0)
  mygrid$icolor=as.factor(mygrid$icolor)
  ggplot(mygrid,aes(x=mygrid[,ivar[1]],y=mygrid[,ivar[2]]))+
    geom_raster(aes(fill=icolor))+
    scale_fill_manual(values=c("white","green","palegreen1","pink1","red"),guide=FALSE)+
    xlab(names(overall_data[,ivar[1]]))+
    ylab(names(overall_data[,ivar[2]]))

}



