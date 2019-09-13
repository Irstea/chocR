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
  mygrid$color<-ifelse(mygrid$tau>0,ifelse(mygrid$tau>mygrid$bsup,1,2),ifelse(mygrid$tau<mygrid$binf,4,3))
  mygrid$color<-ifelse(mygrid$inhull,mygrid$color,0)
  mygrid$color=as.factor(mygrid$color)
  ggplot(mygrid,aes(x=mygrid[,ivar[1]],y=mygrid[,ivar[2]]))+
    geom_raster(aes(fill=color))+
    scale_fill_manual(values=c("white","green","palegreen1","pink1","red"),guide=FALSE)+
    xlab(names(overall_data[,ivar[1]]))+
    ylab(names(overall_data[,ivar[2]]))

}



