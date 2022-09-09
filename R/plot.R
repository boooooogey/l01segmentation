#' @title plotsegments
#' @description
#' Plots the segmentation over the data.
#' @param segments The segmentation outputed by fusedsegmentation (format = "compressed").
#' @param data The original signal vector. If NULL, only the segmentation is plotted.
#' @param label Labels for each location. The data points are colored by this vector.
#' @param ylab Y-axis label.
#' @param xlab X-axis label.
plotsegments <- function(segments,data=NULL,label=NULL,title="",ylab="",xlab="",show="fused"){
    if(!(show %in% c("fused","breakpoints"))){
        stop("Unknown data type (Options: \"fused\" - \"breakpoints\").")
    }
    if(!is.null(data)){
        if (is.null(label)){
            plot(data,ylab = ylab,xlab = xlab,main = title, col = rgb(0,0,0,alpha = 0.5),pch=19)
        }
        else{
            plot(data,col=label+2,ylab = ylab,xlab = xlab,main = title, col = rgb(0,0,0,alpha = 0.5),pch=19)
        }
    }
    else{
        plot(0, xlab="", ylab="", xlim=c(segments[,1][1], tail(segments[,2],n=1)), ylim=c(min(segments[,3]), max(segments[,3])),col="white")
    }
    if(show == "fused"){
        for(i in 1:dim(segments)[1]){
            s = segments[,1][i]
            e = segments[,2][i]
            v = segments[,3][i]
            if(s != e){
                lines(c(s,e), c(v,v),col=rgb(1,0,0,alpha = 1),lwd = 3)
                points(c(s,e),c(v,v),col=rgb(1,0,0,alpha = 1),pch=18, cex = 1)
            }
            else{
                points(s,v,col=2,pch=18, cex = 1)
            }
            if( i < dim(segments)[1]){
                lines(c(e,segments[,1][i+1]),c(v,segments[,3][i+1]) ,col=rgb(1,0,0,alpha = 1),lwd = 3)
            }
        }
    }
    if(show == "breakpoints"){
        abline(v=segments[,1],col=2,lty=2)
        abline(v=segments[,2],col=2,lty=2)
    }
}

plotgroupsegments <- function(segments,data = NULL,...){
    p = dim(segments)[2] - 2
    par(mfrow=c(p,1))
    for(i in 1:p){
        if(!is.null(data)){
            plotsegments(segments[,c(1,2,2+i)],data=data[i,],...)
        }
        else{
            plotsegments(segments[,c(1,2,2+i)],data=NULL,...)
        }
    }
    par(mfrow=c(1,1))
}
