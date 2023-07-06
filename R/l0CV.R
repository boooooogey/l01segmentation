blockedCVpoisson=function(x, block.length=30,l.try=NULL, n.cv=10,  objective="KL", make.plot=TRUE){
  
  run.length=n.cv*block.length
  num.block.lengths=length(x)/run.length+1
  ids=rep(rep((1:n.cv), each=block.length), num.block.lengths)
  ids=ids[1:length(x)]
  if(make.plot){
    par(mfrow=c(2,1))
  }
  
  if(is.null(l.try)){
    l.try=10^( seq(0,4, 0.1))
  }
  #  stop()
  out=matrix(nrow=length(l.try), ncol=n.cv)
  
  x[x<0]=0
  segStart=0:(ceiling(length(x)/block.length)-1)*block.length-block.length+1
  segStart=rep(segStart, each=block.length)
  segEnd=segStart+2*block.length
  segStart=segStart[1:length(x)]
  segEnd=segEnd[1:length(x)]
  #handle edge cases
  ii=which(segStart<1)
  segStart[ii]=segEnd[ii]
  ii=which(segEnd>length(x))
  segEnd[ii]=segStart[ii]
  for(j in 1:n.cv){
    train.ids=(ids!=j)
    
    test.ids=(ids==j)
    
    yout=double(length(x))
    for(i in 1:length(l.try)){
      l=l.try[i]
      tmp=x[train.ids]
      #      tmp[test.ids]=mean(tmp)
      
      tmpseg=l01segmentation::fusedsegmentation(tmp, l, format = "full", objective = "poisson")
      yout[]=NA
      yout[train.ids]=tmpseg #put back
      y=x[test.ids]
   
      ypred=(yout[segStart[test.ids]]+yout[segEnd[test.ids]])/2
      #KL divergence
      offset=.01
      out[i,j]=mean(y*log((y+offset)/(ypred+offset))-y+ypred)
      
    }
  }
  
  
  #  out.m.raw=out.m=apply(out,1,mean, trim=0.05)
  out.m.raw=out.m=apply(out,1,median)
  out.se=apply(out,1,var)/sqrt(ncol(out))
  ll=log(l.try)
  sres=smooth.spline(ll, out.m, df = 5)
  out.m=sres$y
  
  plot(l.try, out.m, type="l", log="x", xlab="lambda", ylab="objective")
  points(l.try, out.m.raw)
  
  
  i.best=which.min(out.m)
  
  best.lam=l.try[i.best]
  
  # show(best.lam)
  
  tmpseg=l01segmentation::fusedsegmentation(x, best.lam, format = "compressed", objective = "poisson")
  plotsegments(tmpseg, data=x)
  
  rownames(out)=l.try
  return(list(lbest=best.lam, cv.result=out, cv.result.smoothed=out.m, lambdas=l.try, segmented=tmpseg))
}

findLambda=function(x, nSegments, lStart=10, ...){
  lUp=lStart
  lDown=lStart
  ntry=0
  tmpres=l01segmentation::fusedsegmentation(x, lambda = lUp, format = "compressed", ...)
  while(ntry<10){
    if(length(tmpres$start)>nSegments){ #too few decrease
      lUp=lUp*10
      tmpres=l01segmentation::fusedsegmentation(x, lambda = lUp, format = "compressed", ...)
    }
    else{
      break
    }
  }
  ntry=10
  tmpres=l01segmentation::fusedsegmentation(x, lambda = lDown, format = "compressed", ...)
  while(ntry<10){
    if(length(tmpres$start)<nSegments){ #too few decrease
      lDown=lDown/10
      tmpres=l01segmentation::fusedsegmentation(x, lambda = lDown, format = "compressed", ...)
    }
    else{
      break
    }
  }
  
  lCur=sqrt(lUp*lDown)

  tmpres=l01segmentation::fusedsegmentation(x, lambda = lCur, format = "compressed", ...)
  
  ntry=0
  
  
  while(ntry<10){
    lCur=sqrt(lUp*lDown)
    ntry=ntry+1
    tmpres=l01segmentation::fusedsegmentation(x, lambda = lCur, format = "compressed", ...)
  

    if(length(tmpres$start)<nSegments){ #too few decrease
        lUp=lCur 
      
    }
    if(length(tmpres$start)>nSegments){ #too many increase
   
        lDown=lCur 
      
    }

  }
  return(tmpres)
}
