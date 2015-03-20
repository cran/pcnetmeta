ci.plot <- function(nma.obj,alphabetic=TRUE){
  par(tck=-0.02,mgp=c(1.5,0.5,0),mar=c(3,4,1,4),cex=0.85)
  ## set parameters
  if(!is.null(nma.obj$AbsoluteRisk)){
    ci<-nma.obj$AbsoluteRisk$Median_CI
    armparam<-"Absolute Risk"
  }
  if(!is.null(nma.obj$TrtEffect)){
    ci<-nma.obj$TrtEffect$Median_CI
    armparam<-"Treatment Effect"
  }
  if(!is.null(nma.obj$HazardRate)){
    ci<-nma.obj$HazardRate$Median_CI
    armparam<-"Hazard Rate"
  }
  n.tr<-dim(ci)[1]
  xx<-1:n.tr
  trtname<-rownames(ci)

  med<-low<-upp<-numeric(n.tr)
  for(i in 1:n.tr){
    str<-ci[i,1]
    split1<-strsplit(str,split=" \\(")
    med[i]<-as.numeric(split1[[1]][1])
    str2<-split1[[1]][2]
    split2<-strsplit(str2,split=", ")
    low[i]<-as.numeric(split2[[1]][1])
    upp[i]<-as.numeric(gsub("\\)","",split2[[1]][2]))
  }

  if(alphabetic){
    od<-order(trtname)
    trtname<-trtname[od]
    med<-med[od]
    low<-low[od]
    upp<-upp[od]
  }

  ## plot
  graph.range<-max(upp)-min(low)
  incr<-signif(0.2*graph.range,1)
  y.lim<-c(min(low)-0.1*graph.range,max(upp)+0.1*graph.range)

  plot(xx,med,xlim=c(0.5,n.tr+0.5),ylim=y.lim,type="n",xaxs="i",
       yaxs="i",xaxt="n",yaxt="n",lwd=2,ylab=armparam,xlab="",cex=1)
  abline(h=seq(0,1,incr),lty=4,lwd=0.5,col="grey")
  for(i in 1:n.tr){
    lines(c(xx[i],xx[i]),c(low[i],upp[i]),lty=1,lwd=2,col="black")
    points(xx[i],med[i],pch=15,cex=1,col="black")
    points(xx[i],low[i],pch=8,cex=1,col="black")
    points(xx[i],upp[i],pch=8,cex=1,col="black")
  }
  axis(1,at=seq(1,n.tr,1),labels=trtname,tick=TRUE)
  axis(2,at=seq(0,1,incr),labels=seq(0,1,incr),tick=TRUE)
}
