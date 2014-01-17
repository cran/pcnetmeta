ci.plot <-
function(summary.stat,trtname,graphtitle=""){
  par(tck=-0.02,mgp=c(1.5,0.5,0),mar=c(5,4,1,4),cex=0.85)
  ## set parameters
  r.probt<-grep("probt",row.names(summary.stat))
  n.tr<-length(r.probt)
  xx<-1:n.tr

  ## set default arguments
  if(missing(trtname)){
    trtname<-paste("trt",1:n.tr,sep="")}else{
      if(length(trtname)!=n.tr) stop("the length of trtname is not equal to treatment number.")
    }

  ## plot
  max.probt<-max(summary.stat[r.probt,7])
  min.probt<-min(summary.stat[r.probt,3])
  graph.range<-max.probt-min.probt
  incr<-signif(0.2*graph.range,1)
  y.lim<-c(min.probt-0.1*graph.range,max.probt+0.1*graph.range)

  plot(xx,summary.stat[r.probt,5],xlim=c(0.5,n.tr+0.5),ylim=y.lim,type="n",xaxs="i",
       yaxs="i",xaxt="n",yaxt="n",lwd=2,ylab="Event Rate",xlab=graphtitle,cex=1)
  for(i in 1:n.tr){
    lines(c(xx[i],xx[i]),summary.stat[r.probt[i],c(3,7)],lty=1,lwd=2,col="green")
    points(xx[i],summary.stat[r.probt[i],5],pch=19,cex=1,col="red")
    points(xx[i],summary.stat[r.probt[i],3],pch=8,cex=1,col="blue")
    points(xx[i],summary.stat[r.probt[i],7],pch=8,cex=1,col="blue")
  }
  axis(1,at=seq(1,n.tr,1),labels=trtname,tick=TRUE)
  axis(2,at=seq(0,1,incr),labels=seq(0,1,incr),tick=TRUE)
  abline(h=seq(0,1,incr),lty=4,lwd=0.5)
}
