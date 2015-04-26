rank.prob <- function(nma.obj){
  if(is.null(nma.obj$TrtRankProb)) stop("users do not specify rank.prob in the argument param of the function which produces nma.obj.")
  if(!is.null(nma.obj$TrtRankProb)){
    rank.prob<-nma.obj$TrtRankProb
  }
  ntrt<-dim(rank.prob)[1]
  trtname<-rownames(rank.prob)
  plot(9/20/ntrt,0.5,col="white",xlim=c(0,1),ylim=c(0,1),xaxt="n",yaxt="n",xlab="Treatments",ylab="Rank Probabilities",frame.plot=FALSE)
  axis(2,at=seq(0,1,0.1),labels=TRUE)
  axis(1,at=9/20/ntrt+seq(from=0,by=1/ntrt,length.out=ntrt),labels=trtname,tick=FALSE)
  rgb.val<-seq(from=0,by=0.8/(ntrt-1),length.out=ntrt)
  for(i in 1:ntrt){
    trt.i.rank.prob<-rank.prob[i,]
    cum.prob<-c(0,cumsum(trt.i.rank.prob))
    xleft.i<-(i-1)/ntrt
    xright.i<-i/ntrt-1/(10*ntrt)
    for(j in 1:ntrt){
      rect(xleft=xleft.i,ybottom=cum.prob[j],xright=xright.i,ytop=cum.prob[j+1],col=rgb(rgb.val[j],rgb.val[j],rgb.val[j]),border="black")
    }
  }
}