contrast.plot <- function(nma.obj,effect.size,reference,left.margin=0.5,width=10,height){
  if(missing(nma.obj)) stop("nma.obj is not specified.")
  if(missing(effect.size)){
    if(!is.null(nma.obj$OddsRatio) | !is.null(nma.obj$LogOddsRatio)){
      effect.size<-"OR"
    }else{
      if(!is.null(nma.obj$RelativeRisk) | !is.null(nma.obj$LogRelativeRisk)){
        effect.size<-"RR"
      }else{
        if(!is.null(nma.obj$RiskDifference)){
          effect.size<-"RD"
        }else{
          effect.size<-""
        }
      }
    }
    if(!is.null(nma.obj$EffectDiff)) effect.size<-"diff"
    if(!is.null(nma.obj$RateRatio) | !is.null(nma.obj$LogRateRatio)) effect.size<-"ratio"
  }
  if(effect.size=="") stop("users do not specify a relative effect size in the argument param of the function which produces nma.obj.")
  if(!is.element(effect.size,c("OR","LOR","RR","LRR","RD","diff","ratio","logratio"))) stop("the input effect size is wrong.")

  if(effect.size=="OR"){
    if(!is.null(nma.obj$OddsRatio)){
      ef<-nma.obj$OddsRatio$Median_CI
      trans<-"id"
    }else{
      if(!is.null(nma.obj$LogOddsRatio)){
        ef<-nma.obj$LogOddsRatio$Median_CI
        trans<-"exp"
      }else{
        stop("the specified effect size is not estimated in nma.obj.")
      }
    }
  }

  if(effect.size=="LOR"){
    if(!is.null(nma.obj$LogOddsRatio)){
      ef<-nma.obj$LogOddsRatio$Median_CI
      trans<-"id"
    }else{
      if(!is.null(nma.obj$OddsRatio)){
        ef<-nma.obj$OddsRatio$Median_CI
        trans<-"log"
      }else{
        stop("the specified effect size is not estimated in nma.obj.")
      }
    }
  }

  if(effect.size=="RR"){
    if(!is.null(nma.obj$RelativeRisk)){
      ef<-nma.obj$RelativeRisk$Median_CI
      trans<-"id"
    }else{
      if(!is.null(nma.obj$LogRelativeRisk)){
        ef<-nma.obj$LogRelativeRisk$Median_CI
        trans<-"exp"
      }else{
        stop("the specified effect size is not estimated in nma.obj.")
      }
    }
  }

  if(effect.size=="LRR"){
    if(!is.null(nma.obj$LogRelativeRisk)){
      ef<-nma.obj$LogRelativeRisk$Median_CI
      trans<-"id"
    }else{
      if(!is.null(nma.obj$RelativeRisk)){
        ef<-nma.obj$RelativeRisk$Median_CI
        trans<-"log"
      }else{
        stop("the specified effect size is not estimated in nma.obj.")
      }
    }
  }

  if(effect.size=="RD"){
    if(!is.null(nma.obj$RiskDifference)){
      ef<-nma.obj$RiskDifference$Median_CI
      trans<-"id"
    }else{
      stop("the specified effect size is not estimated in nma.obj.")
    }
  }

  if(effect.size=="diff"){
    if(!is.null(nma.obj$EffectDiff)){
      ef<-nma.obj$EffectDiff$Median_CI
      trans<-"id"
    }else{
      stop("the specified effect size is not estimated in nma.obj.")
    }
  }

  if(effect.size=="ratio"){
    if(!is.null(nma.obj$RateRatio)){
      ef<-nma.obj$RateRatio$Median_CI
      trans<-"id"
    }else{
      if(!is.null(nma.obj$LogRateRatio)){
        ef<-nma.obj$LogRateRatio$Median_CI
        trans<-"exp"
      }else{
        stop("the specified effect size is not estimated in nma.obj.")
      }
    }
  }

  if(effect.size=="logratio"){
    if(!is.null(nma.obj$LogRateRatio)){
      ef<-nma.obj$LogRateRatio$Median_CI
      trans<-"id"
    }else{
      if(!is.null(nma.obj$RateRatio)){
        ef<-nma.obj$RateRatio$Median_CI
        trans<-"log"
      }else{
        stop("the specified effect size is not estimated in nma.obj.")
      }
    }
  }

  trtname<-rownames(ef)
  ntrt<-length(trtname)
  if(missing(reference)) reference<-trtname[1]
  if(!is.element(reference,trtname)) stop("the reference treatment name is not found.")
  contrast<-ef[,reference]
  contrast<-contrast[contrast!="--"]
  contrast.name<-paste(trtname[trtname!=reference], "vs.", reference)

  med<-low<-upp<-numeric(ntrt-1)
  for(i in 1:(ntrt-1)){
    str<-contrast[i]
    split1<-strsplit(str,split=" \\(")
    med[i]<-as.numeric(split1[[1]][1])
    str2<-split1[[1]][2]
    split2<-strsplit(str2,split=", ")
    low[i]<-as.numeric(split2[[1]][1])
    upp[i]<-as.numeric(gsub("\\)","",split2[[1]][2]))
  }
  if(trans=="log"){
    med<-log(med)
    low<-log(low)
    upp<-log(upp)
  }
  if(trans=="exp"){
    med<-exp(med)
    low<-exp(low)
    upp<-exp(upp)
  }

  if(effect.size=="OR"){
    efname<-"Odds Ratio"
    nullval<-1
  }
  if(effect.size=="LOR"){
    efname<-"Log Odds Ratio"
    nullval<-0
  }
  if(effect.size=="RR"){
    efname<-"Relative Risk"
    nullval<-1
  }
  if(effect.size=="LRR"){
    efname<-"Log Relative Risk"
    nullval<-0
  }
  if(effect.size=="RD"){
    efname<-"Risk Difference"
    nullval<-0
  }
  if(effect.size=="diff"){
    efname<-"Effect Difference"
    nullval<-0
  }
  if(effect.size=="ratio"){
    efname<-"Rate Ratio"
    nullval<-1
  }
  if(effect.size=="logratio"){
    efname<-"Log Rate Ratio"
    nullval<-0
  }

  rg<-max(upp)-min(low)
  incr<-signif(rg/5,1)
  len.out<-ceiling(max(c(abs(max(upp)-nullval),abs(min(low)-nullval)))/incr)+100
  xticks<-c(seq(from=nullval,by=incr,length.out=len.out),seq(from=nullval-incr,by=-incr,length.out=len.out))
  x.left<-min(low)-rg/10
  x.right<-max(upp)+rg/10

  if(missing(height)) height<-ntrt-1
  pdf(paste("ContrastPlot_",effect.size,".pdf",sep=""),width=width,height=height)
  par(mar=c(4.1,max(nchar(contrast.name))*left.margin,1,1))
  plot(c(low,upp),rep(rev(1:(ntrt-1)),2),xlim=c(x.left,x.right),ylim=c(0.8,ntrt-1),frame.plot=FALSE,xaxt="n",yaxt="n",xlab=efname,ylab="",cex=0.1,col="white")
  points(med,rev(1:(ntrt-1)),pch=15)
  axis(1,at=xticks,labels=xticks)
  for(i in 1:(ntrt-1)){
    lines(x=c(low[i],upp[i]),y=c(ntrt-i,ntrt-i))
  }
  abline(v=nullval,col="grey")
  for(i in 1:(ntrt-1)){
    mtext(contrast.name[i],side=2,at=ntrt-i,las=1)
  }
  par(mar=c(5,4,4,2)+0.1)
  garbage<-dev.off()
}