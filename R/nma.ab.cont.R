nma.ab.cont <-
function(s.id,t.id,mean,sd,total.n,trtname,model="het_cor",prior.type,a=0.001,b=0.001,c=10,param=c("mu","diff","best"),higher.better=FALSE,digits=4,n.adapt=5000,n.iter=100000,n.burnin=floor(n.iter/2),n.chains=3,n.thin=max(1,floor((n.iter-n.burnin)/100000)),conv.diag=FALSE,trace="",dic=FALSE,postdens=FALSE){
  ## check the input parameters
  options(warn=1)
  if(missing(s.id)) stop("need to specify study id.")
  if(missing(t.id)) stop("need to specify treatment.")
  if(missing(mean) | missing(sd)) stop("need to specify mean and sd of the continuous outcomes.")
  if(missing(total.n)) stop("need to specify total number.")
  if(length(s.id)!=length(t.id) | length(t.id)!=length(mean) | length(mean)!=length(sd) | length(sd)!=length(total.n) | length(total.n)!=length(s.id)){
    stop("the data input do not have the same length.")
  }
  if(!all(total.n>0)) stop("total number must be positive.")
  if(!all(total.n%%1==0)) warning("at least one event number or total number is not integer.")
  if(!is.element(model,c("hom","het_ind","het_cor"))) stop("model should be specified as \"hom\", \"het_ind\", or \"het_cor\".")

  ## make ids continuous
  s.id.o<-s.id
  t.id.o<-t.id
  s.label<-sort(unique(s.id.o))
  t.label<-sort(unique(t.id.o))
  nstudy<-length(s.label) # total number of studies
  ntrt<-length(t.label) # total number of treatments
  len<-length(s.id)
  s.id<-numeric(nstudy)
  for(i in 1:nstudy) {s.id[which(s.id.o==s.label[i])]<-i}
  t.id<-numeric(ntrt)
  for(i in 1:ntrt) {t.id[which(t.id.o==t.label[i])]<-i}

  if(missing(trtname)){
    if(is.numeric(t.id.o)) {trtname<-paste("Trt",t.label,sep="")}
    if(is.character(t.id.o)) {trtname<-t.label}
  }
  if(length(trtname)!=length(unique(t.id))) stop("the number of treatment names does not match for specified treatment id.")
  if(missing(prior.type)) prior.type<-ifelse(model=="het_cor","invwishart","unif")

  ## jags model
  if(model=="hom"){
    model.cont.hom(prior.type)
  }
  if(model=="het_ind"){
    model.cont.het.ind(prior.type)
  }
  if(model=="het_cor"){
    I <- diag(ntrt)
    model.cont.het.cor(prior.type)
  }

  ## jags data
  if(model == "hom"| model == "het_ind"){
    if(prior.type == "unif"){
      data.jags<-list(s=s.id,t=t.id,mean=mean,sd=sd,n=total.n,len=len,nstudy=nstudy,ntrt=ntrt,c=c,higher.better=higher.better)
    }
    if(prior.type == "invgamma"){
      data.jags<-list(s=s.id,t=t.id,mean=mean,sd=sd,n=total.n,len=len,nstudy=nstudy,ntrt=ntrt,a=a,b=b,higher.better=higher.better)
    }
  }
  if(model=="het_cor"){
    data.jags<-list(s=s.id,t=t.id,mean=mean,sd=sd,n=total.n,len=len,nstudy=nstudy,ntrt=ntrt,zeros=rep(0,ntrt),I=I,higher.better=higher.better)
  }

  ## jags initial value
  rng.seeds<-sample(1000000,n.chains)
  mu.init<-numeric(ntrt)
  for(i in 1:ntrt){
    mu.init[i]<-sum(mean[t.id==t.id[i]]*total.n[t.id==t.id[i]])/sum(total.n[t.id==t.id[i]])
  }
  init.jags<-list(NULL)
  if(model=="hom"){
    if(prior.type=="unif"){
      for(ii in 1:n.chains){
        init.jags[[ii]]<-list(mu=mu.init,vi=rep(0,nstudy),sigma=c/2,.RNG.name="base::Wichmann-Hill",.RNG.seed=rng.seeds[ii])
      }
    }
    if(prior.type=="invgamma"){
      for(ii in 1:n.chains){
        init.jags[[ii]]<-list(mu=mu.init,vi=rep(0,nstudy),inv.sig.sq=a/b,.RNG.name="base::Wichmann-Hill",.RNG.seed=rng.seeds[ii])
      }
    }
  }
  if(model=="het_ind"){
    if(prior.type=="unif"){
      for(ii in 1:n.chains){
        init.jags[[ii]]<-list(mu=mu.init,vi=rep(0,nstudy),sigma=rep(c/2,ntrt),.RNG.name="base::Wichmann-Hill",.RNG.seed=rng.seeds[ii])
      }
    }
    if(prior.type=="invgamma"){
      for(ii in 1:n.chains){
        init.jags[[ii]]<-list(mu=mu.init,vi=rep(0,nstudy),inv.sig.sq=rep(a/b,ntrt),.RNG.name="base::Wichmann-Hill",.RNG.seed=rng.seeds[ii])
      }
    }
  }
  if(model=="het_cor"){
    for(ii in 1:n.chains){
      init.jags[[ii]]<-list(mu=mu.init,vi=matrix(0,nstudy,ntrt),T=(ntrt+1)*I,.RNG.name="base::Wichmann-Hill",.RNG.seed=rng.seeds[ii])
    }
  }

  ## parameters to be monitored in jags
  if(!is.element("mu",param)) param<-c("mu",param)
  monitor<-param[!is.element(param,c("diff"))]
  if(is.element("diff",param)){
    for(ii in 1:ntrt){
      for(jj in 1:ntrt){
        if(ii < jj) monitor<-c(monitor,paste("diff[",ii,",",jj,"]",sep=""))
      }
    }
  }

  ## run jags
  cat("start running MCMC...\n")
  jags.m<-jags.model(file="tempmodel.txt",data=data.jags,inits=init.jags,n.chains=n.chains,n.adapt=n.adapt)
  unlink("tempmodel.txt")
  update(jags.m,n.iter=n.burnin)
  jags.out<-coda.samples(model=jags.m,variable.names=monitor,n.iter=n.iter,thin=n.thin)
  smry<-summary(jags.out)
  smry<-cbind(smry$statistics[,c("Mean","SD")],smry$quantiles[,c("2.5%","50%","97.5%")])
  smry<-signif(smry,digits=digits)

  out<-NULL
  mu.id<-grep("mu",rownames(smry))
  mu.stat<-array(paste(format(smry[mu.id,"Mean"],digits=digits)," (",format(smry[mu.id,"SD"],digits=digits),")",sep=""),dim=c(ntrt,1))
  colnames(mu.stat)<-"Mean (SD)"
  rownames(mu.stat)<-trtname
  mu.quan<-array(paste(format(smry[mu.id,"50%"],digits=digits)," (",format(smry[mu.id,"2.5%"],digits=digits),", ",format(smry[mu.id,"97.5%"],digits=digits),")",sep=""),dim=c(ntrt,1))
  colnames(mu.quan)<-"Median (95% CI)"
  rownames(mu.quan)<-trtname
  out$TrtEffect<-list(Mean_SD=noquote(mu.stat),Median_CI=noquote(mu.quan))

  if(is.element("diff",param)){
    diff.stat<-diff.quan<-array(NA,dim=c(ntrt,ntrt))
    colnames(diff.stat)<-colnames(diff.quan)<-rownames(diff.stat)<-rownames(diff.quan)<-trtname
    for(i in 1:ntrt){
      diff.stat[i,i]<-diff.quan[i,i]<-"--"
      for(j in 1:ntrt){
        if(i < j){
          diff.ij<-paste("diff[",i,",",j,"]",sep="")
          diff.stat[i,j]<-paste(format(smry[diff.ij,"Mean"],digits=digits,nsmall=digits)," (",format(smry[diff.ij,"SD"],digits=digits,nsmall=digits),")",sep="")
          diff.stat[j,i]<-paste(format(-smry[diff.ij,"Mean"],digits=digits,nsmall=digits)," (",format(smry[diff.ij,"SD"],digits=digits,nsmall=digits),")",sep="")
          diff.quan[i,j]<-paste(format(smry[diff.ij,"50%"],digits=digits,nsmall=digits)," (",format(smry[diff.ij,"2.5%"],digits=digits,nsmall=digits),", ",format(smry[diff.ij,"97.5%"],digits=digits,nsmall=digits),")",sep="")
          diff.quan[j,i]<-paste(format(-smry[diff.ij,"50%"],digits=digits,nsmall=digits)," (",format(-smry[diff.ij,"97.5%"],digits=digits,nsmall=digits),", ",format(-smry[diff.ij,"2.5%"],digits=digits,nsmall=digits),")",sep="")
        }
      }
    }
    out$EffectDiff<-list(Mean_SD=noquote(diff.stat),Median_CI=noquote(diff.quan))
  }

  if(is.element("best",param)){
    best.id<-grep("best",rownames(smry))
    best.stat<-array(format(round(smry[best.id,"Mean"],digits=digits),digits=digits),dim=c(ntrt,1))
    colnames(best.stat)<-""
    rownames(best.stat)<-trtname
    out$ProbOfBestTrt<-noquote(best.stat)
  }

  if(conv.diag){
    cat("start calculating MCMC convergence diagnostic statistics...\n")
    conv.out<-gelman.diag(jags.out,multivariate=FALSE)
    conv.out<-conv.out$psrf
    if(is.element("best",param)){
      best.id<-grep("best",rownames(conv.out))
      conv.out<-conv.out[-best.id,]
    }
    write.table(conv.out,"ConvergenceDiagnostic.txt",row.names=rownames(conv.out),col.names=TRUE)
  }

  if(dic){
    cat("start calculating deviance information criterion statistics...\n")
    dic.out<-dic.samples(model=jags.m,n.iter=n.iter,thin=n.thin)
    dev<-sum(dic.out$deviance)
    pen<-sum(dic.out$penalty)
    pen.dev<-dev+pen
    dic.stat<-rbind(dev,pen,pen.dev)
    rownames(dic.stat)<-c("Mean deviance","Penalty","Penalized deviance")
    colnames(dic.stat)<-""
    out$DIC<-dic.stat
  }

  if(!all(trace=="")){
    cat("start saving trace plots...\n")
  }

  if(is.element("mu",trace)){
    for(i in 1:ntrt){
      png(paste("TracePlot_mu_",trtname[i],".png",sep=""),res=600,height=8.5,width=11,units="in")
      par(mfrow=c(n.chains,1))
      for(j in 1:n.chains){
        temp<-as.vector(jags.out[[j]][,paste("mu[",i,"]",sep="")])
        plot(temp,type="l",col="red",ylab="Treatment Effect",xlab="Iterations",main=paste("Chain",j))
      }
      dev.off()
    }
  }
  if(is.element("diff",trace)){
    for(i in 1:ntrt){
      for(k in 1:ntrt){
        if(i < k){
          png(paste("TracePlot_diff_",trtname[i],"_",trtname[k],".png",sep=""),res=600,height=8.5,width=11,units="in")
          par(mfrow=c(n.chains,1))
          for(j in 1:n.chains){
            temp<-as.vector(jags.out[[j]][,paste("diff[",i,",",k,"]",sep="")])
            plot(temp,type="l",col="red",ylab="Effect Difference",xlab="Iterations",main=paste("Chain",j))
          }
          dev.off()
        }
      }
    }
  }

  if(postdens){
    cat("start saving posterior density plot for treatment effect...\n")
    mcmc<-NULL
    dens<-matrix(0,ntrt,3)
    colnames(dens)<-c("ymax","xmin","xmax")
    for(i in 1:ntrt){
      temp<-NULL
      for(j in 1:n.chains){
        temp<-c(temp,as.vector(jags.out[[j]][,paste("mu[",i,"]",sep="")]))
      }
      mcmc[[i]]<-temp
      tempdens<-density(temp)
      dens[i,]<-c(max(tempdens$y),quantile(temp,0.001),quantile(temp,0.999))
    }
    ymax<-max(dens[,"ymax"])
    xmin<-min(dens[,"xmin"])
    xmax<-max(dens[,"xmax"])
    cols<-rainbow(ntrt,s=1,v=0.6)
    pdf("TreatmentEffectDensityPlot.pdf")
    par(mfrow=c(1,1))
    plot(density(mcmc[[1]]),xlim=c(xmin,xmax),ylim=c(0,ymax),xlab="Treatment Effect",ylab="Density",main="",col=cols[1],lty=1,lwd=2)
    for(i in 2:ntrt){
      lines(density(mcmc[[i]]),col=cols[i],lty=i,lwd=2)
    }
    legend("topright",legend=trtname,col=cols,lty=1:ntrt)
    dev.off()
  }
  return(out)
  options(warn=0)
}