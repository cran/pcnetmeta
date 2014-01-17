nma.ab <-
function(s.id,t.id,event.n,total.n,o.path=getwd(),f.name="",hom=FALSE,R,param=c("probt","RR","RD","OR","rk","best"),ity="estimate",n.iter=200000,n.burnin=floor(n.iter/2),n.chains=3,n.thin=max(1,floor((n.iter-n.burnin)/100000)),dic=TRUE,trace=FALSE,postdens=FALSE){
  ## check the input parameters
  options(warn=1)
  if(missing(s.id)) stop("need to specify study id.")
  if(missing(t.id)) stop("need to specify treatment.")
  if(missing(event.n)) stop("need to specify event number.")
  if(missing(total.n)) stop("need to specify total number.")
  if(ity!="same" & ity!="estimate") stop("initial value type (ity) must be set as \"same\" or \"estimate\".")
  if(length(s.id)!=length(t.id) | length(t.id)!=length(event.n) | length(event.n)!=length(total.n) | length(total.n)!=length(s.id)){
    stop("the data input do not have the same length.")
  }
  if(min(event.n)<0 | min(total.n)<0) stop("event number and total number must be greater than zero.")
  if(!all(event.n<=total.n)) stop("total number must be greater than event number.")
  if(!all(total.n>0)) stop("total number must be positive.")
  if(!all(event.n>=0)) stop("event number must be non-negative.")
  if(!all(event.n%%1==0) | !all(total.n%%1==0)) warning("at least one event number or total number is not integer.")
  tN<-max(t.id) # total number of treatments
  if(!all(t.id %in% 1:tN) | !all(1:tN %in% t.id)) stop("treatments must be labeled as consecutive integers starting from 1 to length(t.id).")

  ## set variables
  tS<-max(s.id) # total number of studies
  sN<-length(s.id)
  if(hom==FALSE) if(missing(R)) R<-matrix(0.005,ncol=tN,nrow=tN)+diag(1-0.005,ncol=tN,nrow=tN)

  ## jags codes
  setwd(o.path)
  cat(ifelse(hom,
"model{
 for(i in 1:sN){
  p[i]<-phi(mu[t[i]]+sigma*vi[s[i]])
  r[i]~dbin(p[i],totaln[i])
 }
 for(j in 1:tS){
  vi[j]~dnorm(0,1)
 }
 sigma<-1/sqrt(tau)
 tau~dgamma(1,1)
 for(j in 1:tN){
  mu[j]~dnorm(0,0.001)
  probt[j]<-phi(mu[j]/sqrt(1+1/tau))
 }
 for(j in 1:tN){
  for(k in 1:tN){
   RR[j,k]<-probt[j]/probt[k]
   RD[j,k]<-probt[j]-probt[k]
   OR[j,k]<-probt[j]/(1-probt[j])/probt[k]*(1-probt[k])
  }
 }
 rk[1:tN]<-tN+1-rank(probt[])
 best[1:tN]<-equals(rk[],1)
}",
"model{
 for(i in 1:sN){
  p[i]<-phi(mu[t[i]]+vi[s[i],t[i]])
  r[i]~dbin(p[i],totaln[i])
 }
 for(j in 1:tS){
  vi[j,1:tN]~dmnorm(mn[1:tN],T[1:tN,1:tN])
 }
 invT[1:tN,1:tN]<-inverse(T[,])
 for(j in 1:tN){
  mu[j]~dnorm(0,0.001)
  sigma[j]<-sqrt(invT[j,j])
  probt[j]<-phi(mu[j]/sqrt(1+invT[j,j]))
 }
 T[1:tN,1:tN]~dwish(R[1:tN,1:tN],tN)
 for(k in 1:tN){        
  for(j in 1:tN){
   RR[k,j]<-probt[k]/probt[j]
   RD[k,j]<-probt[k]-probt[j]
   OR[k,j]<-probt[k]/(1-probt[k])/probt[j]*(1-probt[j])
  }
 }
 rk[1:tN]<-tN+1-rank(probt[])
 best[1:tN]<-equals(rk[],1)
}"), file="jags-model.txt")

  ## run jags
  if(hom){
    data.jags<-list(s=s.id,t=t.id,r=event.n,totaln=total.n,sN=sN,tS=tS,tN=tN)}else{
      data.jags<-list(s=s.id,t=t.id,r=event.n,totaln=total.n,sN=sN,tS=tS,tN=tN,mn=rep(0,tN),R=R)}
  if(ity=="same"){
    if(hom) {init.jags<-function(){list(mu=rep(0,tN),tau=1,vi=rep(0,tS))}}else{
      init.jags<-function(){list(mu=rep(0,tN))}}}else{
      init<-NULL
      for(trt in 1:tN){
        init<-c(init,sum(event.n[t.id==trt])/sum(total.n[t.id==trt]))
      }
      if(hom) {init.jags<-function(){list(mu=qnorm(init),tau=1,vi=rep(0,tS))}}else{
        init.jags<-function(){list(mu=qnorm(init))}}
    }

  jags.result<-jags(data=data.jags,inits=init.jags,parameters.to.save=param,model.file="jags-model.txt",n.chains=n.chains,n.iter=n.iter,n.burnin=n.burnin,n.thin=n.thin,DIC=dic)
  gc()
  summary.stat<-signif(jags.result$BUGSoutput$summary[,1:7],4)
  gc()
  rnames<-rownames(summary.stat)
  delete<-NULL
  for(k in 1:tN){
   delete<-c(delete,paste("OR","[",as.character(k),",",as.character(k),"]",sep=""),
                    paste("RD","[",as.character(k),",",as.character(k),"]",sep=""),
                    paste("RR","[",as.character(k),",",as.character(k),"]",sep=""))
  }
  compare<-rnames %in% delete
  smry<-summary.stat[which(compare==FALSE),]
  setwd(o.path)
  write.table(smry,file=file.path(getwd(),paste(f.name,"Summary.stat",sep="")))
  gc()
  if(dic){
    dic.stat<-rbind(jags.result$BUGSoutput$pD,jags.result$BUGSoutput$DIC)
    rownames(dic.stat)<-c("pD","DIC")
    write.table(dic.stat,file=file.path(getwd(),paste(f.name,"DIC.stat",sep="")),col.names=FALSE)
  }
  gc()
  if(trace|postdens) {pn<-paste("probt[",1:tN,"]",sep="");mcmc<-as.mcmc(jags.result)}
  rm(jags.result)
  gc()
  if(trace){
    pdf(paste(f.name,"TracePlot.pdf",sep=""))
    for(i in 1:tN){
      for(j in 1:n.chains){
        tempmcmc<-mcmc[,pn[i]][[j]]
        plot(1:length(tempmcmc),tempmcmc,type="l",col="red",xlab="Iterations",ylab=pn[i],main=paste("Trace plot for chain",j,"of",pn[i]))
      }
    }
    dev.off()
  }
  gc()
  if(postdens){
    cn<-paste("chain",1:n.chains)
    pdf(paste(f.name,"DensPlot.pdf",sep=""))
    for(i in 1:tN){
      maxs<-NULL
      for(j in 1:n.chains){
        maxs<-c(maxs,max(density(mcmc[,pn[i]][[j]])$y))
        tempi<-which(maxs==max(maxs))
      }
      plot(density(mcmc[,pn[i]][[tempi]]),xlab="Event rate",ylab="Density",main=paste("Density plot for",pn[i]),type="l",col=tempi,lwd=2)
      for(j in 1:n.chains){
        if(j!=tempi) lines(density(mcmc[,pn[i]][[j]]),type="l",col=j,lwd=2)
      }
    }
    dev.off()
  }
  unlink("jags-model.txt")
  options(warn=0)
}
