nma.ab <-
function(s.id,t.id,event.n,total.n,o.path=getwd(),f.name="",model="het2",sigma2.a=0.001,sigma2.b=0.001,mu.prec=0.001,R,param=c("probt","RR","RD","OR","rk","best"),ity="estimate",n.iter=200000,n.burnin=floor(n.iter/2),n.chains=3,n.thin=max(1,floor((n.iter-n.burnin)/100000)),dic=TRUE,trace=FALSE,postdens=FALSE){
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
  if(!is.element(model,c("hom","het1","het2"))) stop("model should be specified as \"hom\", \"het1\", or \"het2\".")

  ## make ids consecutive
  s.id.o<-s.id
  t.id.o<-t.id
  s.label<-sort(unique(s.id.o))
  t.label<-sort(unique(t.id.o))
  tS<-length(s.label) # total number of studies
  tN<-length(t.label) # total number of treatments
  sN<-length(s.id)
  s.id<-numeric(sN)
  for(i in 1:tS) {s.id[which(s.id.o==s.label[i])]<-i}
  t.id<-numeric(sN)
  for(i in 1:tN) {t.id[which(t.id.o==t.label[i])]<-i}

  ## jags codes
  setwd(o.path)
  if(model=="hom") cat(
"model{
 for(i in 1:sN){
  p[i]<-phi(mu[t[i]]+sigma*vi[s[i]])
  r[i]~dbin(p[i],totaln[i])
 }
 for(j in 1:tS){
  vi[j]~dnorm(0,1)
 }
 sigma<-1/sqrt(tau)
 tau~dgamma(sigma2.a,sigma2.b)
 for(j in 1:tN){
  mu[j]~dnorm(0,mu.prec)
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
}",file="jags-model.txt")

  if(model=="het1") cat(
"model{
 for(i in 1:sN){
  p[i]<-phi(mu[t[i]]+sigma[t[i]]*vi[s[i]])
  r[i]~dbin(p[i],totaln[i])
 }
 for(j in 1:tS){
  vi[j]~dnorm(0,1)
 }
 for(j in 1:tN){
  mu[j]~dnorm(0,mu.prec)
  tau[j]~dgamma(sigma2.a,sigma2.b)
  sigma[j]<-1/sqrt(tau[j])
  probt[j]<-phi(mu[j]/sqrt(1+1/tau[j]))
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
}",file="jags-model.txt")

  if(model=="het2"){
    cat(
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
  mu[j]~dnorm(0,mu.prec)
  sigma[j]<-sqrt(invT[j,j])
  probt[j]<-phi(mu[j]/sqrt(1+invT[j,j]))
 }
 T[1:tN,1:tN]~dwish(R[1:tN,1:tN],tN)
 for(j in 1:tN){        
  for(k in 1:tN){
   RR[j,k]<-probt[j]/probt[k]
   RD[j,k]<-probt[j]-probt[k]
   OR[j,k]<-probt[j]/(1-probt[j])/probt[k]*(1-probt[k])
  }
 }
 rk[1:tN]<-tN+1-rank(probt[])
 best[1:tN]<-equals(rk[],1)
}",file="jags-model.txt")
    if(missing(R)) R<-matrix(0.005,ncol=tN,nrow=tN)+diag(1-0.005,ncol=tN,nrow=tN)
  }

  ## run jags
  if(is.element(model,c("hom","het1"))) data.jags<-list(s=s.id,t=t.id,r=event.n,totaln=total.n,sN=sN,tS=tS,tN=tN,sigma2.a=sigma2.a,sigma2.b=sigma2.b,mu.prec=mu.prec)
  if(model=="het2") data.jags<-list(s=s.id,t=t.id,r=event.n,totaln=total.n,sN=sN,tS=tS,tN=tN,mn=rep(0,tN),R=R,mu.prec=mu.prec)
  if(ity=="same") init.jags<-function(){list(mu=rep(0,tN))} else{
      init<-NULL
      for(trt in 1:tN){
        init<-c(init,sum(event.n[t.id==trt])/sum(total.n[t.id==trt]))
      }
      init.jags<-function(){list(mu=qnorm(init))}
    }

  jags.result<-jags(data=data.jags,inits=init.jags,parameters.to.save=param,model.file="jags-model.txt",n.chains=n.chains,n.iter=n.iter,n.burnin=n.burnin,n.thin=n.thin,DIC=dic)
  gc()
  summary.stat<-signif(jags.result$BUGSoutput$summary[,1:7],4)
  gc()
  rnames<-rownames(summary.stat)
  delete<-NULL
  for(k in 1:tN){
   delete<-c(delete,paste("OR","[",k,",",k,"]",sep=""),
                    paste("RD","[",k,",",k,"]",sep=""),
                    paste("RR","[",k,",",k,"]",sep=""))
  }
  compare<-rnames %in% delete
  smry<-summary.stat[which(compare==FALSE),]
  temp.rnames<-rnames<-rownames(smry)
  for(i in 1:tN){
    temp.rep<-grep(paste("[",i,"]",sep=""),rnames,fixed=TRUE)
    temp.rnames[temp.rep]<-gsub(as.character(i),as.character(t.label[i]),temp.rnames[temp.rep],fixed=TRUE)
    for(j in 1:tN){
      temp.rep<-grep(paste("[",i,",",j,"]",sep=""),rnames,fixed=TRUE)
      temp.rnames[temp.rep]<-gsub(paste("[",i,",",j,"]",sep=""),paste("[",t.label[i],",",t.label[j],"]",sep=""),temp.rnames[temp.rep],fixed=TRUE)
    }
  }
  rownames(smry)<-temp.rnames
  setwd(o.path)
  write.table(smry,file=file.path(getwd(),paste(f.name,"Summary.stat",sep="")))
  gc()
  if(dic){
    dic.stat<-rbind(jags.result$BUGSoutput$pD,jags.result$BUGSoutput$DIC)
    dic.stat<-round(dic.stat,4)
    rownames(dic.stat)<-c("pD","DIC")
    write.table(dic.stat,file=file.path(getwd(),paste(f.name,"DIC.stat",sep="")),col.names=FALSE)
  }
  gc()
  if(trace){
    pdf(paste(f.name,"TracePlot.pdf",sep=""),width=15,height=5)
    traceplot(jags.result,ask=FALSE,varname="probt")
    dev.off()
  }
  gc()
  if(postdens){
    mcmc<-combine.mcmc(as.mcmc(jags.result))
    mcmc<-mcmc[,grep("probt",colnames(mcmc))]
    mcmc<-data.frame(dens=c(mcmc),lines=rep(paste("probt[",1:tN,"]",sep=""),each=dim(mcmc)[1]))
    pdf(paste(f.name,"DensPlot.pdf",sep=""))
    print(densityplot(~dens,data=mcmc,groups=lines,plot.points=FALSE,ref=TRUE,auto.key=list(space="right"),main="Density plot for population-averaged event rate",xlab="Range",ylab="Density"))
    dev.off()
  }
  unlink("jags-model.txt")
  options(warn=0)
}
