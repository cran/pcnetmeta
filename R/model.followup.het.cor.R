model.followup.het.cor <- function(prior.type="invwishart",rank.prob=TRUE){
if(prior.type=="invwishart" & rank.prob){
modelstring<-"
model{
 for(i in 1:len){
  y[i]~dbin(p[i],totaln[i])
  cloglog(p[i])<-log(f[i])+log(lambda[i])
  log(lambda[i])<-mu[t[i]]+vi[s[i],t[i]]
 }
 for(j in 1:nstudy){
  vi[j,1:ntrt]~dmnorm(zeros[1:ntrt],T[1:ntrt,1:ntrt])
 }
 for(j in 1:ntrt){
  rate[j]<-exp(mu[j]+invT[j,j]/2)
  lograte[j]<-log(rate[j])
  mu[j]~dnorm(0,0.001)
  sigma[j]<-sqrt(invT[j,j])
 }
 invT[1:ntrt,1:ntrt]<-inverse(T[,])
 T[1:ntrt,1:ntrt]~dwish(I[1:ntrt,1:ntrt],ntrt+1)
 for(j in 1:ntrt){        
  for(k in 1:ntrt){
   ratio[j,k]<-rate[j]/rate[k]
   logratio[j,k]<-log(ratio[j,k])
  }
 }
 rk[1:ntrt]<-(ntrt+1-rank(rate[]))*ifelse(higher.better,1,0)+(rank(rate[]))*ifelse(higher.better,0,1)
 for(i in 1:ntrt){
  rank.prob[1:ntrt,i]<-equals(rk[],i)
 }
}
"
}

if(prior.type=="invwishart" & !rank.prob){
modelstring<-"
model{
 for(i in 1:len){
  y[i]~dbin(p[i],totaln[i])
  cloglog(p[i])<-log(f[i])+log(lambda[i])
  log(lambda[i])<-mu[t[i]]+vi[s[i],t[i]]
 }
 for(j in 1:nstudy){
  vi[j,1:ntrt]~dmnorm(zeros[1:ntrt],T[1:ntrt,1:ntrt])
 }
 for(j in 1:ntrt){
  rate[j]<-exp(mu[j]+invT[j,j]/2)
  lograte[j]<-log(rate[j])
  mu[j]~dnorm(0,0.001)
  sigma[j]<-sqrt(invT[j,j])
 }
 invT[1:ntrt,1:ntrt]<-inverse(T[,])
 T[1:ntrt,1:ntrt]~dwish(I[1:ntrt,1:ntrt],ntrt+1)
 for(j in 1:ntrt){        
  for(k in 1:ntrt){
   ratio[j,k]<-rate[j]/rate[k]
   logratio[j,k]<-log(ratio[j,k])
  }
 }
}
"
}

if(!is.element(prior.type,c("invwishart"))){
  stop("specified prior type is wrong.")
}

return(modelstring)
}