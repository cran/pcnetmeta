model.py.hom.ind <- function(prior.type="unif",rank.prob=TRUE){
if(prior.type=="unif" & rank.prob){
cat(
"model{
 for(i in 1:len){
  y[i]~dpois(py[i]*lambda[i])
  lambda[i]<-exp(mu[t[i]]+sigma*vi[s[i]])
 }
 for(j in 1:nstudy){
  vi[j]~dnorm(0,1)
 }
 for(j in 1:ntrt){
  rate[j]<-exp(mu[j]+pow(sigma,2)/2)
  lograte[j]<-log(rate[j])
  mu[j]~dnorm(0,0.001)
 }
 sigma~dunif(0,c)
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
}",file="tempmodel.txt")
}

if(prior.type=="unif" & !rank.prob){
cat(
"model{
 for(i in 1:len){
  y[i]~dpois(py[i]*lambda[i])
  lambda[i]<-exp(mu[t[i]]+sigma*vi[s[i]])
 }
 for(j in 1:nstudy){
  vi[j]~dnorm(0,1)
 }
 for(j in 1:ntrt){
  rate[j]<-exp(mu[j]+pow(sigma,2)/2)
  lograte[j]<-log(rate[j])
  mu[j]~dnorm(0,0.001)
 }
 sigma~dunif(0,c)
 for(j in 1:ntrt){        
  for(k in 1:ntrt){
   ratio[j,k]<-rate[j]/rate[k]
   logratio[j,k]<-log(ratio[j,k])
  }
 }
}",file="tempmodel.txt")
}

if(prior.type=="invgamma" & rank.prob){
cat(
"model{
 for(i in 1:len){
  y[i]~dpois(py[i]*lambda[i])
  lambda[i]<-exp(mu[t[i]]+sigma*vi[s[i]])
 }
 for(j in 1:nstudy){
  vi[j]~dnorm(0,1)
 }
 for(j in 1:ntrt){
  rate[j]<-exp(mu[j]+pow(sigma,2)/2)
  lograte[j]<-log(rate[j])
  mu[j]~dnorm(0,0.001)
 }
 sigma<-1/sqrt(inv.sig.sq)
 inv.sig.sq~dgamma(a,b)
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
}",file="tempmodel.txt")
}

if(prior.type=="invgamma" & !rank.prob){
cat(
"model{
 for(i in 1:len){
  y[i]~dpois(py[i]*lambda[i])
  lambda[i]<-exp(mu[t[i]]+sigma*vi[s[i]])
 }
 for(j in 1:nstudy){
  vi[j]~dnorm(0,1)
 }
 for(j in 1:ntrt){
  rate[j]<-exp(mu[j]+pow(sigma,2)/2)
  lograte[j]<-log(rate[j])
  mu[j]~dnorm(0,0.001)
 }
 sigma<-1/sqrt(inv.sig.sq)
 inv.sig.sq~dgamma(a,b)
 for(j in 1:ntrt){        
  for(k in 1:ntrt){
   ratio[j,k]<-rate[j]/rate[k]
   logratio[j,k]<-log(ratio[j,k])
  }
 }
}",file="tempmodel.txt")
}

if(!is.element(prior.type,c("unif","invgamma"))){
  stop("specified prior type are wrong.")
}
}