model.followup.hom <- function(prior.type="unif"){
if(prior.type=="unif"){
cat(
"model{
 for(i in 1:len){
  y[i]~dbin(p[i],totaln[i])
  cloglog(p[i])<-log(f[i])+log(lambda[i])
  log(lambda[i])<-mu[t[i]]+sigma*vi[s[i]]
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
 best[1:ntrt]<-equals(rk[],1)
}",file="tempmodel.txt")
}

if(prior.type=="invgamma"){
cat(
"model{
 for(i in 1:len){
  y[i]~dbin(p[i],totaln[i])
  cloglog(p[i])<-log(f[i])+log(lambda[i])
  log(lambda[i])<-mu[t[i]]+sigma*vi[s[i]]
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
 best[1:ntrt]<-equals(rk[],1)
}",file="tempmodel.txt")
}

if(!is.element(prior.type,c("unif","invgamma"))){
  stop("specified prior type are wrong.")
}
}