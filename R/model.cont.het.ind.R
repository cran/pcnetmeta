model.cont.het.ind <- function(prior.type="unif"){
if(prior.type=="unif"){
cat(
"model{
 for(i in 1:len){
  mean[i]~dnorm(theta[i],n[i]/pow(sd[i],2))
  theta[i]<-mu[t[i]]+sigma[t[i]]*vi[s[i]]
 }
 for(j in 1:nstudy){
  vi[j]~dnorm(0,1)
 }
 for(j in 1:ntrt){
  mu[j]~dnorm(0,0.001)
  sigma[j]~dunif(0,c)
 }
 for(i in 1:ntrt){
  for(j in 1:ntrt){
   diff[i,j]<-mu[i]-mu[j]
  }
 }
 rk[1:ntrt]<-(ntrt+1-rank(mu[]))*ifelse(higher.better,1,0)+(rank(mu[]))*ifelse(higher.better,0,1)
 best[1:ntrt]<-equals(rk[],1)
}",file="tempmodel.txt")
}

if(prior.type=="invgamma"){
cat(
"model{
 for(i in 1:len){
  mean[i]~dnorm(theta[i],n[i]/pow(sd[i],2))
  theta[i]<-mu[t[i]]+sigma[t[i]]*vi[s[i]]
 }
 for(j in 1:nstudy){
  vi[j]~dnorm(0,1)
 }
 for(j in 1:ntrt){
  mu[j]~dnorm(0,0.001)
  sigma[j]<-1/sqrt(inv.sig.sq[j])
  inv.sig.sq[j]~dgamma(a,b)
 }
 for(i in 1:ntrt){
  for(j in 1:ntrt){
   diff[i,j]<-mu[i]-mu[j]
  }
 }
 rk[1:ntrt]<-(ntrt+1-rank(mu[]))*ifelse(higher.better,1,0)+(rank(mu[]))*ifelse(higher.better,0,1)
 best[1:ntrt]<-equals(rk[],1)
}",file="tempmodel.txt")
}

if(!is.element(prior.type,c("unif","invgamma"))){
  stop("specified prior type are wrong.")
}
}