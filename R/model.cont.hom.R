model.cont.hom <- function(prior.type="unif"){
if(prior.type=="unif"){
cat(
"model{
 for(i in 1:len){
  mean[i]~dnorm(theta[i],n[i]/pow(sd[i],2))
  theta[i]<-mu[t[i]]+sigma*vi[s[i]]
 }
 for(j in 1:nstudy){
  vi[j]~dnorm(0,1)
 }
 for(j in 1:ntrt){
  mu[j]~dnorm(0,0.001)
 }
 for(i in 1:ntrt){
  for(j in 1:ntrt){
   diff[i,j]<-mu[i]-mu[j]
  }
 }
 sigma~dunif(0,c)
 rk[1:ntrt]<-(ntrt+1-rank(mu[]))*ifelse(higher.better,1,0)+(rank(mu[]))*ifelse(higher.better,0,1)
 best[1:ntrt]<-equals(rk[],1)
}",file="tempmodel.txt")
}

if(prior.type=="invgamma"){
cat(
"model{
 for(i in 1:len){
  mean[i]~dnorm(theta[i],n[i]/pow(sd[i],2))
  theta[i]<-mu[t[i]]+sigma*vi[s[i]]
 }
 for(j in 1:nstudy){
  vi[j]~dnorm(0,1)
 }
 for(j in 1:ntrt){
  mu[j]~dnorm(0,0.001)
 }
 for(i in 1:ntrt){
  for(j in 1:ntrt){
   diff[i,j]<-mu[i]-mu[j]
  }
 }
 sigma<-1/sqrt(inv.sig.sq)
 inv.sig.sq~dgamma(a,b)
 rk[1:ntrt]<-(ntrt+1-rank(mu[]))*ifelse(higher.better,1,0)+(rank(mu[]))*ifelse(higher.better,0,1)
 best[1:ntrt]<-equals(rk[],1)
}",file="tempmodel.txt")
}

if(!is.element(prior.type,c("unif","invgamma"))){
  stop("specified prior type are wrong.")
}
}