model.binary.hom <- function(prior.type="unif"){
if(prior.type=="unif"){
cat(
"model{
 for(i in 1:len){
  p[i]<-phi(mu[t[i]]+sigma*vi[s[i]])
  r[i]~dbin(p[i],totaln[i])
 }
 for(j in 1:nstudy){
  vi[j]~dnorm(0,1)
 }
 for(j in 1:ntrt){
  mu[j]~dnorm(0,0.001)
  AR[j]<-phi(mu[j]/sqrt(1+pow(sigma,2)))
 }
 sigma~dunif(0,c)
 for(j in 1:ntrt){
  for(k in 1:ntrt){
   LRR[j,k]<-log(RR[j,k])
   LOR[j,k]<-log(OR[j,k])
   RR[j,k]<-AR[j]/AR[k]
   RD[j,k]<-AR[j]-AR[k]
   OR[j,k]<-AR[j]/(1-AR[j])/AR[k]*(1-AR[k])
  }
 }
 rk[1:ntrt]<-(ntrt+1-rank(AR[]))*ifelse(higher.better,1,0)+(rank(AR[]))*ifelse(higher.better,0,1)
 best[1:ntrt]<-equals(rk[],1)
}",file="tempmodel.txt")
}

if(prior.type=="invgamma"){
cat(
"model{
 for(i in 1:len){
  p[i]<-phi(mu[t[i]]+sigma*vi[s[i]])
  r[i]~dbin(p[i],totaln[i])
 }
 for(j in 1:nstudy){
  vi[j]~dnorm(0,1)
 }
 for(j in 1:ntrt){
  mu[j]~dnorm(0,0.001)
  AR[j]<-phi(mu[j]/sqrt(1+1/inv.sig.sq))
 }
 sigma<-1/sqrt(inv.sig.sq)
 inv.sig.sq~dgamma(a,b)
 for(j in 1:ntrt){
  for(k in 1:ntrt){
   LRR[j,k]<-log(RR[j,k])
   LOR[j,k]<-log(OR[j,k])
   RR[j,k]<-AR[j]/AR[k]
   RD[j,k]<-AR[j]-AR[k]
   OR[j,k]<-AR[j]/(1-AR[j])/AR[k]*(1-AR[k])
  }
 }
 rk[1:ntrt]<-(ntrt+1-rank(AR[]))*ifelse(higher.better,1,0)+(rank(AR[]))*ifelse(higher.better,0,1)
 best[1:ntrt]<-equals(rk[],1)
}",file="tempmodel.txt")
}

if(!is.element(prior.type,c("unif","invgamma"))){
  stop("specified prior type are wrong.")
}
}