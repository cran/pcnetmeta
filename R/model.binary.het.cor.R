model.binary.het.cor <- function(prior.type="invwishart"){
if(prior.type=="invwishart"){
cat(
"model{
 for(i in 1:len){
  p[i]<-phi(mu[t[i]]+vi[s[i],t[i]])
  r[i]~dbin(p[i],totaln[i])
 }
 for(j in 1:nstudy){
  vi[j,1:ntrt]~dmnorm(zeros[1:ntrt],T[1:ntrt,1:ntrt])
 }
 for(j in 1:ntrt){
  AR[j]<-phi(mu[j]/sqrt(1+invT[j,j]))
  mu[j]~dnorm(0,0.001)
  sigma[j]<-sqrt(invT[j,j])
 }
 invT[1:ntrt,1:ntrt]<-inverse(T[,])
 T[1:ntrt,1:ntrt]~dwish(I[1:ntrt,1:ntrt],ntrt+1)
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

if(!is.element(prior.type,c("invwishart"))){
  stop("specified prior type are wrong.")
}
}