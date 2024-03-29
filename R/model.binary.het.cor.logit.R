model.binary.het.cor.logit <- function(prior.type = "chol", rank.prob = TRUE){
if(prior.type == "invwishart" & rank.prob){
modelstring<-"
model{
 for(i in 1:len){
  logit(p[i]) <- mu[t[i]] + vi[s[i], t[i]]
  r[i] ~ dbin(p[i], totaln[i])
  rhat[i] <- p[i]*totaln[i]
  dev[i] <- 2*(r[i]*(log(r[i]) - log(rhat[i])) +
    (totaln[i] - r[i])*(log(totaln[i] - r[i]) - log(totaln[i] - rhat[i])))
 }
 totresdev <- sum(dev[])
 for(j in 1:nstudy){
  vi[j, 1:ntrt] ~ dmnorm(zeros[1:ntrt], T[1:ntrt, 1:ntrt])
 }
 for(j in 1:ntrt){
  AR[j] <- 1/(1 + exp(-mu[j]/sqrt(1 + (16*sqrt(3)/(15*3.1415926))^2*invT[j,j])))
  mu[j] ~ dnorm(0,0.001)
  sigma[j] <- sqrt(invT[j,j])
 }
 invT[1:ntrt, 1:ntrt] <- inverse(T[,])
 T[1:ntrt, 1:ntrt] ~ dwish(I[1:ntrt, 1:ntrt], ntrt + 1)
 for(j in 1:ntrt){        
  for(k in 1:ntrt){
   LRR[j,k] <- log(RR[j,k])
   LOR[j,k] <- log(OR[j,k])
   LOR.med[j,k] <- mu[j] - mu[k]
   RR[j,k] <- AR[j]/AR[k]
   RD[j,k] <- AR[j] - AR[k]
   OR[j,k] <- AR[j]/(1 - AR[j])/AR[k]*(1 - AR[k])
   OR.med[j,k] <- exp(LOR.med[j,k])
  }
 }
 rk[1:ntrt] <- (ntrt + 1 - rank(AR[]))*ifelse(higher.better, 1, 0) + (rank(AR[]))*ifelse(higher.better, 0, 1)
 rk.med[1:ntrt] <- (ntrt + 1 - rank(mu[]))*ifelse(higher.better, 1, 0) + (rank(mu[]))*ifelse(higher.better, 0, 1)
 for(i in 1:ntrt){
  rank.prob[1:ntrt, i] <- equals(rk[], i)
  rank.prob.med[1:ntrt, i] <- equals(rk.med[], i)
 }
}
"
}

if(prior.type == "invwishart" & !rank.prob){
modelstring<-"
model{
 for(i in 1:len){
  logit(p[i]) <- mu[t[i]] + vi[s[i], t[i]]
  r[i] ~ dbin(p[i], totaln[i])
  rhat[i] <- p[i]*totaln[i]
  dev[i] <- 2*(r[i]*(log(r[i]) - log(rhat[i])) +
    (totaln[i] - r[i])*(log(totaln[i] - r[i]) - log(totaln[i] - rhat[i])))
 }
 totresdev <- sum(dev[])
 for(j in 1:nstudy){
  vi[j, 1:ntrt] ~ dmnorm(zeros[1:ntrt], T[1:ntrt, 1:ntrt])
 }
 for(j in 1:ntrt){
  AR[j] <- 1/(1 + exp(-mu[j]/sqrt(1 + (16*sqrt(3)/(15*3.1415926))^2*invT[j,j])))
  mu[j] ~ dnorm(0, 0.001)
  sigma[j] <- sqrt(invT[j,j])
 }
 invT[1:ntrt, 1:ntrt] <- inverse(T[,])
 T[1:ntrt, 1:ntrt] ~ dwish(I[1:ntrt, 1:ntrt], ntrt + 1)
 for(j in 1:ntrt){        
  for(k in 1:ntrt){
   LRR[j,k] <- log(RR[j,k])
   LOR[j,k] <- log(OR[j,k])
   LOR.med[j,k] <- mu[j] - mu[k]
   RR[j,k] <- AR[j]/AR[k]
   RD[j,k] <- AR[j] - AR[k]
   OR[j,k] <- AR[j]/(1 - AR[j])/AR[k]*(1 - AR[k])
   OR.med[j,k] <- exp(LOR.med[j,k])
  }
 }
}
"
}

if(prior.type == "chol" & rank.prob){
modelstring<-"
model{
 for(i in 1:len){
  logit(p[i]) <- mu[t[i]] + vi[s[i], t[i]]
  r[i] ~ dbin(p[i], totaln[i])
  rhat[i] <- p[i]*totaln[i]
  dev[i] <- 2*(r[i]*(log(r[i]) - log(rhat[i])) +
    (totaln[i] - r[i])*(log(totaln[i] - r[i]) - log(totaln[i] - rhat[i])))
 }
 totresdev <- sum(dev[])
 for(j in 1:nstudy){
  vi[j, 1:ntrt] ~ dmnorm(zeros[1:ntrt], invSig[1:ntrt, 1:ntrt])
 }
 for(j in 1:ntrt){
  AR[j] <- 1/(1 + exp(-mu[j]/sqrt(1 + (16*sqrt(3)/(15*3.1415926))^2*pow(sigma[j], 2))))
  mu[j] ~ dnorm(0, 0.001)
 }
 invSig[1:ntrt, 1:ntrt] <- inverse(Sig[,])
 for(i in 1:ntrt){
  for(j in 1:ntrt){
    Sig[i,j] <- sigma[i]*sigma[j]*R[i,j]
  }
 }
 R[1:ntrt, 1:ntrt] <- L[1:ntrt, 1:ntrt] %*% t(L[1:ntrt, 1:ntrt])
 L[1,1] <- 1
 for(j in 2:ntrt){
  L[1,j] <- 0
 }
 for(i in 2:(ntrt - 1)){
  L[i,1] <- cos(psi[i-1,1])
  for(j in 2:(i - 1)){
   L[i,j] <- prod(sin(psi[i-1, 1:(j - 1)]))*cos(psi[i - 1, j])
  }
  L[i,i] <- prod(sin(psi[i-1, 1:(i - 1)]))
  for(j in (i + 1):ntrt){
   L[i,j] <- 0
  }
 }
 L[ntrt,1] <- cos(psi[ntrt-1,1])
 for(j in 2:(ntrt - 1)){
  L[ntrt,j] <- prod(sin(psi[ntrt-1, 1:(j-1)]))*cos(psi[ntrt - 1, j])
 }
 L[ntrt, ntrt] <- prod(sin(psi[ntrt - 1, 1:(ntrt - 1)]))
 for(i in 1:(ntrt - 1)){
  for(j in 1:(ntrt - 1)){
   psi[i, j] ~ dunif(0.01, 3.13)
  }
 }
 for(i in 1:ntrt){
  sigma[i] ~ dunif(0.0001, c)
 }
 for(j in 1:ntrt){        
  for(k in 1:ntrt){
   LRR[j,k] <- log(RR[j,k])
   LOR[j,k] <- log(OR[j,k])
   LOR.med[j,k] <- mu[j] - mu[k]
   RR[j,k] <- AR[j]/AR[k]
   RD[j,k] <- AR[j] - AR[k]
   OR[j,k] <- AR[j]/(1 - AR[j])/AR[k]*(1 - AR[k])
   OR.med[j,k] <- exp(LOR.med[j,k])
  }
 }
 rk[1:ntrt] <- (ntrt + 1 - rank(AR[]))*ifelse(higher.better, 1, 0) + (rank(AR[]))*ifelse(higher.better, 0, 1)
 rk.med[1:ntrt] <- (ntrt + 1 - rank(mu[]))*ifelse(higher.better, 1, 0) + (rank(mu[]))*ifelse(higher.better, 0, 1)
 for(i in 1:ntrt){
  rank.prob[1:ntrt, i] <- equals(rk[], i)
  rank.prob.med[1:ntrt, i] <- equals(rk.med[], i)
 }
}
"
}

if(prior.type == "chol" & !rank.prob){
modelstring<-"
model{
 for(i in 1:len){
  logit(p[i]) <- mu[t[i]] + vi[s[i], t[i]]
  r[i] ~ dbin(p[i], totaln[i])
  rhat[i] <- p[i]*totaln[i]
  dev[i] <- 2*(r[i]*(log(r[i]) - log(rhat[i])) +
    (totaln[i] - r[i])*(log(totaln[i] - r[i]) - log(totaln[i] - rhat[i])))
 }
 totresdev <- sum(dev[])
 for(j in 1:nstudy){
  vi[j, 1:ntrt] ~ dmnorm(zeros[1:ntrt], invSig[1:ntrt, 1:ntrt])
 }
 for(j in 1:ntrt){
  AR[j] <- 1/(1 + exp(-mu[j]/sqrt(1 + (16*sqrt(3)/(15*3.1415926))^2*pow(sigma[j], 2))))
  mu[j] ~ dnorm(0, 0.001)
 }
 invSig[1:ntrt, 1:ntrt] <- inverse(Sig[,])
 for(i in 1:ntrt){
  for(j in 1:ntrt){
    Sig[i,j] <- sigma[i]*sigma[j]*R[i,j]
  }
 }
 R[1:ntrt, 1:ntrt] <- L[1:ntrt, 1:ntrt] %*% t(L[1:ntrt, 1:ntrt])
 L[1,1] <- 1
 for(j in 2:ntrt){
  L[1,j] <- 0
 }
 for(i in 2:(ntrt - 1)){
  L[i,1] <- cos(psi[i - 1, 1])
  for(j in 2:(i - 1)){
   L[i,j] <- prod(sin(psi[i - 1, 1:(j - 1)]))*cos(psi[i - 1, j])
  }
  L[i,i] <- prod(sin(psi[i - 1, 1:(i - 1)]))
  for(j in (i + 1):ntrt){
   L[i,j] <- 0
  }
 }
 L[ntrt,1] <- cos(psi[ntrt - 1, 1])
 for(j in 2:(ntrt - 1)){
  L[ntrt,j] <- prod(sin(psi[ntrt - 1, 1:(j - 1)]))*cos(psi[ntrt - 1, j])
 }
 L[ntrt,ntrt] <- prod(sin(psi[ntrt - 1, 1:(ntrt - 1)]))
 for(i in 1:(ntrt - 1)){
  for(j in 1:(ntrt - 1)){
   psi[i, j] ~ dunif(0.01, 3.13)
  }
 }
 for(i in 1:ntrt){
  sigma[i] ~ dunif(0.0001, c)
 }
 for(j in 1:ntrt){        
  for(k in 1:ntrt){
   LRR[j,k] <- log(RR[j,k])
   LOR[j,k] <- log(OR[j,k])
   LOR.med[j,k] <- mu[j] - mu[k]
   RR[j,k] <- AR[j]/AR[k]
   RD[j,k] <- AR[j] - AR[k]
   OR[j,k] <- AR[j]/(1 - AR[j])/AR[k]*(1 - AR[k])
   OR.med[j,k] <- exp(LOR.med[j,k])
  }
 }
}
"
}

if(!is.element(prior.type,c("invwishart", "chol"))){
  stop("specified prior type is wrong.")
}

return(modelstring)
}