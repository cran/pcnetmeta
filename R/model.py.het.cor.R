model.py.het.cor <- function(prior.type = "invwishart", rank.prob = TRUE){
if(prior.type == "invwishart" & rank.prob){
modelstring<-"
model{
 for(i in 1:len){
  y[i] ~ dpois(py[i]*lambda[i])
  lambda[i] <- exp(mu[t[i]] + vi[s[i], t[i]])
 }
 for(j in 1:nstudy){
  vi[j, 1:ntrt] ~ dmnorm(zeros[1:ntrt], T[1:ntrt, 1:ntrt])
 }
 for(j in 1:ntrt){
  rate[j] <- exp(mu[j] + invT[j,j]/2)
  lograte[j] <- log(rate[j])
  mu[j] ~ dnorm(0, 0.001)
  sigma[j] <- sqrt(invT[j,j])
 }
 invT[1:ntrt, 1:ntrt] <- inverse(T[,])
 T[1:ntrt, 1:ntrt] ~ dwish(I[1:ntrt, 1:ntrt], ntrt + 1)
 for(j in 1:ntrt){        
  for(k in 1:ntrt){
   ratio[j,k] <- rate[j]/rate[k]
   logratio[j,k] <- log(ratio[j,k])
  }
 }
 rk[1:ntrt] <- (ntrt + 1 - rank(rate[]))*ifelse(higher.better, 1, 0) + (rank(rate[]))*ifelse(higher.better, 0, 1)
 for(i in 1:ntrt){
  rank.prob[1:ntrt,i] <- equals(rk[], i)
 }
}
"
}

if(prior.type == "invwishart" & !rank.prob){
modelstring<-"
model{
 for(i in 1:len){
  y[i] ~ dpois(py[i]*lambda[i])
  lambda[i] <- exp(mu[t[i]] + vi[s[i], t[i]])
 }
 for(j in 1:nstudy){
  vi[j,1:ntrt] ~ dmnorm(zeros[1:ntrt], T[1:ntrt, 1:ntrt])
 }
 for(j in 1:ntrt){
  rate[j] <- exp(mu[j] + invT[j,j]/2)
  lograte[j] <- log(rate[j])
  mu[j] ~ dnorm(0, 0.001)
  sigma[j] <- sqrt(invT[j,j])
 }
 invT[1:ntrt, 1:ntrt] <- inverse(T[,])
 T[1:ntrt, 1:ntrt] ~ dwish(I[1:ntrt, 1:ntrt], ntrt + 1)
 for(j in 1:ntrt){        
  for(k in 1:ntrt){
   ratio[j,k] <- rate[j]/rate[k]
   logratio[j,k] <- log(ratio[j,k])
  }
 }
}
"
}

if(prior.type == "chol" & rank.prob){
modelstring<-"
model{
 for(i in 1:len){
  y[i] ~ dpois(py[i]*lambda[i])
  lambda[i] <- exp(mu[t[i]] + vi[s[i], t[i]])
 }
 for(j in 1:nstudy){
  vi[j,1:ntrt] ~ dmnorm(zeros[1:ntrt], invSig[1:ntrt, 1:ntrt])
 }
 for(j in 1:ntrt){
  rate[j] <- exp(mu[j] + pow(sigma[j], 2)/2)
  lograte[j] <- log(rate[j])
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
   ratio[j,k] <- rate[j]/rate[k]
   logratio[j,k] <- log(ratio[j,k])
  }
 }
 rk[1:ntrt] <- (ntrt + 1 - rank(rate[]))*ifelse(higher.better, 1, 0) + (rank(rate[]))*ifelse(higher.better, 0, 1)
 for(i in 1:ntrt){
  rank.prob[1:ntrt,i] <- equals(rk[], i)
 }
}
"
}

if(prior.type == "chol" & !rank.prob){
modelstring<-"
model{
 for(i in 1:len){
  y[i] ~ dpois(py[i]*lambda[i])
  lambda[i] <- exp(mu[t[i]] + vi[s[i], t[i]])
 }
 for(j in 1:nstudy){
  vi[j,1:ntrt] ~ dmnorm(zeros[1:ntrt], invSig[1:ntrt, 1:ntrt])
 }
 for(j in 1:ntrt){
  rate[j] <- exp(mu[j] + pow(sigma[j], 2)/2)
  lograte[j] <- log(rate[j])
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
   ratio[j,k] <- rate[j]/rate[k]
   logratio[j,k] <- log(ratio[j,k])
  }
 }
}
"
}

if(!is.element(prior.type, c("invwishart", "chol"))){
  stop("specified prior type is wrong.")
}

return(modelstring)
}