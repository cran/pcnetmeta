\name{nma.ab}
\alias{nma.ab}
\title{
Arm-based method in network meta-analysis
}
\description{
\code{nma.ab} generates a summary file, containing estimated event rate ("probt"), risk ratio ("RR"), risk difference ("RD"), odds ratio ("OR") etc., using the arm-based approach, proposed by Zhang et al (2013) for network meta-analysis. Also, it can provide DIC statistics for model selection, and give trace plot and posterior density plot for checking if the Markov chain is convergent.
}
\usage{
nma.ab(s.id, t.id, event.n, total.n, o.path = getwd(), f.name = "",
       hom = FALSE, R, param = c("probt","RR","RD","OR","rk","best"),
       ity = "estimate", n.iter = 200000, n.burnin = floor(n.iter/2),
       n.chains = 3, n.thin = max(1, floor((n.iter - n.burnin)/100000)),
       dic = TRUE, trace = FALSE, postdens = FALSE)
}
\arguments{
  \item{s.id}{
  a numeric vector of natural numbers, indicating study ID.
}
  \item{t.id}{
  a numeric vector of natural numbers, indicating treatment ID.
}
  \item{event.n}{
  a numeric vector of non-negative numbers, indicating event number for corresponding treatment in a particular study.
}
  \item{total.n}{
  a numeric vector of non-negative numbers, indicating total number of participants for corresponding treatment in a particular study.
}
  \item{o.path}{
  a character string indicating output path. If not specified, it would be set as the working directory before runing this function.
}
  \item{f.name}{
  a character string indicating the file name of output results, including network meta-analysis summary, DIC statistics (if \code{dic} = \code{TRUE}), trace plot (if \code{trace} = \code{TRUE}) and posterior density plot (if \code{postdens} = \code{TRUE}). The result file would be named as \code{f.name} followed by \code{"Summary"} for the network meta-analysis summary, by \code{"DIC"} for DIC statistics. The default is an empty character.
}
  \item{hom}{
  logical. If \code{FALSE} (the default), this function will use a heterogeneous model; otherwise, a homogeneous model would be applied. See "Details" for the JAGS models.
}
  \item{R}{
  a \code{tN} by \code{tN} covariance matrix for Wishart prior in the heterogeneous model (\code{hom} = \code{FALSE}), where \code{tN} is the number of treatments. The default is a matrix with diagonal elements being 1 and off-diagonal elements being 0.005.
}
  \item{param}{
  a vector of character string indicating the node(s) to be included in the summary result file. The default is "probt" (estimated event rate),"RR","RD","OR","rk" (rank of treatment),"best" (probability of being best treatment). In addition to the default parameters, "mu" (fixed effect), "vi" (random effect), and "sigma" (homogeneous variance for random effect, only if homogeneous model is used, i.e., \code{hom} = \code{TRUE}) could be selected into \code{param}.
}
  \item{ity}{
  initial value type for the fixed effect \code{mu} in the JAGS model. It can be set as \code{"estimate"} (the default) and \code{"same"}. See "Details".
}
  \item{n.iter}{
  the total number of iterations in each chain. The default is 200,000.
}
  \item{n.burnin}{
  the number of iterations for burn-in. The default is the largest integer not greater than \code{n.iter/2}.
}
  \item{n.chains}{
  the number of parallel chains for the model. The default is 3.
}
  \item{n.thin}{
  a positive integer indicating thinning rate.
}
  \item{dic}{
  logical. If \code{TRUE} (the default), the function would generate a file containing the DIC statistics, and a node named "deviance" would be contained in the summary result file; otherwise, the DIC statistics would not be calculated.
}
  \item{trace}{
  logical. If \code{FALSE} (the default), the function would not give the trace plot; otherwise, trace plots for estimated event rate ("probt") would be given, and they would be saved in the working directory where \code{o.path} indicates.
}
  \item{postdens}{
  logical. If \code{FALSE} (the default), the function would not give the posterior density plot; otherwise, posterior density plots for estimated event rate ("probt") would be given, and they would be saved in the working directory where \code{o.paht} indicates.
}
}
\details{
If \code{ity} is specified as \code{"estimate"}, the initial value for fixed effect would be estimated (in probit scale) as the proportion of sum of event numbers for a treatment over sum of total numbers in corresponding study; otherwise, the initial value of \code{mu} in the JAGS model would be set as 0. Setting \code{ity} as \code{"estimate"} may be helpful for the Markov chain converging more quickly.

The heterogeneous model (\code{hom} = \code{FALSE}) defines the following JAGS model:
\preformatted{model{
 for(i in 1:sN){
  p[i] <- phi(mu[t[i]] + vi[s[i], t[i]])
  r[i] ~ dbin(p[i], totaln[i])
 }
 for(j in 1:tS){
  vi[j, 1:tN] ~ dmnorm(mn[1:tN], T[1:tN, 1:tN])
 }
 invT[1:tN, 1:tN] <- inverse(T[,])
 for(j in 1:tN){
  mu[j] ~ dnorm(0, 0.001)
  sigma[j] <- sqrt(invT[j, j])
  probt[j] <- phi(mu[j]/sqrt(1 + invT[j, j]))
 }
 T[1:tN, 1:tN] ~ dwish(R[1:tN, 1:tN], tN)
 for(j in 1:tN){
  for(k in 1:tN){
   RR[j, k] <- probt[j]/probt[k]
   RD[j, k] <- probt[j] - probt[k]
   OR[j, k] <- probt[j]/(1 - probt[j])/probt[k]*(1 - probt[k])
  }
 }
 rk[1:tN] <- tN + 1 - rank(probt[])
 best[1:tN] <- equals(rk[], 1)
}}

The homogeneous model (\code{hom} = \code{TRUE}) defines the following JAGS model:
\preformatted{model{
 for(i in 1:sN){
  p[i] <- phi(mu[t[i]] + sigma*vi[s[i]])
  r[i] ~ dbin(p[i], totaln[i])
 }
 for(j in 1:tS){
  vi[j] ~ dnorm(0, 1)
 }
 sigma <- 1/sqrt(tau)
 tau ~ dgamma(1, 1)
 for(j in 1:tN){
  mu[j] ~ dnorm(0, 0.001)
  probt[j] <- phi(mu[j]/sqrt(1 + 1/tau))
 }
 for(j in 1:tN){
  for(k in 1:tN){
   RR[j, k] <- probt[j]/probt[k]
   RD[j, k] <- probt[j] - probt[k]
   OR[j, k] <- probt[j]/(1 - probt[j])/probt[k]*(1 - probt[k])
  }
 }
 rk[1:tN] <- tN + 1 - rank(probt[])
 best[1:tN] <- equals(rk[], 1)
}}
}
\value{
\code{nma.ab} generates a summary statistics file using the arm-based method. Furthermore, this function would give a DIC statistics file if \code{dic} = \code{TRUE}, a trace plot file if \code{trace} = \code{TRUE}, a posterior density file if \code{postdens} = \code{TRUE}.

In the summary file, each row contains statistics for corresponding OR, RD, RR, estimated event rate ("probt"), rank of treatment ("rk"), probability of being the best treatment ("best"), etc. Note that RR[i, j], RD[i, j] or OR[i, j] means that treatment i is compared with treatment j, e.g., RD[i,j] = probt[i] - probt[j]. The columns show the statistics of these nodes, including mean, standard deviance, 2.5\% percentile, median, and 97.5\% percentile.

The DIC file contains the value of pD and DIC; the trace plot and posterior density plot file contain corresponding plots for the node "probt" (estimated event rate).
}
\references{
Zhang J, Carlin BP, Neaton JD, Soon GG, Nie L, Kane R, Virnig BA and Chu H. "Network meta-analysis of randomized clinical trials: Reporting the proper summaries." \emph{Clinical Trials} 2013.
}
\author{
Lifeng Lin, Jing Zhang, and Haitao Chu.
}
\note{
You may get warning about "adaptation incomplete" if the model adaptation is not enough. In this case, you need to increase the value of \code{n.burnin} to avoid this warning. This function adapts the model for \code{n.burnin} times if \code{n.burnin} > 0; if \code{n.burnin} = 0, model would be adapted for 100 times.

If there exists treatment with all of event numbers equal to 0 or corresponding total participant number, JAGS may come to an error or give incorrect results. To avoid this problem, you need to add (if event rate is 0\%) or deduct (if event rate is 100\%) a small value (e.g., 0.5) to/from those event numbers.
}

%%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%%}
\examples{
data(Ara09)
attach(Ara09)
nma.ab(s.id = Study.ID, t.id = Treatment, event.n = r, total.n = n,
       f.name = "Ara09_", n.iter = 2000, dic = FALSE)
detach(Ara09)

data(Lam07)
attach(Lam07$data)
nma.ab(s.id = Study.ID, t.id = Treatment, event.n = r, total.n = n, hom = TRUE,
       f.name = "Lam07_", ity = "same", n.iter = 2000, n.chains = 5,
       param = c("probt", "RR", "RD", "OR", "rk", "best", "mu", "sigma", "vi"),
       trace = TRUE, postdens = TRUE)
detach(Lam07$data)
}

