\name{nma.ab.bin}
\alias{nma.ab.bin}
\title{
Arm-Based Network Meta-Analysis for Binary Outcomes
}
\description{
Performs the arm-based network meta-analysis proposed by Zhang et al.{\sspace}(2014), including estimating the absolute risk (AR), risk difference (RD), odds ratio (OR), risk ratio (RR), log odds ratio (LOR), and log risk ratio (LRR).
}
\usage{
nma.ab.bin(s.id, t.id, event.n, total.n, data, trtname,
           param = c("AR", "LOR", "LRR", "RD", "rank.prob"),
           model = "het_cor", link = "probit",
           prior.type, a = 0.001, b = 0.001, c = 10,
           higher.better = FALSE, digits = 4, n.adapt = 5000,
           n.iter = 100000, n.burnin = floor(n.iter/2), n.chains = 3,
           n.thin = max(1, floor((n.iter - n.burnin)/100000)),
           conv.diag = FALSE, trace = NULL, dic = FALSE, postdens = FALSE,
           mcmc.samples = FALSE)
}
\arguments{
  \item{s.id}{
  a numeric or character vector specifying study ID, or the corresponding column name in the argument \code{data}.
}
  \item{t.id}{
  a numeric or character vector specifying treatment ID, or the corresponding column name in the argument \code{data}.
}
  \item{event.n}{
  a numeric vector of non-negative integers, specifying event number for a certain treatment in the corresponding study, or the corresponding column name in the argument \code{data}.
}
  \item{total.n}{
  a numeric vector of positive integers, specifying total number of participants for a certain treatment in the corresponding study, or the corresponding column name in the argument \code{data}.
}
  \item{data}{
  an optional data frame containing the dataset for network meta-analysis. If \code{data} is specified, the previous arguments, \code{s.id}, \code{t.id}, \code{event.n}, and \code{total.n}, should be specified as the corresponding column names in \code{data}; otherwise, the previous arguments use environment variables.
}
  \item{trtname}{
  a vector of character strings specifying the treatment names for the corresponding treatment IDs according to their order in \code{t.id}. If not specified, \code{t.id} is used as treatment names.
}
  \item{param}{
  a vector of character strings specifying the effect measures to be estimated. The default includes \code{"AR"} (absolute risk), \code{"LOR"} (log odds ratio), \code{"LRR"} (log risk ratio), \code{"RD"} (risk difference), \code{"rank.prob"} (treatment rank probability). \code{"AR"} is automatically added into \code{param} even if it is not specified. In addition to the default measures, \code{"OR"} (odds ratio) and \code{"RR"} (risk ratio) can be added into \code{param}. If the logit link is used (\code{link} = \code{"logit"}), \code{"OR.med"} (overall median odds ratio), \code{"LOR.med"} (overall median log odds ratio), and \code{"rank.prob.med"} (treatment rank probability based on overall median odds ratio) can be also added; see "Details". If \code{model} is \code{"hom_eqcor"} or \code{"het_eqcor"}, \code{param} can also include \code{"\rho"} (the common correlation coefficient between treatments).
}
  \item{model}{
  a character string specifying which Bayesian hierarchical model to be applied in the arm-based network meta-analysis. This argument can be set as \code{"hom_eqcor"}, \code{"het_eqcor"}, or \code{"het_cor"} (the default). See "Details" for the models.
}
  \item{link}{
  a character string specifying the link function in the Bayesian hierarchical model for binary outcomes. It can be either \code{"probit"} (the default) or \code{"logit"}.
}
  \item{prior.type}{
  prior distribution of variances and/or covariances of random effects. If \code{model} is \code{"hom_eqcor"} or \code{"het_eqcor"}, it can be set as \code{"unif"} (uniform prior for standard deviations, the default) or \code{"invgamma"} (inverse gamma prior for variances). If \code{model} is \code{"het_cor"}, \code{prior.type} can be \code{"invwishart"} (the default) or \code{"chol"}. Specifying \code{"invwishart"} yields an inverse-Wishart prior for the variance-covariance matrix of random effects; by specifying \code{"chol"}, non-informative priors are assigned to variance and correlation components using the separation strategy by Cholesky decomposition. See "Details".
}
  \item{a, b}{
  positive numbers, specifying the shape and scale parameters of inverse gamma priors for variance(s) of random effects if using \code{prior.type} as \code{"invgamma"} for model \code{"hom_eqcor"} or \code{"het_eqcor"}. The default values for both parameters are 0.001.
}
  \item{c}{
  positive number, specifying the upper bound of uniform prior for standard deviation(s) of random effects if using \code{prior.type} as \code{"unif"} for model \code{"hom_eqcor"} or \code{"het_eqcor"}. The default is 10.
}
  \item{higher.better}{
  an optional logical value which needs to be specified when estimating the treatment rank probabilities (i.e., \code{"rank.prob"} is included in the argument \code{param}). \code{TRUE} indicates that a higher event rate implies a better treatment, and vice versa. The default is \code{FALSE}.
}
  \item{digits}{
  a positive integer specifying the digits after the decimal point for the effect measure estimates. The default is 4.
}
  \item{n.adapt}{
  the number of iterations for adaptation in Markov chain Monte Carlo (MCMC) algorithm. The default is 5,000. If a warning "adaptation incomplete" appears, users may increase \code{n.adapt}. This argument and the following \code{n.iter}, \code{n.burnin}, \code{n.chains}, and \code{n.thin} are passed to the functions in R package \pkg{rjags}.
}
  \item{n.iter}{
  the total number of iterations in each MCMC chain. The default is 100,000.
}
  \item{n.burnin}{
  the number of iterations for burn-in period. The default is \code{n.iter/2}.
}
  \item{n.chains}{
  the number of MCMC chains. The default is 3.
}
  \item{n.thin}{
  a positive integer specifying the thinning rate. The default is the thinning rate which yields no more than 100,000 iterations remaining in each chain.
}
  \item{conv.diag}{
  a logical value indicating whether to perform MCMC convergence diagnostic. The default is \code{FALSE}. If \code{TRUE}, \code{n.chains} must be greater than 1, and a .txt file, which contains the point estimates of the potential scale reduction factor and their upper confidence limits (see Gelman and Rubin 1992), will be written in users' current working directory.
}
  \item{trace}{
  a vector of character strings of effect measures. The character strings should be selected from those specified in \code{param} (except \code{"rank.prob"}), and trace plots of the specified effect measures will be saved in users' current working directory. The default is not drawing trace plots (\code{NULL}).
}
  \item{dic}{
  a logical value indicating whether to calculate the deviance information criterion (DIC) value. The default is \code{FALSE}.
}
  \item{postdens}{
  a logical value indicating whether to draw the posterior density plots for treatment-specific absolue risks (ARs). If \code{TRUE}, a .pdf file containing the density plot will be saved in users' current working directory. The default is \code{FALSE}.
}
  \item{mcmc.samples}{
  a logical value indicating whether to save MCMC posterior samples in the output object. The default is \code{FALSE}.
}
}
\details{
Suppose that a network meta-analysis collects \eqn{I} studies on \eqn{K} treatments, where each study investigates a subset of the \eqn{K} treatments. Label the studies from \eqn{i = 1} to \eqn{I} and the treatments from \eqn{k = 1} to \eqn{K}. Let \eqn{T_{i}} be the subset of the \eqn{K} treatments that is compared in the \eqn{i}th study. Also, in the \eqn{i}th study, let \eqn{n_{ik}} be the number of participants allocated to treatment group \eqn{k} (\eqn{k \in T_{i}}), and \eqn{y_{ik}} be the number of events. The arm-based model is constructed as (Zhang et al.{\sspace}2014):
\deqn{y_{ik} \sim Bin (n_{ik}, p_{ik}) \qquad k \in T_{i};}
\deqn{\Phi^{-1} (p_{ik}) = \mu_{k} + \nu_{ik};}
\deqn{(\nu_{i1}, \nu_{i2}, \ldots, \nu_{iK})^{T} \sim N (\boldsymbol{0}, \mathbf{\Sigma}_{K}),}
where \eqn{\Phi (\cdot)} is the standard normal cumulative distribution function, and \eqn{\mathbf{\Sigma}_{K}} is a \eqn{K \times K} positive definite covariance matrix. The \eqn{\mu_{k}}'s are treatment-specific fixed effects, and the random effects \eqn{\nu_{ik}} are correlated within each study with the covariance matrix \eqn{\mathbf{\Sigma}_{K}}. The marginal absolute risk of treatment \eqn{k} is \eqn{p_{k} = E[p_{ik}]}; other effect measures are calculated based on these absolute risks.

An unstructured covariance matrix \eqn{\mathbf{\Sigma}_{K}} in the model above corresponds to \code{model} = \code{"het_cor"}. The inverse-Wishart prior can be assigned to \eqn{\mathbf{\Sigma}_{K}}. Alternatively, using the separation strategy by Cholesky decomposition (\code{prior.type} = \code{"chol"}), uniform priors \eqn{U(0, c)} are assigned to the standard deviations in \eqn{\mathbf{\Sigma}_{K}} and non-informative priors are assigned to the correlation components (Barnard et al., 2000; Lu and Ades, 2009; Wei and Higgins, 2013; Lin and Chu, 2018). Denote \eqn{\sigma_{k}} as the standard deviation of \eqn{\nu_{ik}} and \eqn{\mathbf{D} = diag(\sigma_{1}, \ldots, \sigma_{K})}, then the correlation matrix is \eqn{\mathbf{R}_{K} = \mathbf{D}^{-1} \mathbf{\Sigma}_{K} \mathbf{D}^{-1}}.  If we assume that all of the off-diagonal elements in \eqn{\mathbf{R}_{K}} are equal, say to \eqn{\rho}, then this model corresponds to \code{model} = \code{"het_eqcor"}. If we further assume the homogeneity of variances of the random effects, that is, \eqn{\sigma_{k} = \sigma} for \eqn{k = 1, 2, \ldots, K}, then the model is \code{"hom_eqcor"}. In addition, for the models \code{"hom_eqcor"} and \code{"het_eqcor"}, setting \code{prior.type} as \code{"invgamma"} implies using inverse-gamma priors with shape and scale parameters, \eqn{a} and \eqn{b}, for \eqn{\sigma_{k}^2} or \eqn{\sigma^2}, and \code{"unif"} implies uniform priors \eqn{U(0, c)} for \eqn{\sigma_{k}} or \eqn{\sigma}.

In addition to the probit link as used in the model above, one may also use the logit link (Chu et al., 2012), which is adopted more commonly in contrast-based models. Using the logit link, \eqn{\mu_k} represents the overall median log odds across studies; thus, \eqn{\mu_k - \mu_h} is the overall median log odds ratio between treatments \eqn{h} and \eqn{k}, \eqn{\exp(\mu_k - \mu_h)} is the overall median odds ratio, and they futher yield the rank probabilities. The contrast-based models usually report such overall median (log) odds ratios, rather than the marginal effect measures obtained from the arm-based models. When the argument \code{link} is \code{"logit"}, users can additionally specify \code{"LOR.med"}, \code{"OR.med"}, and \code{"rank.prob.med"} to obtain these results based on the overall median odds ratios.
}
\value{
\code{nma.ab.bin} returns a list with estimates of effect measures specified in \code{param}. If the argument \code{dic} = \code{TRUE}, the deviance information criterion (DIC) statistic will be returned in the output list. In addition, if \code{conv.diag} = \code{TRUE}, a .txt file containing the point estimates of the potential scale reduction factor and their upper confidence limits by Gelman and Rubin (1992) will be saved in users' current working directory. If \code{postdens} = \code{TRUE}, the posterior densities of treatment-specific absolute risks will be saved as a .pdf file. If \code{trace} is specified, the trace plots are saved as .png files.
}
\references{
Barnard J, McCulloch R, Meng XL (2000). "Modeling covariance matrices in terms of standard deviations and correlations, with application to shrinkage." \emph{Statistica Sinica}, \bold{10}(4), 1281--1311.

Chu H, Nie L, Chen Y, Huang Y, Sun W (2012). "Bivariate random effects models for meta-analysis of comparative studies with binary outcomes: methods for the absolute risk difference and relative risk." \emph{Statistical Methods in Medical Research}, \bold{21}(6), 621--633. <\doi{10.1177/0962280210393712}>

Gelman A, Rubin DB (1992). "Inference from iterative simulation using multiple sequences." \emph{Statistical Science}, \bold{7}(4), 457--472. <\doi{10.1214/ss/1177011136}>

Lin L, Chu H (2018). "Bayesian multivariate meta-analysis of multiple factors." \emph{Research Synthesis Methods}, \bold{9}(2), 261--272. <\doi{10.1002/jrsm.1293}>

Lin L, Zhang J, Hodges JS, Chu H (2017). "Performing arm-based network meta-analysis in R with the pcnetmeta package." \emph{Journal of Statistical Software}, \bold{80}(5), 1--25. <\doi{10.18637/jss.v080.i05}>

Lu G, Ades AE (2004). "Combination of direct and indirect evidence in mixed treatment comparisons." \emph{Statistics in Medicine}, \bold{23}(20), 3105--3124. <\doi{10.1002/sim.1875}>

Lu G, Ades AE (2009). "Modeling between-trial variance structure in mixed treatment comparisons." \emph{Biostatistics}, \bold{10}(4), 792--805. <\doi{10.1093/biostatistics/kxp032}>

Spiegelhalter DJ, Best NG, Carlin BP, Van Der Linde A (2002). "Bayesian measures of model complexity and fit." \emph{Journal of the Royal Statistical Society, Series B (Statistical Methodology)}, \bold{64}(4), 583--639. <\doi{10.1111/1467-9868.00353}>

Wei Y, Higgins JPT (2013). "Bayesian multivariate meta-analysis with multiple outcomes." \emph{Statistics in Medicine}, \bold{32}(17), 2911--2934. <\doi{10.1002/sim.5745}>

Zhang J, Carlin BP, Neaton JD, Soon GG, Nie L, Kane R, Virnig BA, Chu H (2014). "Network meta-analysis of randomized clinical trials: Reporting the proper summaries." \emph{Clinical Trials}, \bold{11}(2), 246--262. <\doi{10.1177/1740774513498322}>
}

\note{
If there exists a treatment that has no event in all studies, errors may occur. To avoid such errors, users may set the zero event as a small positive number (say 0.5).

Earlier versions (< 4.0.0) of JAGS do not guarantee exact reproducibility of the results. We recommend users to install the latest version (>= 4.0.0) of JAGS so that exact reproducibility can be ensured by specifying certain seeds.
}

\seealso{
\code{\link{nma.ab.cont}}, \code{\link{nma.ab.py}}, \code{\link{nma.ab.followup}}
}
\examples{
data("smoke")

# For the smoke cessation data,
# higher event rate indicates better treatment

# use the model = "het_cor"
#set.seed(1234)
#het.cor.out <- nma.ab.bin(s.id, t.id, r, n, data = smoke,
#  trtname = c("NC", "SH", "IC", "GC"), param = c("AR", "OR", "RR", "LOR",
#  "LRR", "RD", "rank.prob"), model = "het_cor", higher.better = TRUE,
#  n.iter = 200000, n.thin = 1, conv.diag = TRUE, dic = TRUE,
#  trace = c("AR", "LOR"), postdens = TRUE)

# use the model = "hom_eqcor"
# increase n.iter to reach convergence
set.seed(123)
hom.eqcor.out <- nma.ab.bin(s.id, t.id, r, n, data = smoke,
  param = c("AR", "LRR"), model = "hom_eqcor", prior.type = "unif", c = 10,
  higher.better = TRUE, n.adapt = 1000, n.iter = 100, n.chains = 1)

# use the logit link
set.seed(1234)
hom.eqcor.out2 <- nma.ab.bin(s.id, t.id, r, n, data = smoke,
  param = c("AR", "OR", "OR.med", "rank.prob", "rank.prob.med"),
  model = "hom_eqcor", link = "logit", prior.type = "unif", c = 10,
  higher.better = TRUE, n.adapt = 1000, n.iter = 100, n.chains = 1)
}
\keyword{network meta-analysis}