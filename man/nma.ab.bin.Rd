\name{nma.ab.bin}
\alias{nma.ab.bin}
\title{
Arm-Based Network Meta-Analysis for Binary Outcomes
}
\description{
\code{nma.ab.bin} performs arm-based network meta-analysis proposed by Zhang et al (2014). It estimates absolute risk (AR), risk difference (RD), odds ratio (OR), risk ratio (RR), log odds ratio (LOR), and log risk ratio (LRR).
}
\usage{
nma.ab.bin(s.id, t.id, event.n, total.n, data, trtname,
           param = c("AR", "LOR", "LRR", "RD", "rank.prob"),
           model = "het_cor", prior.type, a = 0.001, b = 0.001, c = 10,
           higher.better = FALSE, digits = 4, n.adapt = 5000,
           n.iter = 100000, n.burnin = floor(n.iter/2), n.chains = 3,
           n.thin = max(1, floor((n.iter - n.burnin)/100000)),
           conv.diag = FALSE, trace = NULL, dic = FALSE, postdens = FALSE,
           mcmc.samples = FALSE)
}
\arguments{
  \item{s.id}{
  a numeric or character vector indicating study ID, or the corresponding column name in the argument \code{data}.
}
  \item{t.id}{
  a numeric or character vector indicating treatment ID, or the corresponding column name in the argument \code{data}.
}
  \item{event.n}{
  a numeric vector of non-negative integers, indicating event number for a certain treatment in the corresponding study, or the corresponding column name in the argument \code{data}.
}
  \item{total.n}{
  a numeric vector of positive integers, indicating total number of participants for a certain treatment in the corresponding study, or the corresponding column name in the argument \code{data}.
}
  \item{data}{
  an optional data frame containing the dataset for network meta-analysis. If \code{data} is specified, the previous arguments, \code{s.id}, \code{t.id}, \code{event.n}, and \code{total.n}, should be specified as the corresponding column names in \code{data}; otherwise, the previous arguments use environment variables.
}
  \item{trtname}{
  a vector of character strings indicating the treatment names for the corresponding treatment IDs according their order in \code{t.id}. If not specified, \code{t.id} is used as treatment names.
}
  \item{param}{
  a vector of character strings indicating the effect sizes to be estimated. The default includes \code{"AR"} (absolute risk), \code{"LOR"} (log odds ratio), \code{"LRR"} (log risk ratio), \code{"RD"} (risk difference), \code{"rank.prob"} (treatment rank probability). \code{"AR"} is automatically added into \code{param} even if it is not specified. In addition to the defaults, \code{"OR"} (odds ratio) and \code{"RR"} (risk ratio) can be added into \code{param}.
}
  \item{model}{
  a character string indicating which Bayesian hierarchical model to be applied in the arm-based network meta-analysis. This argument can be set as \code{"hom_eqcor"}, \code{"het_eqcor"}, or \code{"het_cor"} (default). See "Details" for the models.
}
  \item{prior.type}{
  prior distribution of variances/covariances of random effects. If \code{model} is \code{"hom_eqcor"} or \code{"het_eqcor"}, it can be set as \code{"unif"} (uniform prior for standard deviation, the default) or \code{"invgamma"} (inverse gamma prior for variances). If \code{model} is \code{"het_cor"}, \code{prior.type} is automatically set as \code{"invwishart"} to yield an inverse-Wishart prior for the variance-covariance matrix of random effects. See "Details".
}
  \item{a, b}{
  positive numbers, specifying the shape and scale parameters of inverse gamma priors for variance(s) of random effects if using \code{prior.type} as \code{"invgamma"} for model \code{"hom_eqcor"} or \code{"het_eqcor"}. The defaults for both parameters are 0.001.
}
  \item{c}{
  positive number, specifying the upper bound of uniform prior for standard deviations of random effects if using \code{prior.type} as \code{"unif"} for model \code{"hom_eqcor"} or \code{"het_eqcor"}. The default is 10.
}
  \item{higher.better}{
  an optional logical value which needs to be specified when estimating the treatment rank probabilities (i.e., \code{"rank.prob"} is included in the argument \code{param}). \code{TRUE} indicates higher event rate implying better treatment, and vice versa. The default is \code{FALSE}.
}
  \item{digits}{
  a positive integer specifying the digits after the decimal point for the effect size estimates. The default is 4.
}
  \item{n.adapt}{
  the number of iterations for adaptation in Markov chain Monte Carlo (MCMC) algorithm. The default is 5,000. If a warning "adaptation incomplete" appears, users may increase \code{n.adapt}. This argument and the following \code{n.iter}, \code{n.burnin}, \code{n.chains}, \code{n.thin} are passed to the functions in R package \pkg{rjags}.
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
  a positive integer indicating thinning rate. The default is the thinning rate which yields no more than 100,000 iterations remaining in each chain.
}
  \item{conv.diag}{
  a logical value indicating whether to perform MCMC convergence diagnostic. The default is \code{FALSE}. If \code{TRUE}, \code{n.chains} must be greater than 1, and a .txt file, which contains the point estimates of the potential scale reduction factor and their upper confidence limits (see Gelman and Rubin 1992), will be written in users' current working directory.
}
  \item{trace}{
  a vector of character strings of effect sizes. The character strings should be selected from those specified in \code{param} (except \code{"rank.prob"}), and trace plots of the specified effect sizes will be saved in users' current working directory. The default is not drawing trace plots (\code{NULL}).
}
  \item{dic}{
  a logical value indicating whether the deviance information criterion (DIC) to be calculated. The default is \code{FALSE}.
}
  \item{postdens}{
  a logical value indicating whether to draw the posterior density plots for treatment-specific absolue risks (ARs). If \code{TRUE}, a .pdf file containing the density plot will be saved in users' current working directory. The default is \code{FALSE}.
}
  \item{mcmc.samples}{
  a logical value indicating whether to save MCMC posterior samples in the output object. The default is \code{FALSE}.
}
}
\details{
Suppose that a network meta-analysis collects \eqn{I} studies on \eqn{K} treatments, where each study investigates a subset of the \eqn{K} treatments. Label the studies from \eqn{i = 1} to \eqn{I} and the treatments from \eqn{k = 1} to \eqn{K}. Let \eqn{T_{i}} be the subset of the \eqn{K} treatments that is compared in the \eqn{i}th study. Also, in the \eqn{i}th study, let \eqn{n_{ik}} be the number of participants allocated to treatment group \eqn{k} (\eqn{k \in T_{i}}), and \eqn{y_{ik}} be the number of events. The arm-based model is constructed as (Zhang et al 2014):
\deqn{y_{ik} \sim Bin (n_{ik}, p_{ik}) \qquad k \in T_{i}}
\deqn{\Phi^{-1} (p_{ik}) = \mu_{k} + \nu_{ik}}
\deqn{(\nu_{i1}, \nu_{i2}, \ldots, \nu_{iK})^{T} \sim N (\boldsymbol{0}, \mathbf{\Sigma}_{K}),}
where \eqn{\Phi (\cdot)} is the standard normal cumulative distribution function, and \eqn{\mathbf{\Sigma}_{K}} is a \eqn{K \times K} positive definite covariance matrix. \eqn{\mu_{k}}'s are treatment-specific fixed effects, and random effects \eqn{\nu_{ik}}'s are correlated within each study with the covariance matrix \eqn{\mathbf{\Sigma}_{K}}.

Using an inverse-Wishart prior for \eqn{\mathbf{\Sigma}_{K}}, the above model corresponds to \code{model} = \code{"het_cor"}. Denote \eqn{\sigma_{k}} as the standard deviation of \eqn{\nu_{ik}} and \eqn{\mathbf{D} = diag(\sigma_{1}, \ldots, \sigma_{K})}, then the correlation matrix \eqn{\mathbf{R}_{K} = \mathbf{D}^{-1} \mathbf{\Sigma}_{K} \mathbf{D}^{-1}}.  If we assume that all of the off-diagonal elements in \eqn{\mathbf{R}_{K}} are equal, say to \eqn{\rho}, then this model corresponds to \code{model} = \code{"het_eqcor"}. If we further assume the homogeneity of variances of the random effects, that is, \eqn{\sigma_{k} = \sigma} for \eqn{k = 1, 2, \ldots, K}, then the model is \code{"hom_eqcor"}. In addition, for the models \code{"hom_eqcor"} and \code{"het_eqcor"}, setting \code{prior.type} as \code{"invgamma"} implies using inverse-gamma priors with shape and scale parameters \eqn{a, b} for \eqn{\sigma_{k}^2} or \eqn{\sigma^2}, and \code{"unif"} implies uniform priors \eqn{U(0, c)} for \eqn{\sigma_{k}} or \eqn{\sigma}.
}
\value{
\code{nma.ab.bin} returns a list with estimates of effect sizes specified in \code{param}. If the argument \code{dic} = \code{TRUE}, the deviance information criterion (DIC) statistic will be returned in the output list. In addition, if \code{conv.diag} = \code{TRUE}, a .txt file containing the point estimates of the potential scale reduction factor and their upper confidence limits by Gelman and Rubin (1992) will be saved in users' current working directory. If \code{postdens} = \code{TRUE}, the posterior densities of treatment-specific absolute risks will be saved as a .pdf file. If \code{trace} is specified, the trace plots are saved as .png files.
}
\references{
Gelman A and Rubin, DB (1992). "Inference from iterative simulation using multiple sequences." \emph{Stat Sci} \bold{7}(4), 457--472.

Lu G and Ades AE (2004). "Combination of direct and indirect evidence in mixed treatment comparisons." \emph{Stat Med} \bold{23}(20), 3105--24.

Spiegelhalter DJ, Best NG, Carlin BP, and Van Der Linde A (2002). "Bayesian measures of model complexity and fit." \emph{J R Stat Soc Series B Stat Methodol} \bold{64}(4), 583--639.

Zhang J, Carlin BP, Neaton JD, Soon GG, Nie L, Kane R, Virnig BA, Chu H (2014). "Network meta-analysis of randomized clinical trials: Reporting the proper summaries." \emph{Clin Trials} \bold{11}(2), 246--62.
}

\note{
If there exists a treatment with zero event in all studies, this function may give an error. To avoid this, users may set the zero event in any study as a small positive number (say 0.5).

Note that the earlier versions (< 4.0.0) of JAGS does not guarantee exact reproducibility of the results. Therefore, we recommend users install the latest version (>= 4.0.0) of JAGS so that exact reproducibility can be ensured by specifying certain seeds.
}

\seealso{
\code{\link{nma.ab.cont}}, \code{\link{nma.ab.py}}, \code{\link{nma.ab.followup}}
}
\examples{
data(smoke)

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
# increase n.adapt to enhance efficiency
set.seed(1234)
hom.eqcor.out <- nma.ab.bin(s.id, t.id, r, n, data = smoke,
  param = c("AR", "LRR"), model = "hom_eqcor", prior.type = "unif", c = 10,
  higher.better = TRUE, n.adapt = 400, n.iter = 100, n.chains = 1)
}