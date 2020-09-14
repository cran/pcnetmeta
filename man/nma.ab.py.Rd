\name{nma.ab.py}
\alias{nma.ab.py}
\title{
Arm-Based Network Meta-Analysis for Count Datasets with Exposure Time in Person-Years Reported
}
\description{
Performs the arm-based network meta-analysis for count datasets when the exposure times (in person-years) are reported, which estimates the treatment-specific rate, rate ratio between treatments, and their logarithms.
}
\usage{
nma.ab.py(s.id, t.id, event.n, py, data, trtname,
          param = c("lograte", "logratio", "rank.prob"), model = "het_cor",
          prior.type, a = 0.001, b = 0.001, c = 10, higher.better = FALSE,
          digits = 4, n.adapt = 5000, n.iter = 100000,
          n.burnin = floor(n.iter/2), n.chains = 3,
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
  a numeric vector of non-negative integers, specifying number of events in each study's treatment group, or the corresponding column name in the argument \code{data}.
}
  \item{py}{
  a numeric vector of non-negative numbers, specifying exposure time in person-years in each study's treatment group, or the corresponding column name in the argument \code{data}.
}
  \item{data}{
  an optional data frame containing the dataset of the network meta-analysis. If \code{data} is specified, the previous arguments, \code{s.id}, \code{t.id}, \code{event.n}, and \code{py}, should be specified as the corresponding column names in \code{data}; otherwise, the previous arguments use environment variables.
}
  \item{trtname}{
  a vector of character strings specifying the treatment names for the corresponding treatment IDs according to their order in \code{t.id}. If not specified, \code{t.id} is used as treatment names.
}
  \item{param}{
  a vector of character strings specifying the effect measures to be estimated. The default includes log treatment-specific rate (\code{"lograte"}), log rate ratio (\code{"logratio"}), and treatment rank probability (\code{"rank.prob"}). \code{"lograte"} is automatically added into \code{param} even if it is not specified. In addition to the default measures, treatment-specific rate (\code{"rate"}) and rate ratio (\code{"ratio"}) can be added into the argument \code{param}.
}
  \item{model}{
  a character string specifying which Bayesian hierarchical model to be applied in the arm-based network meta-analysis. This argument can be set as \code{"hom_eqcor"}, \code{"het_eqcor"}, or \code{"het_cor"} (the default). See "Details" for the models.
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
  an optional logical value which needs to be specified when estimating the treatment rank probabilities (i.e., "rank.prob" is included in argument \code{param}). \code{TRUE} indicates that a higher treatment-specific rate implies a better treatment, and vice versa. The default is \code{FALSE}.
}
  \item{digits}{
  a positive integer specifying the digits after the decimal point of the effect measures estimates. The default is 4.
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
  a logical value indicating whether to draw the posterior density plots for treatment-specific rates. If \code{TRUE}, a .pdf file containing the density plot will be saved in users' current working directory. The default is \code{FALSE}.
}
  \item{mcmc.samples}{
  a logical value indicating whether to save MCMC posterior samples in the output object. The default is \code{FALSE}.
}
}
\details{
Suppose that a network meta-analysis collects \eqn{I} studies on \eqn{K} treatments, where each study investigates a subset of the \eqn{K} treatments. The exposure time in person-years and the count of events in each treatment group are reported. Label the studies from \eqn{i = 1} to \eqn{I} and the treatments from \eqn{k = 1} to \eqn{K}. Let \eqn{T_{i}} be the subset of the \eqn{K} treatments that is compared in the \eqn{i}th study. Also, in the \eqn{i}th study, let \eqn{y_{ik}} be the number of events in treatment group \eqn{k}, and \eqn{E_{ik}} be the corresponding exposure time in person-years. The arm-based network meta-analysis model for these settings is constructed as:
\deqn{y_{ik} \sim Pois (E_{ik} \lambda_{ik}) \qquad k \in T_{i};}
\deqn{\log (\lambda_{ik}) = \mu_{k} + \nu_{ik};}
\deqn{(\nu_{i1}, \nu_{i2}, \ldots, \nu_{iK})^{T} \sim N (\boldsymbol{0}, \mathbf{\Sigma}_{K}),}
where \eqn{\mathbf{\Sigma}_{K}} is a \eqn{K \times K} positive definite correlation matrix. The \eqn{\mu_{k}}'s are treatment-specific fixed effects, and the random effects \eqn{\nu_{ik}} are correlated within each study with the covariance matrix \eqn{\mathbf{\Sigma}_{K}}.

An unstructured covariance matrix \eqn{\mathbf{\Sigma}_{K}} in the model above corresponds to \code{model} = \code{"het_cor"}. The inverse-Wishart prior can be assigned to \eqn{\mathbf{\Sigma}_{K}}. Alternatively, using the separation strategy by Cholesky decomposition (\code{prior.type} = \code{"chol"}), uniform priors \eqn{U(0, c)} are assigned to the standard deviations in \eqn{\mathbf{\Sigma}_{K}} and non-informative priors are assigned to the correlation components (Barnard et al., 2000; Lu and Ades, 2009; Wei and Higgins, 2013; Lin and Chu, 2018). Denote \eqn{\sigma_{k}} as the standard deviation of \eqn{\nu_{ik}} and \eqn{\mathbf{D} = diag(\sigma_{1}, \ldots, \sigma_{K})}, then the correlation matrix is \eqn{\mathbf{R}_{K} = \mathbf{D}^{-1} \mathbf{\Sigma}_{K} \mathbf{D}^{-1}}.  If we assume that all of the off-diagonal elements in \eqn{\mathbf{R}_{K}} are equal, say to \eqn{\rho}, then this model corresponds to \code{model} = \code{"het_eqcor"}. If we further assume the homogeneity of variances of the random effects, that is, \eqn{\sigma_{k} = \sigma} for \eqn{k = 1, 2, \ldots, K}, then the model is \code{"hom_eqcor"}. In addition, for the models \code{"hom_eqcor"} and \code{"het_eqcor"}, setting \code{prior.type} as \code{"invgamma"} implies using inverse-gamma priors with shape and scale parameters, \eqn{a} and \eqn{b}, for \eqn{\sigma_{k}^2} or \eqn{\sigma^2}, and \code{"unif"} implies uniform priors \eqn{U(0, c)} for \eqn{\sigma_{k}} or \eqn{\sigma}.
}
\value{
\code{nma.ab.py} returns a list with estimates of effect measures specified in \code{param}. If the argument \code{dic} = \code{TRUE}, the deviance information criterion (DIC) statistic will be returned in the output list. In addition, if \code{conv.diag} = \code{TRUE}, a .txt file containing the point estimates of the potential scale reduction factor and their upper confidence limits by Gelman and Rubin (1992) will be saved in users' current working directory. If \code{postdens} = \code{TRUE}, the posterior densities of treatment-specific absolute risks will be saved as a .pdf file. If \code{trace} is specified, the trace plots are saved as .png files.
}
\references{
Barnard J, McCulloch R, Meng XL (2000). "Modeling covariance matrices in terms of standard deviations and correlations, with application to shrinkage." \emph{Statistica Sinica}, \bold{10}(4), 1281--1311.

Dias S, Sutton AJ, Ades AE, Welton NJ (2013). "Evidence synthesis for decision making 2: a generalized linear modeling framework for pairwise and network meta-analysis of randomized controlled trials." \emph{Medical Decision Making}, \bold{33}(5), 607--617. <\doi{10.1177/0272989X12458724}>

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
Earlier versions (< 4.0.0) of JAGS do not guarantee exact reproducibility of the results. We recommend users to install the latest version (>= 4.0.0) of JAGS so that exact reproducibility can be ensured by specifying certain seeds.
}

\seealso{
\code{\link{nma.ab.bin}}, \code{\link{nma.ab.cont}}, \code{\link{nma.ab.followup}}
}
\examples{
#data("dietaryfat")
# increase n.iter to reach convergence of MCMC
# increase n.adapt to enhance efficiency
#set.seed(1234)
#py.out <- nma.ab.py(s.id, t.id, r, py, data = dietaryfat, model = "het_cor",
#  n.adapt = 300, n.iter = 100, n.chains = 1)
}
\keyword{network meta-analysis}