\name{nma.ab.py}
\alias{nma.ab.py}
\title{
Network Meta-Analysis Using Arm-Based Method for Binary Outcomes with Exposure Time in Person-Years
}
\description{
\code{nma.ab.py} conducts an arm-based network meta-analysis for binary outcomes with exposure time in person-years. It can be applied to estimate hazard rate of each treatment, hazard ratio between treatments, and their logarithms. Also, it can provide deviance information criterion (DIC) statistics for evaluating the goodness of fit; give trace plots and Gelman-Rubin diagnostics to check the MCMC convergence; generate posterior density plot for hazard rates of multiple treatments.
}
\usage{
nma.ab.py(s.id, t.id, event.n, py, trtname, model = "het_cor",
          prior.type, a = 0.001, b = 0.001, c = 10,
          param = c("rate", "logratio", "best"), higher.better = FALSE,
          digits = 4, n.adapt = 5000, n.iter = 100000,
          n.burnin = floor(n.iter/2), n.chains = 3,
          n.thin = max(1, floor((n.iter - n.burnin)/100000)),
          conv.diag = FALSE, trace = "", dic = FALSE, postdens = FALSE)
}
\arguments{
  \item{s.id}{
  a numeric or character vector indicating study ID.
}
  \item{t.id}{
  a numeric or character vector indicating treatment ID.
}
  \item{event.n}{
  a numeric vector of non-negative integers, indicating number of events in each study's treatment group.
}
  \item{py}{
  a numeric vector of non-negative numbers, indicating exposure time in person-years in each study's treatment group.
}
  \item{trtname}{
  a vector of character string indicating the treatment names for the corresponding treatment IDs according their order in \code{t.id}. If not specified, \code{t.id} is used as treatment names.
}
  \item{model}{
  a character string indicating which Bayesian hierarchical model to be applied in the arm-based network meta-analysis. This argument can be set as \code{"hom"}, \code{"het_ind"}, or \code{"het_cor"} (default). See "Details" for the models.
}
  \item{prior.type}{
  prior distribution of variances or covariances of random effects. If \code{model = "hom"} or \code{"het_ind"}, it can be set as \code{"unif"} (uniform prior for standard deviation, default) or \code{"invgamma"} (inverse gamma prior for variances). If \code{model = "het_cor"}, \code{prior.type} is automatically set as \code{"invwishart"} to yield inverse-Wishart prior for the variance-covariance matrix of random effects. See "Details".
}
  \item{a, b}{
  positive numbers, specifying the shape and scale parameters of inverse gamma priors for variance(s) of random effects if using \code{prior.type = "invgamma"} for model \code{"hom"} or \code{"het_ind"}. The defaults for both parameters are 0.001.
}
  \item{c}{
  positive number, specifying the upper bound of uniform prior for variance(s) of random effects if using \code{prior.type = "unif"} for model \code{"hom"} or \code{"het_ind"}. The default is 10.
}
  \item{param}{
  a vector of character string indicating the effect sizes to be estimated. The default includes hazard rate ("rate"), log hazard ratio ("logratio"), and probability of being the best treatment ("best"). "rate" is automatically added into \code{param} even if it is not specified. In addition to the defaults, log hazard rate ("lograte") and hazard ratio ("ratio") can be added into the argument \code{param}.
}
  \item{higher.better}{
  a logical value which is in effect when estimating the probabilities of being the best treatment, i.e., "best" is included in argument \code{param}. \code{TRUE} indicates higher sample mean of the continuous outcomes implying better treatment, that is, the event is "good", and vice versa. The default is \code{FALSE}.
}
  \item{digits}{
  a positive integer specifying the digits after the decimal point of the effect sizes estimations. The default is 4.
}
  \item{n.adapt}{
  the number of iterations for adaptation. The default is 5,000. If a warning "adaptation incomplete" appears, one may increase \code{n.adapt}. This argument and the following \code{n.iter}, \code{n.burnin}, \code{n.chains}, \code{n.thin} are passed to the functions in R package \pkg{rjags}.
}
  \item{n.iter}{
  the total number of iterations in each MCMC chain. The default is 100,000.
}
  \item{n.burnin}{
  the number of iterations for burn-in. The default is \code{n.iter/2}.
}
  \item{n.chains}{
  the number of MCMC chains. The default is 3.
}
  \item{n.thin}{
  a positive integer indicating thinning rate. The default is the thinning rate which yields no more than 100,000 iterations remaining in each chain.
}
  \item{conv.diag}{
  a logical value indicating whether to conduct MCMC convergence diagnostic. The default is \code{FALSE}. If \code{TRUE}, a txt file, which contains the point estimates of the potential scale reduction factor (psrf) and their upper confidence limits (see Gelman and Rubin 1992), would be written in the current working directory.
}
  \item{trace}{
  a vector of character string of effect sizes. The character strings should be selected from those specified by \code{param} (except "best"), and trace plots would be drew for the specified effect sizes and saved in the current working directory. The default is not drawing trace plots (an empty character string).
}
  \item{dic}{
  a logical value indicating whether the deviance information criterion (DIC) to be calculated. The default is \code{FALSE}. If \code{TRUE}, \code{n.chains} must be greater than 1.
}
  \item{postdens}{
  a logical value indicating whether to draw the posterior density plot for hazard rates of multiple treatments. If \code{TRUE}, a pdf file containing the plot would be written in the current working directory. The default is \code{FALSE}.
}
}
\details{
Suppose that a network meta-analysis reviews \eqn{I} studies on \eqn{K} treatments, where each study investigates a subset of the \eqn{K} treatments. The outcome is binary, and the exposure times in person-years are available. Denote the studies from \eqn{i = 1} to \eqn{I} and the treatments from \eqn{k = 1} to \eqn{K}. Let \eqn{T_{i}} be the subset of the \eqn{K} treatments that is compared in the \eqn{i}th study. Also, in the \eqn{i}th study, let \eqn{y_{ik}} be the number of events in treatment group \eqn{k} in the \eqn{i}th study, and \eqn{E_{ik}} is the corresponding exposure time in person-years. The arm-based network meta-analysis model for this data can be constructed as
\deqn{y_{ik} \sim Pois (E_{ik} \lambda_{ik}) \qquad k \in T_{i}}
\deqn{\log (\lambda_{ik}) = \mu_{k} + \sigma_{k} \nu_{ik} \qquad k \in T_{i}}
\deqn{(\nu_{i1}, \nu_{i2}, \ldots, \nu_{iK})^{T} \sim N (\boldsymbol{0}, \mathbf{R}_{K}),}
where \eqn{\mathbf{R}_{K}} is a \eqn{K \times K} positive definite correlation matrix. \eqn{\mu_{k}}'s are the fixed effects for different treatments, and \eqn{\nu_{ik}}'s are random effects.

When we use an inverse-Wishart prior for the variance-covariance matrix of random effects, the above model corresponds to setting the function argument \code{model} as \code{"het_cor"}. If we reduce \eqn{\mathbf{R}_{K}} to a \eqn{K \times K} identity matrix, then it corresponds to model \code{"het_ind"}. If we further let \eqn{\sigma_{k} = \sigma} for \eqn{k = 1, 2, \ldots, K}, then it corresponds to model \code{"hom"}. In addition, for the model \code{"het_ind"} and \code{"hom"}, setting \code{prior.type} as \code{"invgamma"} implies using inverse-gamma priors with shape and scale parameters \eqn{a, b} for \eqn{\sigma_{k}} or \eqn{\sigma}, and \code{"unif"} implies uniform priors \eqn{U(0, c)}.
}
\value{
\code{nma.ab.py} returns a list with estimations of effect sizes specified in \code{param}. Also, if the argument \code{dic} is set as \code{TRUE}, the deviance information criterion (DIC) statistics would be returned in the output list. In addition, if \code{conv.diag} is set as \code{TRUE}, a txt file containing the point estimates of the potential scale reduction factor (psrf) and their upper confidence limits by Gelman and Rubin (1992) would be saved in the current working directory. If \code{postdens} is set as \code{TRUE}, the posterior densities of hazard rates of multiple treatments would be saved as a pdf file. If \code{trace} is specified, the trace plots are saved as png files.
}
\references{
Dias S, Welton NJ, Sutton AJ, and Ades AE (2011). "NICE DSU Technical Support Document 2: a generalised linear modelling framework for pairwise and network meta-analysis of randomised controlled trials." \emph{National Institute for Health and Clinical Excellence}, London, UK.

Gelman A and Rubin, DB (1992). "Inference from iterative simulation using multiple sequences." \emph{Statistical science} 457--472.

Lu G and Ades AE (2004). "Combination of direct and indirect evidence in mixed treatment comparisons." \emph{Stat Med} \bold{23}(20), 3105--24.

Spiegelhalter DJ, Best NG, Carlin BP, and Van Der Linde A (2002). "Bayesian measures of model complexity and fit." \emph{JRSSB} \bold{64}(4), 583--639.

Zhang J, Carlin BP, Neaton JD, Soon GG, Nie L, Kane R, Virnig BA, Chu H (2014). "Network meta-analysis of randomized clinical trials: Reporting the proper summaries." \emph{Clin Trials} \bold{11}(2), 246--62.
}
\seealso{
\code{\link{nma.ab}}, \code{\link{nma.ab.cont}}, \code{\link{nma.ab.followup}}
}
\examples{
data(dietaryfat)

# increase n.iter to reach convergence of MCMC
py.out <- nma.ab.py(dietaryfat$sid, dietaryfat$tid, dietaryfat$y,
  dietaryfat$py, model = "het_cor", n.adapt = 1000, n.iter = 100)
}