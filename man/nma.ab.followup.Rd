\name{nma.ab.followup}
\alias{nma.ab.followup}
\title{
Arm-Based Network Meta-Analysis for Binary Outcomes with Follow-Up Times Reported
}
\description{
\code{nma.ab.followup} performs arm-based network meta-analysis for binary outcomes and the follow-up time of each study is collected. It estimates treatment-specific rate, rate ratio between treatments, and their logarithms.
}
\usage{
nma.ab.followup(s.id, t.id, event.n, total.n, followup, data, trtname,
                param = c("lograte", "logratio", "rank.prob"),
                model = "het_cor", prior.type, a = 0.001, b = 0.001,
                c = 10, higher.better = FALSE, digits = 4, n.adapt = 5000,
                n.iter = 100000, n.burnin = floor(n.iter/2), n.chains = 3,
                n.thin = max(1, floor((n.iter - n.burnin)/100000)),
                conv.diag = FALSE, trace = "", dic = FALSE,
                postdens = FALSE, mcmc.samples = FALSE)
}
\arguments{
  \item{s.id}{
  a numeric or character vector indicating study ID, or the corresponding column name in the argument \code{data}.
}
  \item{t.id}{
  a numeric or character vector indicating treatment ID, or the corresponding column name in the argument \code{data}.
}
  \item{event.n}{
  a numeric vector of non-negative integers, indicating number of events in each study's treatment group, or the corresponding column name in the argument \code{data}.
}
  \item{total.n}{
  a numeric vector of non-negative integers, indicating total number of participants in each study's treatment group, or the corresponding column name in the argument \code{data}.
}
  \item{followup}{
  a numeric vector of positive numbers, indicating follow-up times for different studies, or the corresponding column name in the argument \code{data}.
}
  \item{data}{
  an optional data frame containing the dataset for network meta-analysis. If \code{data} is specified, the previous arguments, \code{s.id}, \code{t.id}, \code{event.n}, \code{total.n}, and \code{followup}, should be specified as the corresponding column names in \code{data}; otherwise, the previous arguments use environment variables.
}
  \item{trtname}{
  a vector of character strings indicating the treatment names for the corresponding treatment IDs according their order in \code{t.id}. If not specified, \code{t.id} is used as treatment names.
}
  \item{param}{
  a vector of character strings indicating the effect sizes to be estimated. The default includes log treatment-specific rate (\code{"lograte"}), log rate ratio (\code{"logratio"}), and treatment rank probability (\code{"rank.prob"}). \code{"lograte"} is automatically added into \code{param} even if it is not specified. In addition to the defaults, treatment-specific rate (\code{"rate"}) and rate ratio (\code{"ratio"}) can be added into the argument \code{param}.
}
  \item{model}{
  a character string indicating which Bayesian hierarchical model to be applied in the arm-based network meta-analysis. This argument can be set as \code{"hom_ind"}, \code{"het_ind"}, \code{"hom_eqcor"}, \code{"het_eqcor"}, or \code{"het_cor"} (default). See "Details" for the models.
}
  \item{prior.type}{
  prior distribution of variances/covariances of random effects. If \code{model} is \code{"hom_eqcor"}, \code{"het_eqcor"}, \code{"hom_ind"}, or \code{"het_ind"}, it can be set as \code{"unif"} (uniform prior for standard deviation, the default) or \code{"invgamma"} (inverse gamma prior for variances). If \code{model} is \code{"het_cor"}, \code{prior.type} is automatically set as \code{"invwishart"} to yield an inverse-Wishart prior for the variance-covariance matrix of random effects. See "Details".
}
  \item{a, b}{
  positive numbers, specifying the shape and scale parameters of inverse gamma priors for variance(s) of random effects if using \code{prior.type} as \code{"invgamma"} for model \code{"hom_eqcor"}, \code{"het_eqcor"}, \code{"hom_ind"}, or \code{"het_ind"}. The defaults for both parameters are 0.001.
}
  \item{c}{
  positive number, specifying the upper bound of uniform prior for standard deviations of random effects if using \code{prior.type} as \code{"unif"} for model \code{"hom_eqcor"}, \code{"het_eqcor"}, \code{"hom_ind"}, or \code{"het_ind"}. The default is 10.
}
  \item{higher.better}{
  an optional logical value which needs to be specified when estimating the treatment rank probabilities (i.e., "rank.prob" is included in argument \code{param}). \code{TRUE} indicates higher treatment-specific rate implying better treatment, and vice versa. The default is \code{FALSE}.
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
  a vector of character strings of effect sizes. The character strings should be selected from those specified in \code{param} (except \code{"rank.prob"}), and trace plots of the specified effect sizes will be saved in users' current working directory. The default is not drawing trace plots (an empty character string).
}
  \item{dic}{
  a logical value indicating whether the deviance information criterion (DIC) to be calculated. The default is \code{FALSE}.
}
  \item{postdens}{
  a logical value indicating whether to draw the posterior density plots for treatment-specific rates. If \code{TRUE}, a .pdf file containing the density plot will be saved in users' current working directory. The default is \code{FALSE}.
}
  \item{mcmc.samples}{
  a logical value indicating whether to save MCMC posterior samples in the output object. The default is \code{FALSE}.
}
}
\details{
Suppose that a network meta-analysis collects \eqn{I} studies on \eqn{K} treatments, where each study investigates a subset of the \eqn{K} treatments. The outcome is binary, and the follow-up time for each study is reported. Label the studies from \eqn{i = 1} to \eqn{I} and the treatments from \eqn{k = 1} to \eqn{K}. Let \eqn{T_{i}} be the subset of the \eqn{K} treatments that is compared in the \eqn{i}th study. Also, in the \eqn{i}th study, let \eqn{y_{ik}} and \eqn{n_{ik}} be the number of events and the total number of participants in treatment group \eqn{k}. Denote \eqn{f_{i}} as the follow-up time of the \eqn{i}th study. The arm-based network meta-analysis model for these settings is constructed as:
\deqn{y_{ik} \sim Bin (n_{ik}, p_{ik}) \qquad k \in T_{i}}
\deqn{cloglog (p_{ik}) = \log (f_{i}) + \log (\lambda_{ik})}
\deqn{\log (\lambda_{ik}) = \mu_{k} + \nu_{ik}}
\deqn{(\nu_{i1}, \nu_{i2}, \ldots, \nu_{iK})^{T} \sim N (\boldsymbol{0}, \mathbf{\Sigma}_{K}),}
where \eqn{cloglot (t) = \log (-\log(1 - t))} is the complementary log-log link function, and \eqn{\mathbf{\Sigma}_{K}} is a \eqn{K \times K} positive definite covariance matrix. \eqn{\mu_{k}}'s are treatment-specific fixed effects, and random effects \eqn{\nu_{ik}}'s are correlated within each study with the covariance matrix \eqn{\mathbf{\Sigma}_{K}}.

Using an inverse-Wishart prior for \eqn{\mathbf{\Sigma}_{K}}, the above model corresponds to \code{model} = \code{"het_cor"}. Denote \eqn{\sigma_{k}} as the standard deviation of \eqn{\nu_{ik}} and \eqn{\mathbf{D} = diag(\sigma_{1}, \ldots, \sigma_{K})}, then the correlation matrix \eqn{\mathbf{R}_{K} = \mathbf{D}^{-1} \mathbf{\Sigma}_{K} \mathbf{D}^{-1}}.  If we assume that all of the off-diagonal elements in \eqn{\mathbf{R}_{K}} are equal, say to \eqn{\rho}, then this model corresponds to \code{model} = \code{"het_eqcor"}. If we further assume the homogeneity of variances of the random effects, that is, \eqn{\sigma_{k} = \sigma} for \eqn{k = 1, 2, \ldots, K}, then the model is \code{"hom_eqcor"}. Moreover, if \eqn{\mathbf{R}_{K}} is set to be the \eqn{K \times K} identity matrix, then the model becomes \code{"het_ind"} for heterogeneous variances of random effects, and \code{"hom_ind"} for homogeneous variance. Note that the models \code{"het_ind"} and \code{"hom_ind"} lead to independence between different treatments. Therefore, they may break study randomization, and we do not recommend these two models. In addition, for the models \code{"hom_eqcor"}, \code{"het_eqcor"}, \code{"hom_ind"} and \code{"het_ind"}, setting \code{prior.type} as \code{"invgamma"} implies using inverse-gamma priors with shape and scale parameters \eqn{a, b} for \eqn{\sigma_{k}^2} or \eqn{\sigma^2}, and \code{"unif"} implies uniform priors \eqn{U(0, c)} for \eqn{\sigma_{k}} or \eqn{\sigma}.
}
\value{
\code{nma.ab.followup} returns a list with estimates of effect sizes specified in \code{param}. If the argument \code{dic} = \code{TRUE}, the deviance information criterion (DIC) statistic will be returned in the output list. In addition, if \code{conv.diag} = \code{TRUE}, a .txt file containing the point estimates of the potential scale reduction factor and their upper confidence limits by Gelman and Rubin (1992) will be saved in users' current working directory. If \code{postdens} = \code{TRUE}, the posterior densities of treatment-specific absolute risks will be saved as a .pdf file. If \code{trace} is specified, the trace plots are saved as .png files.
}
\references{
Dias S, Sutton AJ, Ades AE, and Welton NJ (2013). "Evidence synthesis for decision making 2: A generalized linear modeling framework for pairwise and network meta-analysis of randomized controlled trials." \emph{Med Decis Making} \bold{33}(5), 607--617.

Gelman A and Rubin, DB (1992). "Inference from iterative simulation using multiple sequences." \emph{Stat Sci} \bold{7}(4), 457--472.

Lu G and Ades AE (2004). "Combination of direct and indirect evidence in mixed treatment comparisons." \emph{Stat Med} \bold{23}(20), 3105--24.

Spiegelhalter DJ, Best NG, Carlin BP, and Van Der Linde A (2002). "Bayesian measures of model complexity and fit." \emph{J R Stat Soc Series B Stat Methodol} \bold{64}(4), 583--639.

Zhang J, Carlin BP, Neaton JD, Soon GG, Nie L, Kane R, Virnig BA, Chu H (2014). "Network meta-analysis of randomized clinical trials: Reporting the proper summaries." \emph{Clin Trials} \bold{11}(2), 246--62.
}
\seealso{
\code{\link{nma.ab}}, \code{\link{nma.ab.cont}}, \code{\link{nma.ab.py}}
}
\examples{
data(diabetes)

# increase n.iter to reach convergence of MCMC
followup.out <- nma.ab.followup(s.id, t.id, r, n, folup, data = diabetes,
  model = "het_cor", n.adapt = 1000, n.iter = 1000)
}
