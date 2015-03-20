\name{pcnetmeta-package}
\alias{pcnetmeta-package}
\alias{pcnetmeta}
\docType{package}
\title{
Methods for Patient-Centered Network Meta-Analysis
}
\description{
Provides functions to conduct network meta-analysis using arm-based methods for different data types, which include binary and continuous outcomes.
}
\details{
\tabular{ll}{
Package: \tab pcnetmeta\cr
Type: \tab Package\cr
Version: \tab 2.0\cr
Date: \tab 2015-03-20\cr
Depends: \tab R (>= 3.1.3), rjags (>= 3-14)\cr
Imports: \tab coda\cr
License: \tab GPL (>= 2)\cr
}
Currently, much effort in network meta-analysis has been devoted to contrast-based methods, which focus on modeling the relative treatment effects, such as odds ratio. However, in some situations, patient-centered parameters, such as absolute risk, are preferred. For example, when the outcome is binary and there are two scenarios of survival rates on treatment A and B: (i) \eqn{\pi_{A} = 0.8} vs \eqn{\pi_{B} = 0.5} and (ii) \eqn{\pi_{A} = 0.004} vs \eqn{\pi_{B} = 0.001}. Both scenarios yield an odds ratio 4.0, but clearly patients would prefer scenario (i) as its survival rates are much higher. The contrast-based network meta-analysis models require external data sources or separate models to estimate the patient-centered parameters. However, the external data sources may not be representative for some populations.

Alternatively, arm-based method in network meta-analysis focuses on estimating patient-centerd arm parameters, and parameters for relative effects can be subsequently obtained. The arm-based models for binary outcomes are discussed in articles such as Salanti et al (2008) and Zhang et al (2014).

This package provides user-friendly functions for the arm-based network meta-analysis. Function \code{\link{nma.ab}} follows the models in Zhang et al (2014) for binary outcomes, and it can provide estimations for absolute risk (AR), relative risk (RR), risk difference (RD), odds ratio (OR), log RR, log OR, and probability of being the best treatment. Furthermore, this package extends the arm-based models for other data types, such as continuous outcomes. Specifically, \code{\link{nma.ab.cont}} deals with network meta-analysis on studies with continuous outcomes, \code{\link{nma.ab.py}} uses a Poisson likelihood when the exposure times in pearson-years at risk are available, and \code{\link{nma.ab.followup}} considers different follow-up times for different studies. This package also includes several illustrative datasets, which mostly come from Dias et al (2011).

Parameter estimations in the Bayesian hierarchical models are conducted through JAGS. Note that this package does not include a copy of JAGS library, so you must install this software separately. Please refer to the JAGS home page at \url{http://mcmc-jags.sourceforge.net/} for instructions on downloading JAGS.
}
\author{
Lifeng Lin, Jing Zhang, and Haitao Chu

Maintainer: Lifeng Lin <linl@umn.edu>
}
\references{
Dias S, Welton NJ, Sutton AJ, and Ades AE (2011). "NICE DSU Technical Support Document 2: a generalised linear modelling framework for pairwise and network meta-analysis of randomised controlled trials." \emph{National Institute for Health and Clinical Excellence}, London, UK.

Lu G, Ades AE (2004). "Combination of direct and indirect evidence in mixed treatment comparisons." \emph{Stat Med} \bold{23}(20), 3105--24.

Salanti G, Higgins JP, Ades AE, and Ioannidis JP (2008). "Evaluation of networks of randomized trials." \emph{Statistical methods in medical research} \bold{17}(3), 279--301.

Zhang J, Carlin BP, Neaton JD, Soon GG, Nie L, Kane R, Virnig BA, Chu H (2014). "Network meta-analysis of randomized clinical trials: Reporting the proper summaries." \emph{Clin Trials} \bold{11}(2), 246--262.
}
