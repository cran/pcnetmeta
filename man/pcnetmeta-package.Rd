\name{pcnetmeta-package}
\alias{pcnetmeta-package}
\alias{pcnetmeta}
\docType{package}
\title{
Patient-Centered Network Meta-Analysis
}
\description{
Provides functions to perform arm-based network meta-analysis for datasets with binary and continuous outcomes.
}
\details{
Currently, much effort in network meta-analysis has been devoted to contrast-based methods, which focus on modeling the relative treatment effects, such as log odds ratio. However, in some situations, patient-centered treatment-specific parameters, such as absolute risk, are preferred. For example, when the outcome is binary, consider two scenarios comparing treatments A and B according to one-year survival rates: (i) \eqn{\pi_{A} = 0.8} vs \eqn{\pi_{B} = 0.5} and (ii) \eqn{\pi_{A} = 0.004} vs \eqn{\pi_{B} = 0.001}. Both scenarios yield an odds ratio 4.0, but patients would prefer treatment A in scenario (i) more strongly than in scenario (ii). The contrast-based network meta-analysis requires external data sources or separate modeling to estimate treatment-specific parameters.

Alternatively, arm-based network meta-analysis focuses on estimating treatment-specific parameters, and relative effects can be subsequently obtained. The arm-based models for binary outcomes are discussed in Salanti et al (2008) and Zhang et al (2014).

This package provides user-friendly functions for arm-based network meta-analysis. Function \code{\link{nma.ab}} uses the model in Zhang et al (2014) for binary outcomes; it provides estimates for absolute risk (AR), relative risk (RR), risk difference (RD), odds ratio (OR), log RR, and log OR. This package also handles continuous outcomes and count data.

Parameter estimation in arm-based Bayesian hierarchical models are performed through JAGS. Note that this package does not include a copy of JAGS library, so users must install JAGS separately. Please refer to the JAGS home page at \url{http://mcmc-jags.sourceforge.net/} for instructions on downloading and installing JAGS.
}
\author{
Lifeng Lin, Jing Zhang, and Haitao Chu

Maintainer: Lifeng Lin <linl@umn.edu>
}
\references{
Dias S, Sutton AJ, Ades AE, and Welton NJ (2013). "Evidence synthesis for decision making 2: A generalized linear modeling framework for pairwise and network meta-analysis of randomized controlled trials." \emph{Med Decis Making} \bold{33}(5), 607--617.

Lu G, Ades AE (2004). "Combination of direct and indirect evidence in mixed treatment comparisons." \emph{Stat Med} \bold{23}(20), 3105--24.

Salanti G, Higgins JP, Ades AE, and Ioannidis JP (2008). "Evaluation of networks of randomized trials." \emph{Stat Methods Med Res} \bold{17}(3), 279--301.

Zhang J, Carlin BP, Neaton JD, Soon GG, Nie L, Kane R, Virnig BA, Chu H (2014). "Network meta-analysis of randomized clinical trials: Reporting the proper summaries." \emph{Clin Trials} \bold{11}(2), 246--262.
}
