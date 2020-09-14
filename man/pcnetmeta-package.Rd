\name{pcnetmeta-package}
\alias{pcnetmeta-package}
\alias{pcnetmeta}
\docType{package}
\title{
Patient-Centered Network Meta-Analysis
}
\description{
Provides functions to perform Bayesian arm-based network meta-analysis for datasets with binary, continuous, and count outcomes.
}
\details{
Much effort in the current literature of network meta-analysis has been devoted to contrast-based methods, which focus on modeling the relative treatment effects, such as log odds ratios. In some situations, patient-centered treatment-specific parameters, such as absolute risks, may be preferred. For example, when the outcome is binary, consider two scenarios comparing treatments A and B according to one-year survival rates: (i) \eqn{\pi_{A} = 0.8} vs.{\sspace}\eqn{\pi_{B} = 0.5} and (ii) \eqn{\pi_{A} = 0.004} vs.{\sspace}\eqn{\pi_{B} = 0.001}. Both scenarios yield an odds ratio of 4.0, but patients would prefer treatment A in scenario (i) more strongly than in scenario (ii). The contrast-based network meta-analysis requires external data sources or separate modeling to estimate treatment-specific parameters.

Alternatively, the arm-based network meta-analysis focuses on estimating treatment-specific parameters, and relative effects can be subsequently obtained. The arm-based models for binary outcomes are discussed in Salanti et al.{\sspace}(2008) and Zhang et al.{\sspace}(2014).

This package provides user-friendly functions for the arm-based network meta-analysis. The function \code{\link{nma.ab.bin}} uses the model in Zhang et al.{\sspace}(2014) for binary outcomes; it provides estimates for the absolute risk (AR), risk ratio (RR), risk difference (RD), odds ratio (OR), log RR, log OR, as well as rank probabilities. This package also handles continuous outcomes and count data. A detailed instruction of using this package is provided in Lin et al.{\sspace}(2017).

Parameter estimation in the arm-based Bayesian hierarchical models is performed through JAGS. Note that this package does not include a copy of JAGS library, so users must install JAGS separately. Please refer to the JAGS home page at \url{http://mcmc-jags.sourceforge.net/} for instructions on downloading and installing JAGS.
}
\author{
Lifeng Lin, Jing Zhang, and Haitao Chu

Maintainer: Lifeng Lin <linl@stat.fsu.edu>
}
\references{
Dias S, Sutton AJ, Ades AE, Welton NJ (2013). "Evidence synthesis for decision making 2: a generalized linear modeling framework for pairwise and network meta-analysis of randomized controlled trials." \emph{Medical Decision Making}, \bold{33}(5), 607--617. <\doi{10.1177/0272989X12458724}>

Lin L, Zhang J, Hodges JS, Chu H (2017). "Performing arm-based network meta-analysis in R with the pcnetmeta package." \emph{Journal of Statistical Software}, \bold{80}(5), 1--25. <\doi{10.18637/jss.v080.i05}>

Lu G, Ades AE (2004). "Combination of direct and indirect evidence in mixed treatment comparisons." \emph{Statistics in Medicine}, \bold{23}(20), 3105--3124. <\doi{10.1002/sim.1875}>

Salanti G, Higgins JPT, Ades AE, Ioannidis JPA (2008). "Evaluation of networks of randomized trials." \emph{Statistical Methods in Medical Research} \bold{17}(3), 279--301. <\doi{10.1177/0962280207080643}>

Zhang J, Carlin BP, Neaton JD, Soon GG, Nie L, Kane R, Virnig BA, Chu H (2014). "Network meta-analysis of randomized clinical trials: Reporting the proper summaries." \emph{Clinical Trials}, \bold{11}(2), 246--262. <\doi{10.1177/1740774513498322}>
}
