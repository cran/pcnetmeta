\name{pcnetmeta-package}
\alias{pcnetmeta-package}
\alias{pcnetmeta}
\docType{package}
\title{
Methods for patient-centered network meta-analysis
}
\description{
Provides functions to do network meta-analysis using arm-based method, which was proposed by Zhang et al (2013); contains functions to generate summary table for the results, plot network comparisons, and plot credible interval for estimated event rate.
}
\details{
\tabular{ll}{
Package: \tab pcnetmeta\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2014-01-17\cr
Depends: \tab R (>= 3.0.2), R2jags, network\cr
License: \tab GPL (>=2)\cr
}
Most popular approaches in network meta-analysis are contrast-based, and they focus on modeling the relative treatment effects, such as ORs. However, the arm-based method, proposed by Zhang et al (2013), treats the analysis from the perspective of missing data, and it focuses on estimating the event rate for each treatment. This package conducts network meta-analysis by the arm-based method, and \code{\link{nma.ab}} is the most important function, which generates summary result file containing estimated event rates, RRs, RDs, ORs, etc.

Analyses in this package are conducted by JAGS. Note that this package does not include a copy of JAGS library, so you must install this software separately. Please refer to the JAGS home page at \url{http://mcmc-jags.sourceforge.net/} for instructions on downloading JAGS.
}
\author{
Lifeng Lin, Jing Zhang <zhan0979@umn.edu>, and Haitao Chu <chux0051@umn.edu>

Maintainer: Lifeng Lin <lifeng.lin@outlook.com>
}
\references{
Zhang J, Carlin BP, Neaton JD, Soon GG, Nie L, Kane R, Virnig BA and Chu H. "Network meta-analysis of randomized clinical trials: Reporting the proper summaries." \emph{Clinical Trials} 2013.
}
