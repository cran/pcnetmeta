\name{pcnetmeta-package}
\alias{pcnetmeta-package}
\alias{pcnetmeta}
\docType{package}
\title{
Methods for patient-centered network meta-analysis
}
\description{
Provides functions to conduct network meta-analysis using arm-based method, which was proposed by Zhang et al (2014); contains functions to generate summary table for effect sizes (e.g., odds ratio, population-averaged event rate), plot network comparisons, and plot credible interval for population-averaged event rate.
}
\details{
\tabular{ll}{
Package: \tab pcnetmeta\cr
Type: \tab Package\cr
Version: \tab 1.2\cr
Date: \tab 2014-07-18\cr
Depends: \tab R (>= 3.0.2), R2jags, network, runjags\cr
License: \tab GPL (>=2)\cr
}
Currently, most popular approaches in network meta-analysis are contrast-based, and they focus on modeling the relative treatment effects, such as odds ratio. However, the arm-based method, proposed by Zhang et al (2014), considers modeling from the perspective of missing data analysis, and it focuses on estimating the absolute risk for each treatment. This package conducts network meta-analysis by the arm-based method, and \code{\link{nma.ab}} is the most important function, which generates summary result file containing estimated absolute risk (AR), relative risk (RR), risk difference (RD), odds ratio (OR), etc.

Analysis in this package is conducted through JAGS. Note that this package does not include a copy of JAGS library, so you must install this software separately. Please refer to the JAGS home page at \url{http://mcmc-jags.sourceforge.net/} for instructions on downloading JAGS.
}
\author{
Lifeng Lin, Jing Zhang <zhan0979@umn.edu>, and Haitao Chu <chux0051@umn.edu>

Maintainer: Lifeng Lin <linl@umn.edu>
}
\references{
Zhang J, Carlin BP, Neaton JD, Soon GG, Nie L, Kane R, Virnig BA, Chu H (2014). "Network meta-analysis of randomized clinical trials: Reporting the proper summaries." \emph{Clin Trials} \bold{11}(2), 246--262.

Lu G, Ades AE (2004). "Combination of direct and indirect evidence in mixed treatment comparisons." \emph{Stat Med} \bold{23}(20), 3105--24.

Butts CT (2008). "network: A Package for Managing Relational Data in R." \emph{J Stat Softw} \bold{24}(2), 1--36.

Middleton LJ et al. (2010) "Hysterectomy, endometrial destruction, and levonorgestrel releasing intrauterine system (Mirena) for heavy menstrual bleeding: systematic review and meta-analysis of data from individual patients." \emph{BMJ} \bold{341}, c3929.
}
