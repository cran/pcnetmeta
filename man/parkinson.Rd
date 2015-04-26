\name{parkinson}
\alias{parkinson}
\docType{data}
\title{
Network Meta-Ananlysis on Parkinson's Disease
}
\description{
An example of network meta-analysis for continuous outcomes.
}
\usage{data("parkinson")}
\format{
  A data frame containing 7 studies which compare 5 treatments.
  \describe{
    \item{\code{s.id}}{a numeric vector indicating study IDs.}
    \item{\code{t.id}}{a numeric vector indicating treatment IDs.}
    \item{\code{mean}}{a numeric vector indicating the mean of continuous outcomes in each treatment group in each study.}
    \item{\code{sd}}{a numeric vector indicating the standard deviance of continuous outcomes in each treatment group in each study.}
    \item{\code{n}}{a numeric vector indicating the total number of participants in each treatment group in each study.}
  }
}
\details{
The continuous outcome measures the off-time reduction in patients given dopamine agonists as adjunct therapy in Parkinson's disease. Treatment 1 is placebo and Treatments 2 to 5 are active drugs.
}
\source{
Dias S, Sutton AJ, Ades AE, and Welton NJ (2013). "Evidence synthesis for decision making 2: A generalized linear modeling framework for pairwise and network meta-analysis of randomized controlled trials." \emph{Med Decis Making} \bold{33}(5), 607--617.
}
\keyword{datasets}
