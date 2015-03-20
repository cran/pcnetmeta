\name{parkinson}
\alias{parkinson}
\docType{data}
\title{
Network Meta-Ananlysis on Parkinson's Disease
}
\description{
An example of network meta-analysis on studies with continuous outcomes.
}
\usage{data("parkinson")}
\format{
  A frame containing 7 studies which compare 5 treatments. The outcomes are the off-time reduction in patients given dopamine agonists as adjunct therapy in Parkinson's disease.
  \describe{
    \item{\code{sid}}{a numeric vector indicating study IDs.}
    \item{\code{tid}}{a numeric vector indicating treatment IDs.}
    \item{\code{mean}}{a numeric vector indicating the mean of continuous outcomes in each treatment group in each study.}
    \item{\code{sd}}{a numeric vector indicating the standard deviance of continuous outcomes in each treatment group in each study.}
    \item{\code{n}}{a numeric vector indicating the total number of participants in each treatment group in each study.}
  }
}
\details{
Treatment 1 is placebo and treatments 2 to 5 are active drugs.
}
\source{
Dias S, Welton NJ, Sutton AJ, and Ades AE (2011). "NICE DSU Technical Support Document 2: a generalised linear modelling framework for pairwise and network meta-analysis of randomised controlled trials." \emph{National Institute for Health and Clinical Excellence}, London, UK.
}
\keyword{datasets}
