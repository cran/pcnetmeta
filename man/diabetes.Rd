\name{diabetes}
\alias{diabetes}
\docType{data}
\title{
Network Meta-Ananlysis on Diabetes
}
\description{
An example of network meta-analysis for binary outcomes with follow-up times.
}
\usage{data("diabetes")}
\format{
  A frame containing 10 studies which compare 6 treatments. The outcome is new case of diabetes, and the follow-up time in each study are available.
  \describe{
    \item{\code{sid}}{a numeric vector indicating study IDs.}
    \item{\code{folup}}{a numeric vector indicating follow-up time (in years) for each study.}
    \item{\code{tid}}{a numeric vector indicating treatment IDs.}
    \item{\code{y}}{a numeric vector indicating the total mortality in each treatment group in each study.}
    \item{\code{n}}{a numeric vector indicating the total number of participants in each treatment group in each study.}
  }
}
\details{
Treatment IDs represent 1: diuretic; 2: placebo; 3: \eqn{\beta}-blocker; 4: CCB; 5: ACE inhibitor; 6: ARB.
}
\source{
Dias S, Welton NJ, Sutton AJ, and Ades AE (2011). "NICE DSU Technical Support Document 2: a generalised linear modelling framework for pairwise and network meta-analysis of randomised controlled trials." \emph{National Institute for Health and Clinical Excellence}, London, UK.
}
\keyword{datasets}
