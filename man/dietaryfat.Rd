\name{dietaryfat}
\alias{dietaryfat}
\docType{data}
\title{
Network Meta-Ananlysis on Dietary Fat
}
\description{
An example of network meta-analysis for binary outcomes with exposure times in person-years.
}
\usage{data("dietaryfat")}
\format{
  A frame containing 10 studies which compare 2 treatments. The outcome is mortality, and the total number of person-years at risk are available for each treatment group in each study.
  \describe{
    \item{\code{sid}}{a numeric vector indicating study IDs.}
    \item{\code{tid}}{a numeric vector indicating treatment IDs.}
    \item{\code{py}}{a numeric vector indicating the total person-years in each treatment group in each study.}
    \item{\code{y}}{a numeric vector indicating the total mortality in each treatment group in each study.}
    \item{\code{n}}{a numeric vector indicating the total number of participants in each treatment group in each study.}
  }
}
\details{
Treatment 1 is control and treatment 2 in a diet. The original study 2 compares three treatments, which include treatment 1 and 2, and another diet. Here, we treat the two different types of diet as the same treatment, but keep the treatment arm separately. Thus, the random effects help to provide two correlated estimates of treatment 2 in this study.
}
\source{
Dias S, Welton NJ, Sutton AJ, and Ades AE (2011). "NICE DSU Technical Support Document 2: a generalised linear modelling framework for pairwise and network meta-analysis of randomised controlled trials." \emph{National Institute for Health and Clinical Excellence}, London, UK.
}
\keyword{datasets}
