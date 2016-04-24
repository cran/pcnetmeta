\name{diabetes}
\alias{diabetes}
\docType{data}
\title{
Network Meta-Ananlysis on Diabetes
}
\description{
An example of network meta-analysis for binary outcomes with follow-up times reported.
}
\usage{data("diabetes")}
\format{
  A data frame containing 22 studies which compare 6 treatments. The outcome is binary; mean follow-up time in each study is collected.
  \describe{
    \item{\code{s.id}}{a numeric vector indicating study IDs.}
    \item{\code{folup}}{a numeric vector indicating follow-up time (in years) for each study.}
    \item{\code{t.id}}{a numeric vector indicating treatment IDs.}
    \item{\code{r}}{a numeric vector indicating the total mortality in each treatment group in each study.}
    \item{\code{n}}{a numeric vector indicating the total number of participants in each treatment group in each study.}
  }
}
\details{
This network meta-analysis is studied by Elliott and Meyer (2007) to assess the effects of anti-hypertensive agents on incident diabetes. Treatment IDs represent 1) diuretic; 2) placebo; 3) \eqn{\beta}-blocker; 4) CCB; 5) ACE inhibitor; and 6) ARB.
}
\source{
Elliott WJ and Meyer PM (2007). "Incident diabetes in clinical trials of antihypertensive drugs: a network meta-analysis." \emph{Lancet} \bold{369}(9557), 201--7.
}
\keyword{datasets}
