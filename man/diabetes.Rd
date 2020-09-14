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
  A data frame containing 22 studies which compare a total of 6 treatments. The outcome is binary; the mean follow-up time in each study is collected.
  \describe{
    \item{\code{s.id}}{a numeric vector of study IDs.}
    \item{\code{folup}}{a numeric vector of follow-up time (in years) for each study.}
    \item{\code{t.id}}{a numeric vector of treatment IDs.}
    \item{\code{r}}{a numeric vector of the total mortality in each treatment group in each study.}
    \item{\code{n}}{a numeric vector of the total number of participants in each treatment group in each study.}
  }
}
\details{
This network meta-analysis is reported by Elliott and Meyer (2007) to assess the effects of anti-hypertensive agents on incident diabetes. Treatment IDs represent 1) diuretic; 2) placebo; 3) \eqn{\beta}-blocker; 4) calcium-channel blocker (CCB); 5) angiotensin-converting-enzyme (ACE) inhibitor; and 6) angiotensin-receptor blocker (ARB).
}
\source{
Elliott WJ, Meyer PM (2007). "Incident diabetes in clinical trials of antihypertensive drugs: a network meta-analysis." \emph{Lancet}, \bold{369}(9557), 201--207. <\doi{10.1016/S0140-6736(07)60108-1}>
}
\keyword{dataset}
