\name{smoke}
\alias{smoke}
\docType{data}
\title{
Network Meta-Analysis on Smoking Cessation Data
}
\description{
An example of network meta-analysis for binary outcomes.
}
\usage{data("smoke")}
\format{
  A data frame containing 24 studies on smoking cessation, comparing four treatments.
  \describe{
    \item{\code{s.id}}{a numeric vector of study IDs.}
    \item{\code{t.id}}{a numeric vector of treatment IDs.}
    \item{\code{r}}{a numeric vector of the event count for a certain treatment in the corresponding study.}
    \item{\code{n}}{a numeric vector of the total number of participants for a certain treatment in the corresponding study.}
  }
}
\details{
Treatment IDs stand for 1) no contact; 2) self-help; 3) individual counseling; and 4) group counseling.
}
\source{
Hasselblad V (1998) "Meta-analysis of multitreatment studies." \emph{Medical Decision Making}, \bold{18}(1), 37--43. <\doi{10.1177/0272989X9801800110}>

Lu G, Ades AE (2006) "Assessing evidence inconsistency in mixed treatment comparisons." \emph{Journal of the American Statistical Association} \bold{101}(474), 447--459. <\doi{10.1198/016214505000001302}>
}
\keyword{dataset}
