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
    \item{\code{s.id}}{a numeric vector indicating study IDs.}
    \item{\code{t.id}}{a numeric vector indicating treatment IDs.}
    \item{\code{r}}{a numeric vector indicating event number for a certain treatment in the corresponding study.}
    \item{\code{n}}{a numeric vector indicating the total number of participants for a certain treatment in the corresponding study.}
  }
}
\details{
Treatment IDs stand for 1) no contact; 2) self-help; 3) individual counselling; and 4) group counselling.
}
\source{
Hasselblad V (1998) "Meta-analysis of multitreatment studies." \emph{Med Decis Making} \bold{18}(1), 37--43.

Lu G and Ades AE (2006) "Assessing evidence inconsistency in mixed treatment comparisons." \emph{J Am Stat Assoc} \bold{101}(474), 447--59.
}
\keyword{datasets}
