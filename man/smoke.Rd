\name{smoke}
\alias{smoke}
\docType{data}
\title{
Network Meta-Analysis on Smoking Cessation Data
}
\description{
An example of network meta-analysis on studies with binary outcomes.
}
\usage{data("smoke")}
\format{
  A frame containing 24 studies on smoking cessation, comparing 4 treatments: no contact, self-help, individual counselling, and group counselling.
  \describe{
    \item{\code{sid}}{a numeric vector indicating study IDs.}
    \item{\code{tid}}{a numeric vector indicating treatment IDs.}
    \item{\code{r}}{a numeric vector indicating event number for a certain treatment in the corresponding study.}
    \item{\code{n}}{a numeric vector indicating the total number of participants for a certain treatment in the corresponding study.}
  }
}
\details{
Treatment IDs stand for 1: no contact; 2: self-help; 3: individual counselling; 4: group counselling.
}
\source{
Hasselblad V (1998) "Meta-analysis of multitreatment studies." \emph{Med Decis Making} 18, 37--43

Lu G and Ades AE (2006) "Assessing evidence inconsistency in mixed treatment comparisons." \emph{JASA} \bold{101}(474), 447--59.
}
\keyword{datasets}
