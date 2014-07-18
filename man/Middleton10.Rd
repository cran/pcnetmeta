\name{Middleton10}
\alias{Middleton10}
\docType{data}
\title{
Network meta-analysis: Middleton 2010
}
\description{
\code{Middleton10} is a network data set served as an example for the function \code{\link{nma.networkplot}} and the main function \code{\link{nma.ab}} in this \pkg{pcnetmeta} package.
}
\usage{data("Middleton10")}
\format{
  A data frame containing 20 randomized controlled trials on 4 treatments.
  \describe{
    \item{\code{sid}}{a numeric vector indicating study IDs.}
    \item{\code{tid}}{a numeric vector indicating treatment IDs.}
    \item{\code{r}}{a numeric vector indicating event number for a certain treatment in the corresponding study.}
    \item{\code{n}}{a numeric vector indicating total number of participants for a certain treatment in the corresponding study.}
  }
}
\details{
Treatment IDs stand for 1: First generation; 2: Hysterectomy; 3: Second generation; 4: Mirena.
}
\source{
Middleton LJ et al. (2010) "Hysterectomy, endometrial destruction, and levonorgestrel releasing intrauterine system (Mirena) for heavy menstrual bleeding: systematic review and meta-analysis of data from individual patients." \emph{BMJ} \bold{341}, c3929.
}
\keyword{datasets}
