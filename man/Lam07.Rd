\name{Lam07}
\alias{Lam07}
\docType{data}
\title{
Network Lam 2007
}
\description{
\code{Lam07} is a list containing two elements: one is a vector of character string indicating treatment names; another is a data frame for a network data set. \code{Lam07} is served as an example for the function \code{\link{nma.networkplot}} and the main function \code{\link{nma.ab}} in this \pkg{pcnetmeta} package.
}
\usage{data(Lam07)}
\format{
  A list with 2 elements: \code{trtname} and \code{data}. \code{trtname} is a vector of character string, while \code{data} is a data frame with 24 rows and 5 columns.
  \describe{
    \item{\code{trtname}}{a vector of character string indicating five treatment names.}
    \item{\code{data}}{
      \itemize{
      \item{\code{Study.Name}: a vector of character string indicating study names.}
      \item{\code{Study.ID}: a numeric vector indicating study IDs.}
      \item{\code{Treatment}: a numeric vector indicating treatment IDs.}
      \item{\code{r}: a numeric vector indicating event number for a certain treatment in the corresponding study.}
      \item{\code{n}: a numeric vector indicating total number of participants for a certain treatment in the corresponding study.}
      }  
    }
  }
}
\source{
Lam SK, Owen A. "Combined resynchronisation and implantable defibrillator therapy in left ventricular dysfunction: Bayesian network meta-analysis of randomised controlled trials." \emph{BMJ} 2007; \bold{335}(7626):925.
}
