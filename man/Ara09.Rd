\name{Ara09}
\alias{Ara09}
\docType{data}
\title{
Network Ara 2009
}
\description{
\code{Ara09} is a network data set served as an example for function \code{\link{nma.networkplot}} and the main function \code{\link{nma.ab}} in this \pkg{pcnetmeta} package.
}
\usage{data(Ara09)}
\format{
  A data frame with 24 observations on the following 5 variables.
  \describe{
    \item{\code{Study.Name}}{a vector of character string indicating study names.}
    \item{\code{Study.ID}}{a numeric vector indicating study IDs.}
    \item{\code{Treatment}}{a numeric vector indicating treatment IDs.}
    \item{\code{r}}{a numeric vector indicating event number for corresponding treatment in a particular study.}
    \item{\code{n}}{a numeric vector indicating total number of participants for corresponding treatment in a particular study.}
  }
}
\details{
The five treatment IDs stand for 1: Atorvastatin 80 mg/day; 2: Simvastatin 40 mg/day; 3: Placebo; 4: Simvastatin 80 mg/day; 5: Rosuvastatin 40 mg/day.
}
\source{
Ara R, Pandor A, Stevens J, Rees A, Rafia R. "Early high-dose lipid-lowering therapy to avoid cardiac events: a systematic review and economic evaluation." \emph{Health Technol Assess} 2009; \bold{13}(34):1-74, 75-118.
}
