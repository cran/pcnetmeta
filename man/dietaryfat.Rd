\name{dietaryfat}
\alias{dietaryfat}
\docType{data}
\title{
Network Meta-Ananlysis on Dietary Fat
}
\description{
An example of network meta-analysis for binary outcomes with exposure times in person-years reported.
}
\usage{data("dietaryfat")}
\format{
  A data frame containing 10 studies which compare 2 treatments. The total number of person-years at risk are reported for each treatment group in each study.
  \describe{
    \item{\code{s.id}}{a numeric vector indicating study IDs.}
    \item{\code{t.id}}{a numeric vector indicating treatment IDs.}
    \item{\code{py}}{a numeric vector indicating the total person-years in each treatment group in each study.}
    \item{\code{r}}{a numeric vector indicating the total mortality in each treatment group in each study.}
    \item{\code{n}}{a numeric vector indicating the total number of participants in each treatment group in each study.}
  }
}
\details{
This network meta-analysis is collected by Hooper et al (2000) to assess the effects of change in dietary fats on cardiovascular mortality. Treatment 1 is a control diet and treatment 2 is a reduced fat diet. The original study 2 compares three treatments, which include treatment 1 and 2, and another diet. Here, we treat the two different types of diet as the same treatment, but keep the treatment arms separate.
}
\source{
Hooper L, Summerbell CD, Higgins JPT, Thompson RL, Clements G, Capps N, Davey Smith G, Riemersma R, and Ebrahim S (2000). "Reduced or modified dietary fat for preventing cardiovascular disease." \emph{The Cochrane Database of Systematic Reviews}.
}
\keyword{datasets}
