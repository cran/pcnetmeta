\name{ci.plot}
\alias{ci.plot}
\title{
Plot credible interval
}
\description{
\code{ci.plot} creates credible interval plot for estimated event rates from the summary result file given by function \code{\link{nma.ab}}.
}
\usage{
ci.plot(summary.stat, trtname, graphtitle = "")
}
\arguments{
  \item{summary.stat}{
  a data frame read from the summary result file by function \code{\link{nma.ab}}. We can use function \code{read.table} to import the summary file.
}
  \item{trtname}{
  a vector of character string indicating treatment names. These will be shown on x-axis in the plot. The default is "trt1", "trt2", and so on.
}
  \item{graphtitle}{
  a charater string indicating the graph title. The default is an empty character.
}
}
\value{
A credible interval plot for estimated event rate is created. Different treatments are placed on x-axis, while the y-axis shows the lower bound, upper bound and median of credible interval.
}
\references{
Zhang J, Carlin BP, Neaton JD, Soon GG, Nie L, Kane R, Virnig BA and Chu H. "Network meta-analysis of randomized clinical trials: Reporting the proper summaries." \emph{Clinical Trials} 2013.
}
\author{
Lifeng Lin, Jing Zhang, and Haitao Chu.
}

\examples{
## CI plot for network Ara 2009
data(Ara09)
attach(Ara09)
nma.ab(s.id = Study.ID, t.id = Treatment, event.n = r, total.n = n,
       f.name = "Ara09_", n.iter = 2000, dic = FALSE)
detach(Ara09)
ci.plot(summary.stat = read.table("Ara09_Summary.stat", header = TRUE),
        graphtitle = "CI plot of estimated event rate for Ara09 network")

## CI plot for network Lam 2007
data(Lam07.Summary)
data(Lam07)
ci.plot(summary.stat = Lam07.Summary, trtname = Lam07$trtname,
        graphtitle = "CI plot of estimated event rate for Lam07 network")
}

