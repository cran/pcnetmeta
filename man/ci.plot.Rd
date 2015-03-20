\name{ci.plot}
\alias{ci.plot}
\title{
Plot Credible Intervals with Medians
}
\description{
\code{ci.plot} generates a plot of 95% credible intervals (CIs) with medians for the arm parameters from objects obtained by \code{\link{nma.ab}}, \code{\link{nma.ab.cont}}, \code{\link{nma.ab.py}}, \code{\link{nma.ab.followup}}.
}
\usage{
ci.plot(nma.obj, alphabetic = TRUE)
}
\arguments{
  \item{nma.obj}{
  a list object obtained by \code{\link{nma.ab}}, \code{\link{nma.ab.cont}}, \code{\link{nma.ab.py}}, or \code{\link{nma.ab.followup}}.
}
  \item{alphabetic}{
  a logical value indicating whether to sort the CIs according to treatment names alphabetically. The default is \code{TRUE}. If \code{FALSE}, the CIs are plotted in the order of treatment IDs specified in \code{\link{nma.ab}}, \code{\link{nma.ab.cont}}, \code{\link{nma.ab.py}}, or \code{\link{nma.ab.followup}}.
}
}

\examples{
## CI plot for network Lu and Ades 2006
data(smoke)
# increase n.iter to reach convergence of MCMC
nma.out <- nma.ab(smoke$sid, smoke$tid, smoke$r, smoke$n,
  model = "het_cor", trtname = c("No contact", "Self-help",
  "Individual counselling", "Group counselling"), param=c("AR"),
  n.adapt = 1000, n.iter = 100)
ci.plot(nma.out)
ci.plot(nma.out, FALSE)
}