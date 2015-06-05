\name{absolute.plot}
\alias{absolute.plot}
\title{
Plotting Treatment-Specific 95\% Credible Intervals
}
\description{
\code{absolute.plot} generates a plot of 95\% credible intervals (CIs) for the treatment-specific effect sizes.
}
\usage{
absolute.plot(nma.obj, alphabetic = TRUE)
}
\arguments{
  \item{nma.obj}{
  a list object obtained by function \code{\link{nma.ab}}, \code{\link{nma.ab.cont}}, \code{\link{nma.ab.py}}, or \code{\link{nma.ab.followup}}.
}
  \item{alphabetic}{
  a logical value indicating whether to sort the CIs according to treatment names alphabetically. The default is \code{TRUE}. If \code{FALSE}, the CIs are plotted in the order of treatment IDs specified in \code{nma.obj}.
}
}

\examples{
data(smoke)
# increase n.iter to reach convergence
# increase n.adapt to enhance efficiency
nma.out <- nma.ab(s.id, t.id, r, n, data = smoke,
  trtname = c("NC", "SH", "IC", "GC"), param= "AR",
  model = "het_cor", n.adapt = 0, n.iter = 100, n.chains = 1)
absolute.plot(nma.out)
absolute.plot(nma.out, FALSE)
}