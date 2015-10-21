\name{absolute.plot}
\alias{absolute.plot}
\title{
Plotting Treatment-Specific 95\% Credible Intervals
}
\description{
\code{absolute.plot} generates a plot of 95\% credible intervals (CIs) for the treatment-specific effect sizes.
}
\usage{
absolute.plot(nma.obj, alphabetic = TRUE, digits = 2, width = 5, height)
}
\arguments{
  \item{nma.obj}{
  a list object obtained by function \code{\link{nma.ab}}, \code{\link{nma.ab.cont}}, \code{\link{nma.ab.py}}, or \code{\link{nma.ab.followup}}.
}
  \item{alphabetic}{
  a logical value indicating whether to sort the CIs according to treatment names alphabetically. The default is \code{TRUE}. If \code{FALSE}, the CIs are plotted in the order of treatment IDs specified in \code{nma.obj}.
}
  \item{digits}{
  an integer indicating the number of decimal places to be used for the point estimates and 95\% CIs. The default is 2.
}
  \item{width}{
  a positive number indicating the plot width. The default is 5.
}
  \item{height}{
  a positive number indicating the plot height. The default is the treatment number minus 1.
}
}

\examples{
data(smoke)
# increase n.iter to reach convergence
# increase n.adapt to enhance efficiency
set.seed(1234)
nma.out <- nma.ab(s.id, t.id, r, n, data = smoke,
  trtname = c("NC", "SH", "IC", "GC"), param= "AR",
  model = "het_cor", n.adapt = 400, n.iter = 100, n.chains = 1)
absolute.plot(nma.out)
absolute.plot(nma.out, alphabetic = FALSE)
}