\name{absolute.plot}
\alias{absolute.plot}
\title{
Plotting Treatment-Specific 95\% Credible Intervals
}
\description{
Generates a plot of 95\% credible intervals (CIs) for the treatment-specific effect sizes.
}
\usage{
absolute.plot(nma.obj, alphabetic = TRUE, digits = 2,
              save = TRUE, width = 5, height, network.name)
}
\arguments{
  \item{nma.obj}{
  a list object obtained by function \code{\link{nma.ab.bin}}, \code{\link{nma.ab.cont}}, \code{\link{nma.ab.py}}, or \code{\link{nma.ab.followup}}.
}
  \item{alphabetic}{
  a logical value indicating whether to sort the CIs according to treatment names alphabetically. The default is \code{TRUE}. If \code{FALSE}, the CIs are plotted in the order of treatment IDs specified in \code{nma.obj}.
}
  \item{digits}{
  an integer specifying the number of decimal places to be used for the point estimates and 95\% CIs. The default is 2.
}
  \item{save}{
  a logical value indicating whether to save the plot as a .pdf file. If \code{FALSE}, a plot window will be displayed. The following arguments are only used when \code{save} = \code{TRUE}.
}
  \item{width}{
  a positive number specifying the plot width. The default is 5.
}
  \item{height}{
  a positive number specifying the plot height. The default is the treatment number minus 1.
}
  \item{network.name}{
  a character string indicating the network name to be used for the produced .pdf file name.
}
}
\value{
None.
}

\examples{
data("smoke")
# increase n.iter to reach convergence
set.seed(1234)
nma.out <- nma.ab.bin(s.id, t.id, r, n, data = smoke,
  trtname = c("NC", "SH", "IC", "GC"), param= "AR",
  model = "het_cor", n.adapt = 1000, n.iter = 100, n.chains = 1)
absolute.plot(nma.out, save = FALSE)
#absolute.plot(nma.out)
absolute.plot(nma.out, alphabetic = FALSE, save = FALSE)
}
\keyword{plot}