\name{contrast.plot}
\alias{contrast.plot}
\title{
Contrast Plot of Relative Effect Measures
}
\description{
Generates a plot that shows 95\% credible intervals (CIs) for relative effect measures.
}
\usage{
contrast.plot(nma.obj, effect.measure, reference, digits = 2,
              save = TRUE, width = 5, height, network.name)
}
\arguments{
  \item{nma.obj}{
  a list object obtained by function \code{\link{nma.ab.bin}}, \code{\link{nma.ab.cont}}, \code{\link{nma.ab.py}}, or \code{\link{nma.ab.followup}}.
}
  \item{effect.measure}{
  a character string specifying the relative effect measure to be shown in the contrast plot. If \code{nma.obj} is obtained from \code{\link{nma.ab.bin}}, this argument can be \code{"OR"} (the default), \code{"LOR"}, \code{"RR"}, \code{"LRR"}, or \code{"RD"}; it can be also \code{"OR.med"} or \code{"LOR.med"} if they are specified in \code{\link{nma.ab.bin}} with \code{link} = \code{"logit"}. If \code{nma.obj} is obtained from \code{\link{nma.ab.cont}}, this argument should be \code{"diff"} (the default). If \code{nma.obj} is obtained from \code{\link{nma.ab.py}} or \code{\link{nma.ab.followup}}, this argument can be \code{"ratio"} (the default) or \code{"logratio"}. Note that the specified effect measures (or its logarithm/exponential) must have been estimated in \code{nma.obj}.
}
  \item{reference}{
  a character string specifying the reference treatment name to be compared against.
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
  trtname = c("NC", "SH", "IC", "GC"), param = "LOR",
  model = "het_cor", n.adapt = 1000, n.iter = 100, n.chains = 1)
contrast.plot(nma.out, save = FALSE)
#contrast.plot(nma.out)
}
\keyword{plot}