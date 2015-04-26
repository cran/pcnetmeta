\name{forest.plot}
\alias{forest.plot}
\title{
Forest Plot of Relative Effect Sizes
}
\description{
\code{forest.plot} generates forest plot, which shows 95\% credible intervals (CIs) for relative effect sizes.
}
\usage{
forest.plot(nma.obj, effect.size, reference, left.margin = 0.5, width = 10, height)
}
\arguments{
  \item{nma.obj}{
  a list object obtained by function \code{\link{nma.ab}}, \code{\link{nma.ab.cont}}, \code{\link{nma.ab.py}}, or \code{\link{nma.ab.followup}}.
}
  \item{effect.size}{
  a character string indicating the relative effect size to be shown in the forest plot. If \code{nma.obj} is obtained from \code{\link{nma.ab}}, this argument can be \code{"OR"} (default), \code{"LOR"}, \code{"RR"}, \code{"LRR"}, or \code{"RD"}. If \code{nma.obj} is obtained from \code{\link{nma.ab.cont}}, this argument should be \code{"diff"} (default). If \code{nma.obj} is obtained from \code{\link{nma.ab.py}} or \code{\link{nma.ab.followup}}, this argument can be \code{"ratio"} (default) or \code{"logratio"}. Note that the specified effect sizes (or its logarithm/exponential) must have been estimated in \code{nma.obj}.
}
  \item{reference}{
  a character string indicating the reference treatment name to be compared against.
}
  \item{left.margin}{
  a positive number to adjust the left margin of the plot. A higher value leads to a wider left margin. The default is 0.5.
}
  \item{width}{
  a positive number indicating the plot width. The default is 10.
}
  \item{height}{
  a positive number indicating the plot height. The default is the treatment number minus 1.
}
}
\value{
A forest plot for relative effect sizes is saved as a .pdf file in users' current working directory.
}

\examples{
data(smoke)
# increase n.iter to reach convergence
nma.out <- nma.ab(s.id, t.id, r, n, data = smoke,
  trtname = c("NC", "SH", "IC", "GC"), param= "LOR",
  model = "het_cor", n.adapt = 1000, n.iter = 100)
forest.plot(nma.out)
# use getwd() to find the path of generated forest plot
}
