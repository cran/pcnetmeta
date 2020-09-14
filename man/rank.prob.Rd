\name{rank.prob}
\alias{rank.prob}
\title{
Plotting Treatment Rank Probabilities
}
\description{
\code{rank.prob} generates a plot for treatment rank probabilities.
}
\usage{
rank.prob(nma.obj, cex.axis = 1, cex.lab = 1, median = FALSE)
}
\arguments{
  \item{nma.obj}{
  a list object obtained by function \code{\link{nma.ab.bin}}, \code{\link{nma.ab.cont}}, \code{\link{nma.ab.py}}, or \code{\link{nma.ab.followup}}.
}
  \item{cex.axis}{
  a numeric value specifying the size of the tick label numbers/text.
}
  \item{cex.lab}{
  a numeric value specifying the size of the axis label text.
}
  \item{median}{
  a logical value indicating whether the rank probabilities are based on the marginal effect measures (\code{FALSE}, the default) or the overall median odds ratios. This argument is only used when \code{link} = \code{"logit"} in \code{\link{nma.ab.bin}}.
}
}
\details{
A plot of treatment rank probabilities is generated. Each vertical bar represents probabilities of being different ranks for a specific treatment. A darker area indicates the probability of being a higher rank. The black area indicates the probability of being the best treatment.
}
\value{
None.
}

\examples{
data("smoke")
# increase n.iter to reach convergence
set.seed(1234)
nma.out <- nma.ab.bin(s.id, t.id, r, n, data = smoke,
  trtname = c("NC", "SH", "IC", "GC"), param = "rank.prob", model = "het_cor",
  higher.better = TRUE, n.adapt = 1000, n.iter = 100, n.chains = 1)
rank.prob(nma.out)
}
\keyword{plot}