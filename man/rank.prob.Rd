\name{rank.prob}
\alias{rank.prob}
\title{
Plotting Treatment Rank Probabilities
}
\description{
\code{rank.prob} generates a plot for treatment rank probabilities.
}
\usage{
rank.prob(nma.obj)
}
\arguments{
  \item{nma.obj}{
  a list object obtained by function \code{\link{nma.ab}}, \code{\link{nma.ab.cont}}, \code{\link{nma.ab.py}}, or \code{\link{nma.ab.followup}}.
}
}
\details{
A plot is generated. Each vertical bar represents probabilities of being different ranks for a specific treatment. A darker area indicates the probability of being a higher rank. The black area indicates the probability of being the best treatment.
}

\examples{
data(smoke)
# increase n.iter to reach convergence
nma.out <- nma.ab(s.id, t.id, r, n, data = smoke,
  trtname = c("NC", "SH", "IC", "GC"), param= "rank.prob",
  model = "het_cor", higher.better = TRUE, n.adapt = 1000, n.iter = 100)
rank.prob(nma.out)
}
