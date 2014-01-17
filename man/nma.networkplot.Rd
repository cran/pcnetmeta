\name{nma.networkplot}
\alias{nma.networkplot}
\title{
Network plot of treatment comparisons
}
\description{
\code{nma.networkplot} plots a network presenting comparisons of treatments in various studies.
}
\usage{
nma.networkplot(c1, c2, percomparison, trtname, weight = FALSE, VAR1,
                graphtitle, thickness, nodetextsize, nodesize)
}
\arguments{
  \item{c1}{
  a numeric vector of natural numbers. If the third argument \code{percomparison} is set as \code{TRUE}, \code{c1} represents treatment ID; otherwise, it represents study ID. 
}
  \item{c2}{
  a numeric vector of natural numbers indicating treatment ID.
}
  \item{percomparison}{
  logical. If \code{TRUE}, \code{c1} would be treated as the treatment in the first arm of a comparison while \code{c2} would be the second arm; otherwise, \code{c1} is the study ID column while \code{c2} is the treatment vector.
}
  \item{trtname}{
  a vector of character string indicating the treatment names. It is optional, and the default is "treat.1", "treat.2", and so on.
}
  \item{weight}{
  logical. If \code{TRUE}, the weights of edges (i.e., total number of corresponding comparison in the network study) would be shown in the network plot. The default is \code{FALSE}.
}
  \item{VAR1}{
  a numeric vector used to plot the size of each node according to a treatment characteristic, e.g., by providing a vector with the sample size randomized in each treatment. It is an optional argument, and the default is to plot the node proportional to the number of trials that include the given treatment.
}
  \item{graphtitle}{
  a character string indicating the gragh title. It is optional, and the default is to plot without a title.
}
  \item{thickness}{
  a numeric value used to change the thickness of the edges. It is optional, and the default is 10.
}
  \item{nodetextsize}{
  a numeric value used to change the text size of the node label. It is optional, and the default is 1.
}
  \item{nodesize}{
  a numeric value used to change the size of nodes. It is optional, and the default is 5.
}
}

\value{
A network plot is created by this function. Each node in the plot represents a treatment, and the links between nodes indicate studies comparing pairs of corresponding treatments.
}
\author{
Antonis Mairgiotis.
}

\section{Warning }{
When large values given to the \code{nodetextsize} argument (e.g., larger than 4), the graph might be forced out of the margins.
}

\examples{
t1 <- c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 7, 1, 1, 1, 1, 1, 1, 3, 2, 3)
t2 <- c(2, 3, 4, 4, 5, 5, 8, 8, 8, 8, 10, 7, 9, 9, 9, 10, 6, 6, 3, 4)
nma.networkplot(c1 = t1, c2 = t2, percomparison = TRUE)

## network plot for network Ara 2009
data(Ara09)
attach(Ara09)
nma.networkplot(c1 = Study.ID, c2 = Treatment, percomparison = FALSE,
                weight = TRUE, graphtitle = "Ara 2009")
detach(Ara09)

## network plot for network Lam 2007
data(Lam07)
attach(Lam07$data)
nma.networkplot(c1 = Study.ID, c2 = Treatment, percomparison = FALSE,
                trtname = Lam07$trtname, graphtitle = "Lam 2007")
detach(Lam07$data)
}