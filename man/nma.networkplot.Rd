\name{nma.networkplot}
\alias{nma.networkplot}
\title{
Plotting the Network
}
\description{
\code{nma.networkplot} plots a visual network presenting treatment comparisons.
}
\usage{
nma.networkplot(s.id, t.id, data, title = "", trtname, alphabetic = TRUE,
                weight.edge = TRUE, adjust.thick = 5,
                weight.node = TRUE, adjust.node.size = 10,
                node.col = "orange", edge.col = "black", text.cex = 1,
                adjust.figsizex = 1.1, adjust.figsizey = 1.1)
}
\arguments{
  \item{s.id}{
  a numeric or character vector indicating study ID, or the corresponding column name in the argument \code{data}.
}
  \item{t.id}{
  a numeric or character vector indicating treatment ID, or the corresponding column name in the argument \code{data}.
}
  \item{data}{
  an optional data frame containing the dataset for network meta-analysis. If \code{data} is specified, the previous arguments, \code{s.id} and \code{t.id}, should be specified as the corresponding column names in \code{data}; otherwise, the previous arguments use environment variables.
}
  \item{title}{
  a character string indicating plot title.
}
  \item{trtname}{
  a vector of character strings indicating the treatment names for the corresponding treatment IDs according their order in \code{t.id}. If not specified, \code{t.id} is used as treatment names.
}
  \item{alphabetic}{
  a logical value indicating whether to sort treatment nodes alphabetically in the network plot. The default is \code{TRUE}. If \code{FALSE}, treatment nodes are sorted by the treatment IDs specified in \code{t.id}.
}
  \item{weight.edge}{
  a logical value indicating whether to draw the edges proportionally to the number of direct treatment comparisons. The default is \code{TRUE}.
}
  \item{adjust.thick}{
  a positive integer indicating the maximum thickness of the edge when \code{weight} is \code{TRUE}. The default is 5.
}
  \item{weight.node}{
  a logical value indicating whether node size is proportional to the number of direct treatment comparisons which include that node. The default is \code{TRUE}.
}
  \item{adjust.node.size}{
  a positive number to adjust the node sizes when \code{weight.node} is \code{TRUE}. The default is 10.
}
  \item{node.col}{
  a character string indicating the color of treatment nodes. The default is \code{"orange"}.
}
  \item{edge.col}{
  a character string indicating the color of edges between treatments nodes. The default is \code{"black"}.
}
  \item{text.cex}{
  a positive integer indicating the sizes of treatment names placed around/on the corresponding nodes. The default is 1.
}
  \item{adjust.figsizex}{
  a positive number used to adjust the plot width. The default is 1.1.
}
  \item{adjust.figsizey}{
  a positive number used to adjust the plot height. The default is 1.1.
}
}

\value{
A network plot is generated. Each node represents a treatment, and the edges indicate the direct comparisons between the two treatments in various studies.
}

\references{
Lin L, Zhang J, Hodges JS, and Chu H (2017). "pcnetmeta: An R package for arm-based network meta-analysis." \emph{J Stat Softw}, \bold{80}(5), 1--25.
}

\examples{
data(smoke)
nma.networkplot(s.id, t.id, data = smoke, title = "Smoke Cessation",
  trtname = c("NC", "SH", "IC", "GC"))
# NC: No contact; SH: Self-help
# IC: individual counselling; GC: group counselling

data(diabetes)
nma.networkplot(s.id, t.id, data = diabetes, title = "Diabetes",
  trtname = c("Diuretic", "Placebo", "b-blocker", "CCB", "ACE inhibitor",
  "ARB"))
}