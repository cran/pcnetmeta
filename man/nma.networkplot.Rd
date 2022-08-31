\name{nma.networkplot}
\alias{nma.networkplot}
\title{
Plotting Treatment Network
}
\description{
Visualizes the network of treatment comparisons.
}
\usage{
nma.networkplot(s.id, t.id, n, data, title = "", title.cex = 1, trtname,
                alphabetic = TRUE, multi.show = FALSE, multi.col,
                weight.edge = TRUE, adjust.thick = 5, weight.node = TRUE,
                weight.node.ss = FALSE, adjust.node.size = 10,
                node.col = "gray", edge.col = "black", text.cex = 1,
                adjust.figsizex = 1.1, adjust.figsizey = 1.1)
}
\arguments{
  \item{s.id}{
  a numeric or character vector specifying study ID, or the corresponding column name in the argument \code{data}.
}
  \item{t.id}{
  a numeric or character vector specifying treatment ID, or the corresponding column name in the argument \code{data}.
}
  \item{n}{
  a numeric or character vector specifying sample size, or the corresponding column name in the argument \code{data}. It is used only when the node size is proportional to the total sample size of the corresponding treatment (\code{weight.node.ss} = \code{TRUE}).
}
  \item{data}{
  an optional data frame containing the dataset of the network meta-analysis. If \code{data} is specified, the previous arguments, \code{s.id}, \code{t.id}, and \code{n}, should be specified as the corresponding column names in \code{data}; otherwise, the previous arguments use environment variables.
}
  \item{title}{
  a character string indicating the plot title.
}
  \item{title.cex}{
  the size of the plot title.
}
  \item{trtname}{
  a vector of character strings specifying the treatment names for the corresponding treatment IDs according to their order in \code{t.id}. If not specified, \code{t.id} is used as treatment names.
}
  \item{alphabetic}{
  a logical value indicating whether to sort treatment nodes alphabetically in the network plot. The default is \code{TRUE}. If \code{FALSE}, treatment nodes are sorted by the treatment IDs specified in \code{t.id}.
}
  \item{multi.show}{
  a logical value indicating whether to show multi-arm studies (if any). Multi-arm studies with different treatment sets are visualized by shaded treatment loops in different colors.
}
  \item{multi.col}{
  an optional vector of colors for visualizing multi-arm studies.
}
  \item{weight.edge}{
  a logical value indicating whether to draw the edges proportionally to the number of direct treatment comparisons. The default is \code{TRUE}.
}
  \item{adjust.thick}{
  a positive integer specifying the maximum thickness of the edge when \code{weight} is \code{TRUE}. The default is 5.
}
  \item{weight.node}{
  a logical value indicating whether the node size is weighted by the total sample size or the total number of direct treatment comparisons of the corresponding treatment (\code{TRUE}, the default) or all treatment nodes have equal sizes (\code{FALSE}).
}
  \item{weight.node.ss}{
  a logical value indicating whether the node size is proportional to the total sample size (\code{TRUE}) or the total number of direct treatment comparisons (\code{FALSE}, the default) of the corresponding treatment.
}
  \item{adjust.node.size}{
  a positive number to adjust the node sizes when \code{weight.node} is \code{TRUE}. The default is 10.
}
  \item{node.col}{
  a character string specifying the color of treatment nodes. The default is \code{"gray"}.
}
  \item{edge.col}{
  a character string specifying the color of edges between treatments nodes. The default is \code{"black"}.
}
  \item{text.cex}{
  a positive integer specifying the sizes of treatment names placed around/on the corresponding nodes. The default is 1.
}
  \item{adjust.figsizex}{
  a positive number used to adjust the plot width. The default is 1.1.
}
  \item{adjust.figsizey}{
  a positive number used to adjust the plot height. The default is 1.1.
}
}

\value{
A network plot is generated. Each node represents a treatment, and the edges indicate the direct comparisons between the two treatments.
}

\references{
Lin L, Zhang J, Hodges JS, Chu H (2017). "Performing arm-based network meta-analysis in R with the pcnetmeta package." \emph{Journal of Statistical Software}, \bold{80}(5), 1--25. <\doi{10.18637/jss.v080.i05}>
}

\examples{
data("smoke")
# NC: No contact; SH: Self-help
# IC: individual counseling; GC: group counseling
nma.networkplot(s.id, t.id, data = smoke, title = "Smoke Cessation",
  trtname = c("NC", "SH", "IC", "GC"))
nma.networkplot(s.id, t.id, n, data = smoke, title = "Smoke Cessation",
  trtname = c("NC", "SH", "IC", "GC"),
  weight.node.ss = TRUE, multi.show = TRUE)

data("diabetes")
nma.networkplot(s.id, t.id, data = diabetes, title = "Diabetes",
  trtname = c("Diuretic", "Placebo", "b-blocker", "CCB", "ACE inhibitor",
  "ARB"))
nma.networkplot(s.id, t.id, data = diabetes, title = "Diabetes",
  trtname = c("Diuretic", "Placebo", "b-blocker", "CCB", "ACE inhibitor",
  "ARB"), multi.show = TRUE)
nma.networkplot(s.id, t.id, n, data = diabetes, title = "",
  weight.node.ss = TRUE, multi.show = TRUE)
mtext(text = "Diabetes", side = 3, line = 0, font = 1, cex = 2)
}
\keyword{plot}