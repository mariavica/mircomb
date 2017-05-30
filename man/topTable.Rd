\name{topTable}
\alias{topTable}
\title{
Print or plot the top connected miRNA/mRNA}
\description{Print or plot the top connected miRNA/mRNA
}
\usage{
topTable(obj, class, pval.cutoff = 0.05, dat.sum =
    obj@info[["dat.sum"]], score.cutoff = NULL, plot = FALSE,
    names = FALSE, n = NULL, remove.names = FALSE)
}
\arguments{
  \item{obj}{
corObject
}
  \item{class}{
\kbd{"miRNA"} or \kbd{"mRNA"}.
}
  \item{pval.cutoff}{
\emph{p} value to cutoff.
}
  \item{dat.sum}{
number of minimum ocurrences across databases.
}
  \item{score.cutoff}{
maximum score allowed
}
  \item{plot}{
FALSE or TRUE.
}
  \item{names}{
FALSE or TRUE. Apart from the frequence, add the names of the hits. (In this case the object returned is a data.frame).
}
  \item{n}{maximum number of pairs to consider.
}
  \item{remove.names}{in case of a plot, omit x-axis labels.
}
  \item{table.items}{when names=TRUE, number of items to be included on the table.
}
}
\details{
If \kbd{plot=FALSE} then the table is displayed, if \kbd{plot=TRUE} the table is not displayed and a barplot is plotted. If names=FALSE a vector is returned, if \kbd{names=TRUE}, a \kbd{data.frame} is returned.
}
\value{A data.frame
}
\seealso{
\code{\link{plotHeatmap}}
}
\examples{
data(data.obj)

# get the names
topTable(data.obj, "miRNA", names=TRUE, plot=FALSE)
topTable(data.obj, "mRNA", names=TRUE, plot=FALSE)

# plot
topTable(data.obj, "miRNA", names=TRUE, plot=TRUE)

}
\keyword{plot}
\keyword{barplot}
