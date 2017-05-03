\name{plotHeatmap}
\alias{plotHeatmap}

\title{
Plot heatmaps}
\description{
Plot heatmaps of the top "n" miRNAs or mRNAs from selected miRNA-mRNA pairs
}
\usage{
plotHeatmap(obj, class, n = 50, col.color = 1, min.exp = NULL, main = NULL, pval.cutoff = NULL, grouping.col = NULL, grouping.row = NULL, order = "pval", cex.lab = 0.75)
}
\arguments{
  \item{obj}{
a \kbd{corObject}
}
  \item{class}{
\kbd{"miRNA"}, \kbd{"mRNA"} or \kbd{"both"}. When \kbd{"both"}, the heatmap is a mixture of miRNAs and mRNAs.
}
  \item{n}{
number of items to plot.
}
  \item{col.color}{character vector, column name or number which contains that will be used to label the samples.}
  \item{min.exp}{
(optional) minimum mean expression of the items to be plotted}
  \item{main}{
(optional) title of the plot}
  \item{pval.cutoff}{
(optional) maximum adjusted pvalue}
  \item{grouping.col}{
(optional) factor or character vector indicating predefined groups}
  \item{grouping.row}{
(optional) factor or character vector indicating predefined groups}
  \item{order}{criteria used to select the "top" miRNAs or mRNAs. \kbd{"pval"} or \kbd{"adj.pval"} (default), or \kbd{"logratio"} or \kbd{"FC"}. Only for class="miRNA" or "mRNA"}
  \item{cex.lab}{cex of the label axis. Default is 0.75}
}
\value{
A plot with the heatmap of the selected miRNAs or mRNAs.
}

\seealso{
\code{\link{plotCorrelation}}
}
\examples{
data(data.obj)
plotHeatmap(data.obj, "miRNA", n=100)
}
\keyword{plot}
\keyword{heatmap}
