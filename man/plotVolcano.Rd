\name{plotVolcano}
\alias{plotVolcano}
\title{
Volcano plot
}
\description{
Volcano plot of miRNA or mRNA differential expression analysis}
\usage{
plotVolcano(obj, subset, FC1 = 1.5, FC2 = 2, FDR = 0.05, cex = 1, cex.lab = 1, cex.axis = 1 )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{a \kbd{"corObject"} with a \kbd{"diffexp.miRNA"} or \kbd{"diffexp.mRNA"} slot already defined}
  \item{subset}{\kbd{"miRNA"} or \kbd{"mRNA"}}
  \item{FC1}{first FoldChange cutoff}
  \item{FC2}{second FoldChange cutoff}
  \item{FDR}{significance cutoff (FDR)}
  \item{cex}{cex value for the dots. Default is 1}
  \item{cex.lab}{cex of the label axis. Default is 1}
  \item{cex.axis}{cex of the axis. Default is 1}

}
\value{A volcano plot of the miRNA or mRNA data from differential expression analysis. Items with FDR lower than the limit are highlighted in yellow. Items with FDR lower than the limit and higher than FC1 are highlighted in orange. Items with FDR lower than the limit and higher than FC2 are highlighted in red.

}

\seealso{
\code{\link{addSurv}}
}

\examples{
data(data.obj)
plotVolcano(data.obj,"miRNA")
plotVolcano(data.obj,"mRNA")
}

\keyword{plot}
\keyword{volcano plot}
