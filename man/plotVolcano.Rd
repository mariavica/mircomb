\name{plotVolcano}
\alias{plotVolcano}
\title{
Volcano plot
}
\description{
Plots a Volcano plot from the differential expression analysis.
}
\usage{
plotVolcano(obj, subset, FC1 = 1.5, FC2 = 2, FDR = 0.05, cex = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{a "corObject"}
  \item{subset}{"miRNA" or "mRNA"}
  \item{FC1}{first FoldChange cutoff}
  \item{FC2}{second FoldChange cutoff}
  \item{FDR}{significance cutoff (FDR)}
  \item{cex}{cex value}
}
\value{A volcano plot.}
\examples{
data(data.obj)
plotVolcano(data.obj,"miRNA")
plotVolcano(data.obj,"mRNA")
}

\keyword{plot}
\keyword{volcano plot}
