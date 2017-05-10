\name{plotMA}
\alias{plotMA}
\title{
MA plot
}
\description{
Plots a MA plot of the 2 selected samples.
}
\usage{
plotMA(obj, subset, sample1 = NULL, sample2 = NULL, cex.lab = 1, cex.axis = 1 )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{a "corObject"}
  \item{subset}{"miRNA" or "mRNA"}
  \item{sample1}{number of the column or column name of sample 1 to use}
  \item{sample2}{number of the column or column name of sample 2 to use}
  \item{cex.lab}{cex of the label axis. Default is 1}
  \item{cex.axis}{cex of the axis. Default is 1}
}
\value{A MA plot of the selected two samples}

\seealso{
\code{\link{plot3d}}, \code{\link{plotPca}}, \code{boxplotSamples}
}
\examples{
data(data.obj)
plotMA(data.obj,"miRNA",sample1=1,sample2=2)
plotMA(data.obj,"mRNA",sample1="H_1",sample2="D_7")
}

\keyword{plot}
\keyword{MA plot}
