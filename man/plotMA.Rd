\name{plotMA}
\alias{plotMA}
\title{
MA plot
}
\description{
Plots a MA plot of the 2 selected samples.
}
\usage{
plotMA(obj, subset, sample1 = NULL, sample2 = NULL )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{a "corObject"}
  \item{subset}{"miRNA" or "mRNA"}
  \item{sample1}{number of the column or column name of sample 1 to use}
  \item{sample2}{number of the column or column name of sample 2 to use}
}
\value{A MA plot.}
\examples{
data(data.obj)
plotMA(data.obj,"miRNA",sample1=1,sample2=2)
plotMA(data.obj,"mRNA",sample1="H_1",sample2="D_7")
}

\keyword{plot}
\keyword{MA plot}
