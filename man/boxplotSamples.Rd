\name{boxplotSamples}
\alias{boxplotSamples}
\title{
Boxplots of samples expression
}
\description{Plot boxplots of the miRNA or mRNA expression, for each sample. It is possible to colour the samples according to a phenotypical description.
}
\usage{
boxplotSamples(obj, subset, col.color = 1, las = 1, colors = c("turquoise", "violet"))
}
\arguments{
  \item{obj}{
a \kbd{corObject}}
  \item{subset}{
 \kbd{"miRNA"} or \kbd{"mRNA"}
}
  \item{col.color}{
number or name of the phenotype column which define the grouping variables.}
  \item{las}{las parameter
}
  \item{colors}{character vector indicating the colors to be used.
}
}
\seealso{
\code{plotCordist}, \code{boxplotCorrelation}
}
\examples{
data(data.obj)
boxplotSamples(data.obj, "miRNA", col.color=1)
boxplotSamples(data.obj, "mRNA", col.color=1)
}

\keyword{boxplot}
\keyword{plot}
