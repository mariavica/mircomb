\name{plotHclust}
\alias{plotHclust}
\title{
Plot hierarchical clustering of miRNA or mRNA samples
}
\description{
Plot hierarchical clustering of miRNA or mRNA samples
}
\usage{
plotHclust(obj, subset)
}

\arguments{
  \item{obj}{
a \kbd{corObject}.
}
  \item{subset}{
 \kbd{"miRNA"} or \kbd{"mRNA"}
}

}
\value{
A plot.
}

\seealso{
\code{\link{plot3d}}, \code{\link{plotPca}}, \code{dist}, \code{\link{plotCordist}}
}
\examples{
data(data.obj)
plotHclust(data.obj, "miRNA")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{plot}
\keyword{cluster}
\keyword{distance}
\keyword{samples}
