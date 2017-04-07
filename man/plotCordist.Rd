\name{plotCordist}
\alias{plotCordist}
\title{
Plot distances/correlation between miRNA or mRNA samples
}
\description{
Plot distances/correlation between miRNA or mRNA samples
}
\usage{
plotCordist(obj, subset, type = "cor", method.cor = "pearson",
 method.dist = "euclidean", hierarchical = FALSE, ...)
}

\arguments{
  \item{obj}{
a \kbd{corObject}.
}
  \item{subset}{
 \kbd{"miRNA"} or \kbd{"mRNA"}
}
  \item{type}{
 \kbd{"cor"} (correlation) or \kbd{"dist"} (distance).
}
  \item{method.cor}{
method used for computing correlation: \kbd{"pearson"} or \kbd{"spearman"}. 
}
  \item{method.dist}{
method used for computing distance from \kbd{dist} function:  This must be one of
          ‘"euclidean"’, ‘"maximum"’, ‘"manhattan"’, ‘"canberra"’,
          ‘"binary"’ or ‘"minkowski"’.  Any unambiguous substring can
          be given. 
}
  \item{hierarchical}{
	TRUE or FALSE. Sort the rows using hierarchical clustering.
}

  \item{\dots}{
	further arguments.
}
}
\value{
A plot.
}
\references{
\url{http://www.phaget4.org/R/image_matrix.html}
}
\seealso{
\code{\link{plot3d}}, \code{\link{plotPca}}, \code{dist}
}
\examples{
data(data.obj)
plotCordist(data.obj, "miRNA", type="cor")
plotCordist(data.obj, "mRNA", type="dist")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{plot}
\keyword{samples}
\keyword{distance}
