\name{plot3d}
\alias{plot3d}

\title{
PCA plot in 3D
}
\description{
PCA plot in 3D}
\usage{
plot3d(obj, subset, col.color = 1, angle = 45, colors = c("violet", "turquoise"), lty = 0, 
	cex.points = 1, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{
a \code{corObject}
}
  \item{subset}{
"miRNA" or "mRNA"
}
  \item{col.color}{
number or name of the column in the pheno slot that will be used
}
  \item{angle}{
angle orientation
}
  \item{colors}{
colors
}
  \item{lty}{
\kbd{lty} of lines. 0 means no line.
}
  \item{cex.points}{
cex
}
  \item{\dots}{
further arguments to be passed}
}

\value{
A 3d pca plot.
}


\note{
Error in factor(as.numeric(as.factor(obj@pheno.miRNA[, col.color])) +  : 
  invalid 'labels'; length 2 should be 1 or 12
}

\examples{
data(data.obj)
plot3d(data.obj,"mRNA")
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{plot}
\keyword{pca}
