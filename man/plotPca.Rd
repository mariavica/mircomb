º\name{plotPca}
\alias{plotPca}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
PCA with miRNA or mRNA data
}
\description{
Plot Principal Components Analysis of miRNA or mRNA data
}
\usage{
plotPca(obj, subset, col.color = 1, colors = c("turquoise", "violet"), pos.leg="topleft",
	names = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{
a \kbd{corObject}.
}
  \item{subset}{
 \kbd{"miRNA"} or \kbd{"mRNA"}.
}
  \item{col.color}{
name of the column in the pheno slot used for coloring the samples.
}
  \item{colors}{
In case of a two-level factor grouping, colors to be used.
}
  \item{pos.leg}{
legend position, \bkd{"bottomright"}, \bkd{"topright"}, \bkd{"bottomleft"} or \bkd{"topleft"}
}
  \item{names}{FALSE or TRUE. Plot names of the samples.
}
  \item{\dots}{
further arguments.
}
}

\value{
A pca plot showing the principal components analysis (PCA) of miRNA or mRNA data.}


\seealso{
\code{\link{plot3d}}
}

\examples{
data(data.obj)
plotPca(data.obj,"miRNA")
plotPca(data.obj,"mRNA")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{plot}
\keyword{pca}
