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
\kbd{"miRNA"} or \kbd{"mRNA"}
}
  \item{col.color}{
number or name of the column in the pheno slot that will be used to color the samples
}
  \item{angle}{
angle orientation
}
  \item{colors}{
character vector with the colors that will be used for each group in col.color.
}
  \item{lty}{
\kbd{lty} of lines that are plotted parallel to the z-axis for each sample. 0 means no line.
}
  \item{cex.points}{
cex scaling of the dots.
}
  \item{\dots}{
further arguments to be passed}
}

\value{
A 3d pca plot.
}

\seealso{
\code{\link{plotPca}}
}



\note{
This error: 

Error in factor(as.numeric(as.factor(obj@pheno.miRNA[, col.color])) +  : 
  invalid 'labels'; length 2 should be 1 or 12
  
Is due to invalid length of \kbd{colors} option. Use character vectors with the same length of the number of groups in \kbd{col.color}.
}

\examples{
data(data.obj)
plot3d(data.obj,"mRNA")
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{plot}
\keyword{pca}
