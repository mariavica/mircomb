\name{plotCircos}
\alias{plotCircos}
\title{
Circos plot}
\description{
Do a standard circos plot
}
\usage{
plotCircos(obj, pval.cutoff = 0.05, dat.sum = obj@info[["dat.sum"]],
 n = NULL, sub.miRNA = NULL, sub.mRNA = NULL)
}
\arguments{
  \item{obj}{
a \kbd{corObject} whith a cytofile slot
}
  \item{pval.cutoff}{
maximum corrected \emph{p} value allowed}
  \item{dat.sum}{
numeric, minimum concordance across databases
}
  \item{n}{
numeric, limit to the first "n" pairs (sorted by corrected \emph{p} value)
}
  \item{sub.miRNA}{
(optional) character vector, include only pairs containing these miRNAs.
}
  \item{sub.mRNA}{
(optional) character vector, include only pairs containing these mRNAs.
}
}
\value{
a plot
}
\references{
\url{www.circos.ca}

\url{http://cran.r-project.org/web/packages/circlize/index.html}
}
\seealso{
\kbd{plot}
}
\examples{
data(data.obj)
plotCircos(data.obj, n="50")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{plot}
\keyword{circos}
