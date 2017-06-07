\name{plotCircos}
\alias{plotCircos}
\title{
Circos plot}
\description{
Plot a circos plot showing the miRNA-mRNA interactions
}
\usage{
plotCircos(obj, pval.cutoff = 0.05, dat.sum = obj@info[["dat.sum"]],
 n = NULL, sub.miRNA = NULL, sub.mRNA = NULL)
}
\arguments{
  \item{obj}{
a \kbd{corObject} whith a \kbd{net} slot already definded
}
  \item{pval.cutoff}{
maximum corrected \emph{p} value of the selected miRNA-mRNA pairs}
  \item{dat.sum}{
numeric, minimum concordance across databases of the selected miRNA-mRNA pairs
}
  \item{n}{
(optional) numeric. If specified, limit to the first "n" pairs (sorted by corrected \emph{p} value)
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
\url{www.circos.ca} : Krzywinski MI, Schein JE, Birol I, Connors J, Gascoyne R, Horsman D, Jones SJ, Marra MA. Circos: An information aesthetic for comparative genomics. Genome Research, 2009.

\url{http://cran.r-project.org/web/packages/circlize/index.html}
}
\seealso{
\code{\link{genes_human_h37}}, \code{\link{mirnas_human_17_h37}}
}
\examples{
data(data.obj)
plotCircos(data.obj, n="50")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{plot}
\keyword{circos}
