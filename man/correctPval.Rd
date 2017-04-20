\name{correctPval}
\alias{correctPval}
\title{
Correct \emph{p} values}
\description{
Correct \emph{p} values from one column of a \kbd{net} slot.}
\usage{
correctPval(obj, method.adj = "BH", pval= "pval")
}

\arguments{
  \item{obj}{
a \kbd{corObject} with a \kbd{net} slot already defined.
}
  \item{method.adj}{
one of \kbd{"holm"}, \kbd{"hochberg"}, \kbd{"hommel"}, \kbd{"bonferroni"}, \kbd{"BH"}, \kbd{"BY"}, \kbd{"fdr"}, \kbd{"none"}.
}
  \item{pval}{name of the column with the \emph{p} values to correct.
}
}

\value{
a \kbd{corObject} in which a column named \kbd{adj.pval}, which contains the corrected \emph{p} values, has been added to the \kbd{net} slot.
}

\seealso{
\code{\link{combinePval}}, \code{\link{p.adjust}} , \code{\link{addNet}} 
}
\examples{
data(data.obj)
data.obj <- correctPval (data.obj, method.adj="BH", pval="pval")
head(data.obj@net)
}

\keyword{Type II error}
