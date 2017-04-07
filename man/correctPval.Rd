\name{correctPval}
\alias{correctPval}
\title{
Correct \emph{p} values}
\description{
Correct \emph{p} values of one column of a \kbd{cytofile} slot.}
\usage{
correctPval(obj, method.adj = "BH", pval= "pval")
}

\arguments{
  \item{obj}{
a \kbd{corObject} with a \kbd{cytofile} slot already defined.
}
  \item{method.adj}{
one of \kbd{"holm"}, \kbd{"hochberg"}, \kbd{"hommel"}, \kbd{"bonferroni"}, \kbd{"BH"}, \kbd{"BY"}, \kbd{"fdr"}, \kbd{"none"}.
}
  \item{pval}{name of the column with the \emph{p} values to correct.
}
}

\value{
\item{corObject}{a \kbd{corObject} in which a column containing the corrected \emph{p} values has been added to the \kbd{cytofile} slot.}
}

\seealso{
\code{\link{combinePval}}, \code{\link{p.adjust}} 
}
\examples{
data(data.obj)
data.obj <- correctPval (data.obj, method.adj="BH", pval="pval")
}

\keyword{Type II error}
