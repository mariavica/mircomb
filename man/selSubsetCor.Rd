\name{selSubsetCor}
\alias{selSubsetCor}

\title{
Select relevant miRNA-mRNA interactions
}
\description{
Select differentially expressed miRNAs or mRNAs, according to your criteria}
\usage{
selSubsetCor(obj, pval.cutoff = 1, dat.sum = 0, sub.miRNA = NULL,
 sub.mRNA = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{
a \kbd{corObject}
}
  \item{pval.cutoff}{
maximum adj.pval of the selected miRNA-mRNA pairs
}
  \item{dat.sum}{
minimum number of concurrences across databases of the selected miRNA-mRNA pairs
}
  \item{sub.miRNA}{
(optional) character vector, limit to these miRNAs
}
  \item{sub.mRNA}{
(optional) character vector, limit to these mRNAs
}

}

\value{
A data.frame with the selected miRNA-mRNA pairs.
}

\seealso{
\code{\link{writeCsv}}, \code{\link{writeExcel}}, \code{\link{writeSif}}
}

\examples{
data(data.obj)
selSubsetCor(data.obj, pval.cutoff=0.05, dat.sum=2)
}

\keyword{selection}
\keyword{subset}
\keyword{expression}
