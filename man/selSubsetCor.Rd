\name{selSubsetCor}
\alias{selSubsetCor}

\title{
Show relevant miRNA-mRNA interactions
}
\description{
Show differentially expressed miRNAs or mRNAs, according to your criteria}
\usage{
selSubsetCor(obj, pval.cutoff = 1, dat.sum = 0, sub.miRNA = NULL,
 sub.mRNA = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{
a \code{corObject}
}
  \item{pval.cutoff}{
minimum adj.pval to consider
}
  \item{dat.sum}{
minimum number of concurrences across databases
}
  \item{sub.miRNA}{
optional character vector, limit to these miRNAs
}
  \item{sub.mRNA}{
optional character vector, limit to these mRNAs
}

}

\value{
A data.frame with the selected miRNA-mRNA pairs.
}

\examples{
data(data.obj)
selSubsetCor(data.obj, pval.cutoff=0.05, dat.sum=2)
}

\keyword{selection}
\keyword{subset}
\keyword{expression}
