\name{addNet}
\alias{addNet}
\title{
Create a net slot
}
\description{Creates and fills the net slot from a \kbd{corObject}.
}
\usage{
addNet(obj)
}
\arguments{
  \item{obj}{a \kbd{corObject} with the slots \kbd{pval} and \kbd{cor} already defined}
}

\details{The \kbd{net} slot is a \kbd{data.frame} where each row represents a potential miRNA-mRNA interaction. The row names are the concatenation of the miRNA and mRNA name (format "miRNA:mRNA") and the columns contain all the available information for each pair. Net is sorted by miRNA name and then by mRNA name.}

\value{
 A \kbd{data.frame} with the \kbd{cytofile} slot filled.
}
\examples{
data(data.obj)
data.obj@net <- data.frame()

## create the minimal net slot
data.obj<-addNet(data.obj)
head(data.obj@net)
}

\keyword{add}
\keyword{net}
