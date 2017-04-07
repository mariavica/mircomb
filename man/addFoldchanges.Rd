\name{addFoldchanges}
\alias{addFoldchanges}
\title{
Add foldchanges to a cytofile
}
\description{
Adds information regarding to miRNA and mRNA differential expression on the cytofile slot of a \kbd{corObject}}
\usage{
addFoldchanges(obj, add.pvals = FALSE)
}
\arguments{
  \item{obj}{
	a \kbd{corObject} with a \kbd{cytofile}, \kbd{diffexp.miRNA} and \kbd{diffexp.mRNA} slots already defined.
}
  \item{add.pvals}{
	TRUE or FALSE. TRUE if p-values should be added to the net slot.
}

}
\value{
	\item{corObject}{The same \kbd{corObject}, but with the following columns added in the {cytofile} slot:
	\itemize{
		\item logratio.miRNA
		\item logratio.mRNA
		\item meanExp.miRNA
		\item meanExp.mRNA
	}
}
}

\seealso{
\code{addDiffexp}
}
\examples{
## obtain minimal net slot
data(data.obj)
data.obj@net <- data.obj@net[,-c(5:ncol(data.obj@net))]
head(data.obj@net)

## add the foldchanges from diffexp.miRNA and diffexp.mRNA slots
data.obj<-addFoldchanges(data.obj)
head(data.obj@net)
}
\keyword{differential expression}
\keyword{add}
\keyword{net}

