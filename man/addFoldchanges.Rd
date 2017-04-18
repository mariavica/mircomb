\name{addFoldchanges}
\alias{addFoldchanges}
\title{
Add foldchanges to the \kbd{net} slot
}
\description{
Adds information regarding to miRNA and mRNA differential expression on the \kbd{net} slot of a \kbd{corObject}}
\usage{
addFoldchanges(obj, add.pvals = FALSE)
}
\arguments{
  \item{obj}{
	a \kbd{corObject} with a \kbd{net}, \kbd{diffexp.miRNA} and \kbd{diffexp.mRNA} slots already defined.
}
  \item{add.pvals}{
	TRUE or FALSE (default). If TRUE, p values are added to the \kbd{net} slot.
}

}
\value{
	a \kbd{corObject} with the following columns added in the {net} slot:
	\itemize{
		\item logratio.miRNA
		\item logratio.mRNA
		\item meanExp.miRNA
		\item meanExp.mRNA
	}
	plus, if add.pvals=TRUE:
	\itemize{
		\item pval.miRNA
		\item pval.mRNA
		\item adj.pval.miRNA
		\item adj.pval.mRNA
	}
}

\seealso{
\code{\link{corObject-class}}, \code{\link{addDiffexp}}, \code{\link{addDiffexp}}, 
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

