\name{writeSif}
\alias{writeSif}

\title{
Write a SIF file
}
\description{Export the network (from \kbd{cytofile} slot) to a SIF file.
}
\usage{
writeSif(obj, file, pval.cutoff = 0.05, cor = NULL, alternative = "less",
        dat.sum = obj@info[["dat.sum"]], add.other = NULL, sub.miRNA = NULL,
        sub.mRNA = NULL, expand = FALSE, vertex.cex = "interact.table")
}
\arguments{
  \item{obj}{a \kbd{corObject}.}
  \item{file}{
file to write.
}
  \item{pval.cutoff}{maximum corrected \emph{p} value to take.
}
  \item{cor}{correlation cutoff.
}
  \item{alternative}{\kbd{"less"}: select correlations lower than correlation cutoff; \kbd{"greater"}: select correlations greater than correlation cutoff; or \kbd{"two.sided"}: select correlations with absolute value greater than correlation cutoff.
}
  \item{dat.sum}{number of minimum ocurrences across databases
}
  \item{add.other}{a character vector. Name of the dataframe that contains extra interactions (usually mRNA-mRNA interactions) that will be added to the network.
}
  \item{sub.miRNA}{character vector. Restrict to these miRNAs.
}
  \item{sub.mRNA}{character vector. Restrict to these mRNAs.
}
  \item{expand}{expand with another table. TRUE or FALSE.
}
  \item{vertex.cex}{table to use to expand
}

}
\value{A sif file, plus a node attributes file (*the node attributes still in preparation)
}
\references{
www.cytoscape.org}
\seealso{
\code{\link{writeCsv}}, \code{\link{writeSif}}, \code{\link{selSubsetCor}}, \code{\link{selSubsetExprs}}, \code{\link{openCytoscape}}
}
\examples{

## do not run
#data(data.obj)
#writeSif(data.obj,"results_sif")
}
\keyword{write}
\keyword{external file}
