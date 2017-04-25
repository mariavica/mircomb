\name{plotNetwork}
\alias{plotNetwork}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot a network
}
\description{
Plot the network of the selected miRNA-mRNa interactions, with selected features}
\usage{
plotNetwork(obj, pval.cutoff = 0.05, score.cutoff = NULL, sub.miRNA = NULL,
 sub.mRNA = NULL, names = TRUE, dat.sum = obj@info[["dat.sum"]],
 add.other = NULL, vertex.cex = NULL, n = NULL, node.size = 1.5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{
a \kbd{corObject}.
}
  \item{pval.cutoff}{
\emph{p} value cutoff.
}
  \item{score.cutoff}{
maximum score.
}
  \item{sub.miRNA}{
Optional, character vector with the names of the miRNAs to correlate.
}
  \item{sub.mRNA}{
Optional, character vector with the names of the miRNAs to correlate.
}
  \item{names}{
Plot the names of the miRNAs and mRNAs. TRUE or FALSE.
}
  \item{dat.sum}{
minimum occurrences between databases.
}
  \item{add.other}{
Optional, character vector: name of the dataframe containing additional interactions (usually mRNA-mRNA interactions) that will also be displayed.
}
  \item{vertex.cex}{
Optional, character vector: name of the dataframe containing the relative size for each node in the network. For example \kbd{"interact.table"}.
}
  \item{n}{
maximum number of interactions.
}
  \item{node.size}{
Size of the node
}
}

\details{The colurs are representative of the interactions. \kbd{"interact.table"} can be loaded from \kbd{miRData} package using \kbd{data(interact.table)}.}

\value{
A network Plot
}
\examples{
data(data.obj)
data(interact.table)
plotNetwork(data.obj, pval.cutoff=0.01, dat.sum=1,
 vertex.cex="interact.table", names=FALSE)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{plot}
\keyword{network}
