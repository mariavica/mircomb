\name{addScore}
\alias{addScore}

\title{
Create a score}
\description{
Create a score for each miRNA:mRNA pair associated to their relative miRNA and mRNA FoldChange}
\usage{
addScore(obj)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{
a \kbd{corObject}, with a \kbd{net} slot containing \kbd{logratio.miRNA} and \kbd{logratio.mRNA} columns.
}
}
\value{
a \kbd{corObject} in which a column containing the score values (\kbd{score}) has been added to the \kbd{net} slot.
}

\note{
This score is needed by \code{\link{plotNetwork}} function.
}

\seealso{
\code{\link{corObject-class}}, \code{\link{addFoldchanges}}, \code{\link{addDiffexp}}, \code{\link{plotNetwork}}
}

\examples{
data(data.obj)
data.obj@net$score<-NULL
head(data.obj@net)

data.obj<-addScore(data.obj)
head(data.obj@net)
}

\keyword{score}
\keyword{differential expression}

