\name{addScore}
\alias{addScore}

\title{
Create a score}
\description{
Create a score, needed for plot.network.corObject function}
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
A \kbd{corObject} in which a column containing the score values (\kbd{score}) has been added to the \kbd{net} slot.
}

\seealso{
\code{\link{addFoldchanges}}, \code{\link{addDiffexp}}
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

