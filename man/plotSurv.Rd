\name{plotSurv}
\alias{plotSurv}
\title{
Plot kaplan-meier curve
}
\description{
Plot kaplan-meier curve
}
\usage{
plotSurv(obj, subset, item, time, event)
}

\arguments{
  \item{obj}{
a \kbd{corObject}.
}
  \item{subset}{
 \kbd{"miRNA"} or \kbd{"mRNA"}
}
  \item{item}{
 name of the miRNA or mRNA
}

  \item{time}{
name or number of the column in the pheno slot that contains the time.
}
  \item{event}{
name or number of the column in the pheno slot that contains the censoring (0=censored, 1=event).
}


}
\value{
A plot.
}
\references{
\code{survival} package: https://cran.r-project.org/web/packages/survival/index.html
}
\seealso{
\code{\link{plot3d}}, \code{\link{plotPca}}, \code{dist}, \code{\link{plotCordist}}
}
\examples{
#data(data.obj)
#plotSurv(data.obj, "miRNA", item = ...)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{plot}
\keyword{survival}
