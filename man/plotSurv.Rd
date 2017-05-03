\name{plotSurv}
\alias{plotSurv}
\title{
Plot Kaplan-Meier curve
}
\description{
Plot Kaplan-Meier curve of miRNA or mRNA
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
A plot with the Kaplan-Meier.
}
\references{
\code{survival} package: https://cran.r-project.org/web/packages/survival/index.html
}
\seealso{
\code{\link{addSurv}}
}
\examples{
#data(data.obj)
#plotSurv(data.obj, "miRNA", item = ...)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{plot}
\keyword{survival}
