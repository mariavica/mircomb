\name{boxplotCorrelation}
\alias{boxplotCorrelation}
\title{
Bopxlot, full
}
\description{
Combines boxplot with correlation}
\usage{
boxplotCorrelation(obj, miRNA, mRNA, col.color = 1, pos.leg = "topright",
 colors = c("turquoise", "violet"), ...)
}
\arguments{
  \item{obj}{a \kbd{corObject}
}
  \item{miRNA}{character with the name of the miRNA to be plotted.
}
  \item{mRNA}{character with the name of the miRNA to be plotted.
}
  \item{col.color}{numeric. Number of the column in the \kbd{pheno.miRNA} slot that will be used for coloring.  
}
  \item{colors}{character vector indicating the colors to be used.
}
 \item{pos.leg}{legend position.}

 \item{\dots}{other parameters}
}

\value{
It returns a plot.
\itemize{
	\item Top-right: boxplot for the miRNA
	\item Bottom-left: boxplot for the mRNA
	\item Bottom-right: plot of the correlation (see \code{\link{plotCorrelation}}.)
}

}

\seealso{
\code{\link{plotCorrelation}}
}
\examples{
data(data.obj)
boxplotCorrelation(data.obj, "hsa-miR-107", "ACPP", pos.leg="topleft")
}
\keyword{plot}
\keyword{correlation}
\keyword{boxplot}
\keyword{regression}
