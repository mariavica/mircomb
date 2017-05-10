\name{boxplotCorrelation}
\alias{boxplotCorrelation}
\title{
Plot boxplot and correlation of a miRNA-mRNA pair
}
\description{
Given a miRNA and a mRNA, plots in the same a frame a boxplot of the miRNA expression, a boxplot of the mRNA expression and a scatterplot of showing the correlation between the two.}
\usage{
boxplotCorrelation(obj, miRNA, mRNA, col.color = 1, pos.leg = "topright",
 colors = c("turquoise", "violet"), ...)
}
\arguments{
  \item{obj}{a \kbd{corObject}
}
  \item{miRNA}{character with the name of the miRNA to plot
}
  \item{mRNA}{character with the name of the miRNA to plot
}
  \item{col.color}{number or name of the column in the \kbd{pheno.miRNA} or \kbd{pheno.mRNA} slot which define the grouping variables  
}
  \item{colors}{character vector indicating the colors that will be used for each variable  from \kbd{col.color}
}
 \item{pos.leg}{legend position: \kbd{"topright"}, \kbd{"bottomright"}, \kbd{"topleft"}, or \kbd{"bottomleft"} }

 \item{\dots}{other parameters}
}

\value{
a plot that includes the following figures:
\itemize{
	\item Top-right: boxplot of the miRNA
	\item Bottom-left: boxplot of the mRNA
	\item Bottom-right: scatter plot of the correlation (see \code{\link{plotCorrelation}}.)
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
