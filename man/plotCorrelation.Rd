\name{plotCorrelation}
\alias{plotCorrelation}
\title{
Plot correlations
}
\description{
Plot the correlation of a miRNA and mRNA and/or the diagnostic plots of the linear regression.}
\usage{
plotCorrelation(obj, miRNA, mRNA, type = "cor", samples = "all",
  col.color = 1, i.legend = !is.na(col.color), pos.legend = "topright",
  sample.names = FALSE, pos.sample.names = 1, cex.main = 1.35,
  alternative = "two.sided", colors = c("turquoise", "violet"))
}
\arguments{
  \item{obj}{
corObject
}
  \item{miRNA}{
character, miRNA selected}
  \item{mRNA}{
character, mRNA selected}
  \item{type}{
"cor" or "residuals"
}
  \item{samples}{
colunm name of the pheno.miRNA where to look for the grouping factor. The column must contain character names.
}
  \item{col.color}{name of the column in the pheno slot used for coloring the samples.
}
  \item{i.legend}{
TRUE or FALSE. Plot legend}
  \item{pos.legend}{
set the legend position from \code{legend}  ‘"bottomright"’, ‘"bottom"’, ‘"bottomleft"’,
     ‘"left"’, ‘"topleft"’, ‘"top"’, ‘"topright"’, ‘"right"’ and ‘"center"’.}
  \item{sample.names}{
TRUE or FALSE. Plot text of sample names.}
  \item{pos.sample.names}{
Position of the text for sample names from \code{\link{text}}. Values of ‘1’, ‘2’, ‘3’ and
          ‘4’, respectively indicate positions below, to the left of,
          above and to the right of the specified coordinates. }

  \item{cex.main}{\kbd{cex} parameter for the main title.}
  \item{alternative}{alternative hypotheses used for computing the \emph{p} value of the correlation. One of \kbd{"two.sided"}, \kbd{"less"},\kbd{"greater"}.
}
  \item{colors}{In case of a two-level factor grouping, colors to be used.
}

}

\value{
A correlation plot.
}

\examples{
data(data.obj)
plotCorrelation(data.obj, miRNA="hsa-miR-107", mRNA="ACPP",
 type="cor", col.color="group", sample.names=TRUE)
plotCorrelation(data.obj, miRNA="hsa-miR-107", mRNA="ACPP",
 type="residuals")
}
\keyword{plot}
\keyword{correlation}
\keyword{residuals}
\keyword{regression}
