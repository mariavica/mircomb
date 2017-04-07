\name{plotDensity}
\alias{plotDensity}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot miRNA or mRNA density
}
\description{
Plot miRNA or mRNA density, groping as you want}
\usage{
plotDensity(obj, subset, col.color=1, colors=c("turquoise", "violet"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{
a \kbd{corObject}.
}
  \item{subset}{
 \kbd{"miRNA"} or \kbd{"mRNA"}.
}
  \item{col.color}{
name of the column in the pheno slot used for coloring the samples.
}
  \item{colors}{
In case of a two-level factor grouping, colors to be used.
}

}
\value{
A density plot
}
\examples{
data(data.obj)
plotDensity(data.obj, "miRNA")
plotDensity(data.obj, "mRNA")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{plot}
\keyword{density}
