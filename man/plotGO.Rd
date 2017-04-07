\name{plotGO}
\alias{plotGO}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot GO enriched terms
}
\description{
From RamiGO package, plots the hierarchy of the significant GO terms
}
\usage{
plotGO(obj, type, ontology, fdr = 0.05, filename = "GO_tree_default")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{
a \kbd{corObject}.
}
  \item{type}{
\kbd{"GO"} or \kbd{"KEGG"}.
}
  \item{ontology}{
For GO terms, one of \kbd{"BP"}, \kbd{"CC"} or \kbd{"MF"}. For KEGG terms, use \kbd{"KEGG"}.
}
  \item{fdr}{
FDR cutoff.
}
  \item{filename}{
name of the TIFF figure.
}
}

\value{
A file with the GO terms plotted.
}
\examples{

### do not run
# data(data.obj)
# plotGO(data.obj,"GO","BP", fdr = 0.05, filename= "GO_example")

### this will create the file "GO_example.png"
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{plot}
\keyword{GO}
\keyword{external file}
