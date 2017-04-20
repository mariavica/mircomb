\name{data.obj}
\alias{data.obj}
\docType{data}
\title{
Example of a \kbd{corObject}
}
\description{
Example of a \kbd{corObject} that was used for the main vignette.
}
\usage{data("data.obj")}
\format{
  The format is:
Formal class 'corObject' [package "miRComb"] with 13 slots
  ..@ dat.miRNA    : numeric matrix
  ..@ dat.mRNA     : numeric matrix
  ..@ pheno.miRNA  :'data.frame'
  ..@ pheno.mRNA   :'data.frame'
  ..@ cor          : numeric matrix
  ..@ pval         : numeric matrix
  ..@ net          :'data.frame'
  ..@ diffexp.miRNA:'data.frame'
  ..@ diffexp.mRNA :'data.frame'
  ..@ sig.miRNA    : character
  ..@ sig.mRNA     : character
  ..@ GO.results   :List 
  ..@ info         :List 
}

\source{
Modified from: \emph{Affò S, Dominguez M, Lozano JJ, Sancho-Bru P et al. Transcriptome analysis identifies TNF superfamily receptors as potential therapeutic targets in alcoholic hepatitis. Gut 2013 Mar;62(3):452-60. PMID: 22637703}
}
\references{
Affò S, Dominguez M, Lozano JJ, Sancho-Bru P et al. Transcriptome analysis identifies TNF superfamily receptors as potential therapeutic targets in alcoholic hepatitis. Gut 2013 Mar;62(3):452-60. PMID: 22637703
}
\examples{
data(data.obj)
str(data.obj)
}
\keyword{datasets}
