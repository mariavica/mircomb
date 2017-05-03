\name{removeSamp}
\alias{removeSamp}
\title{
Remove samples or miRNA/mRNA}
\description{
Remove samples or miRNA/mRNA from a corObject.
}
\usage{
removeSamp(obj, dataset, samples = NA, genes = NA, keep=FALSE)
}
\arguments{
  \item{obj}{
a \kbd{corObject}
  \item{dataset}{
\kbd{"miRNA"} or \kbd{"mRNA"}}
  \item{samples}{
colnames of the samples to be removed.
}
  \item{genes}{
rownames of the genes (mRNA or miRNA) to be removed.
}
  \item{keep}{
TRUE (keep given colnames/rownames) or FALSE. By default, FALSE.
}

}
\details{
Genes are removed from miRNAdat/mRNAdat slots. Samples are removed from both pheno.miRNA/pheno.mRNA and miRNAdat/mRNAdat slots.

Take into account that minimum number of required samples is 2, otherwise the function will give an error.
}
\value{
a \kbd{corObject} with the selected samples already removed.
}

\examples{
data(data.obj)

dim(data.obj@dat.miRNA)
data.obj<-removeSamp(data.obj,"miRNA",samples="D_3",genes="hsa-miR-200c")
dim(data.obj@dat.miRNA)

colnames(data.obj@dat.mRNA)
data.obj<-removeSamp(data.obj,"mRNA",samples=c("D_1","D_2"),keep=TRUE)
colnames(data.obj@dat.mRNA)

}
\keyword{remove}
