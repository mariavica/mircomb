\name{addCorrelation.R}
\alias{addCorrelation.R}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Correlation, old version
}
\description{
Correlation, old version
}
\usage{
addCorrelation.R(obj, method = "pearson", subset.miRNA = obj@sig.miRNA,
 subset.mRNA = obj@sig.mRNA, common = NULL, d.influences = FALSE,
 alternative = "two.sided")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{
a \kbd{corObject}.
}
  \item{method}{
method used for computing correlation: \kbd{"pearson"} or \kbd{"spearman"}. 
}
  \item{subset.miRNA}{
Optional, character vector with the names of the miRNAs to correlate.  It is recommended that miRNAs are added using \kbd{addSig} function.
}
  \item{subset.mRNA}{
Optional, character vector with the names of the mRNAs to correlate.  It is recommended that mRNAs are added using \kbd{addSig} function.
}
  \item{common}{
Optional, character vector with the names of the samples to correlate (the samples must appear in both miRNA and mRNA datasets.)
}
  \item{d.influences}{
compute a matrix with the Cook's Distance of each sample in each miRNA-mRNA correlatio.
}
  \item{alternative}{
specification of the alternative hypothesis: \kbd{"less"} (default), \kbd{"two-sided"} or \kbd{"greater"}.
}
}


\details{
A more complete version of the \kbd{\link{addCorrelation}} function, but significantly slower (specially \kbd{kendall} correlation). Use always \code{\link{addCorrelation}} function whenever it is possible.

If \kbd{d.influences = TRUE}, a 3-dimension matrix is added to the \kbd{info} slot, labeled \kbd{"influenting.sample"}. First dimension: miRNA names; second dimension: mRNA names; third dimension: sample names; fill: Cook's Distance for a specific sample in a specific miRNA-mRNA linear regression (defined by the dimension label-names).}


\value{
A \kbd{corObject} with the slots \kbd{"cor"} and \kbd{"pval"} filled. Optionally, a matrix named \kbd{"influenting.sample"} is added to the \kbd{info} slot.

}

\note{This function can take a long time to complete when is applied to large datasets.
}


\seealso{
\code{corObject-class},\code{\link{cor}}, \code{\link{addCorrelation}}, \code{\link{addSig}}, \code{\link{cooks.distance}}
}

\examples{
data(data.obj)

data.obj<-addSig(data.obj,"miRNA",manual=c("hsa-let-7c","hsa-miR-122"))
data.obj<-addSig(data.obj,"mRNA",manual=c("A1BG","A1CF"))


data.obj.correlated<-addCorrelation.R(data.obj, method="pearson", alternative="less",
 d.influences=TRUE)

data.obj.correlated@cor
data.obj.correlated@pval
data.obj.correlated@info[["influenting.sample"]]
}

\keyword{correlation}
\keyword{pearson}
\keyword{Cook distance}
