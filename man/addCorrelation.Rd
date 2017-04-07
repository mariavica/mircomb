\name{addCorrelation}
\alias{addCorrelation}
\title{
Correlate miRNA and mRNA expression
}
\description{The function correlates miRNA and mRNA expression from a corObject and fills the \kbd{cor} and \kbd{pval} slots.
}
\usage{
addCorrelation(obj, method = "pearson", subset.miRNA = obj@sig.miRNA,
  subset.mRNA = obj@sig.mRNA, common = NULL, alternative = "two.sided")
}

\arguments{
 \item{obj}{
a \kbd{corObject}
}
  \item{method}{
method used for computing correlation: \kbd{"pearson"} or \kbd{"spearman"}. }
  \item{subset.miRNA}{
Optional, character vector with the names of the miRNAs to correlate.
}
  \item{subset.mRNA}{
	Optional, character vector with the names of the mRNAs to correlate.
}
  \item{common}{
Optional, character vector with the names of the samples to correlate (the samples must appear in both miRNA and mRNA datasets.)
}
  \item{alternative}{
specification of the alternative hypothesis: \kbd{"two-sided"}, \kbd{"less"} or \kbd{"greater"}.
}


}


\value{
\code{corObject} the slots \code{cor} and \code{pval} filled.
}

\note{
\code{\link{addCorrelation.R}} is the slow version of this funcion, but has the option to compute if there are any influential samples.
}

\seealso{
\code{corObject-class},\code{\link{cor}}, \code{\link{addCorrelation.R}}
}
\examples{
data(data.obj)

data.obj.correlated<-addCorrelation(data.obj, method="pearson", alternative="less",
 subset.miRNA=c("hsa-let-7e","hsa-miR-122"), subset.mRNA=c("A1BG","A1CF"))

data.obj.correlated@cor
data.obj.correlated@pval
}
\keyword{correlation}
