\name{addCorrelation}
\alias{addCorrelation}
\title{
Correlate miRNA and mRNA expression
}
\description{The function correlates miRNA and mRNA expression from a \kbd{corObject} and fills the \kbd{cor} and \kbd{pval} slots.
}
\usage{
addCorrelation(obj, method = "pearson", subset.miRNA = obj@sig.miRNA,
  subset.mRNA = obj@sig.mRNA, common = NULL, alternative = "less", voom = FALSE)
}

\arguments{
 \item{obj}{
a \kbd{corObject}
}
  \item{method}{
method used for computing correlation: \kbd{"pearson"} or \kbd{"spearman"}. }
  \item{subset.miRNA}{
Optional, character vector with the names of the miRNAs to correlate. It is recommended that miRNAs are added using \kbd{addSig} function.
}
  \item{subset.mRNA}{
	Optional, character vector with the names of the mRNAs to correlate. It is recommended that mRNAs are added using \kbd{addSig} function.
}
  \item{common}{
Optional, character vector with the names of the samples to correlate (the samples must appear in both miRNA and mRNA datasets.)
}
  \item{alternative}{
specification of the alternative hypothesis:  \kbd{"less"} (default), \kbd{"two-sided"} or \kbd{"greater"}.
  \item{voom}{
If \kbd{TRUE}, normalise NGS data according to voom procedure (check \kbd{limma} package), default is \kbd{FALSE}.
}
}


\value{
\code{corObject} with the slots \code{cor} and \code{pval} filled.
}

\note{
\code{\link{addCorrelation.R}} is the slow version of this function, but has the option to compute if there are any influential samples affecting the correlation values.
}

\seealso{
\code{\link{corObject-class}}, \code{\link{cor}}, \code{\link{addCorrelation.R}}, \code{\link{addSig}}
}
\examples{
data(data.obj)

data.obj.correlated<-addCorrelation(data.obj, method="pearson", alternative="less",
 subset.miRNA=c("hsa-let-7e","hsa-miR-122"), subset.mRNA=c("A1BG","A1CF"))

data.obj.correlated@cor
data.obj.correlated@pval
}
\keyword{correlation}
