\name{pearson}
\alias{pearson}
\title{
Pearson correlation with C++ code
}
\description{The function correlates miRNA and mRNA expression from a corObject and fills the \kbd{cor} and \kbd{pval} slots.
}
\usage{
pearson(obj, method = "pearson", subset.miRNA = obj@sig.miRNA, subset.mRNA = obj@sig.mRNA,
  common = NULL, d.influences = FALSE, alternative = "two.sided")
}

\arguments{
 \item{obj}{
a \kbd{corObject}
}
  \item{method}{
method used for computing correlation: \kbd{"pearson"} or \kbd{"spearman"}. }
  \item{subset.miRNA}{	Optional, character vector with the names of the miRNAs to correlate.

}
  \item{subset.mRNA}{
	Optional, character vector with the names of the mRNAs to correlate.
}
  \item{common}{
Optional, character vector with the names of the samples to correlate (the samples must appear in both miRNA and mRNA datasets.)
}
  \item{d.influences}{Compute if there are any influential samples: \kbd{TRUE} or \kbd{FALSE} (default).}
  \item{alternative}{
specification of the alternative hypothesis: \kbd{"two-sided"}, \kbd{"less"} or \kbd{"greater"}.
}


}

\details{
Slow version of the correlation function. Works on large datasets but is much more slow. Moreover, Kendall correlation is specially slow. Use always correlation function if it is possible.

If \kbd{TRUE}, a 3-dimension matrix is added to the \kbd{info} slot, labeled \kbd{"influenting.sample"}. First dimension: miRNA names; second dimension: mRNA names; third dimension: sample names; fill: Cook's Distance for a specific sample in a specific miRNA-mRNA linear regression (defined by the dimension label-names).}
\value{
corObject the slots \kbd{"cor"} and \kbd{"pval"} filled. Optionally, a matrix named \kbd{"influenting.sample"} is added to the \kbd{info} slot.

%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}

\note{This function can take a long time to complete when is applied to large datasets.
}

\seealso{
\code{corObject-class},\code{\link{cor}}, \code{correlation.alternative}, \code{cooks.distance}
}
\examples{
#data.obj<-correlation(data.obj, method = "pearson", alternative="less")
}
\keyword{correlation}
