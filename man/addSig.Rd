\name{addSig}
\alias{addSig}
\title{
Select significant miRNAs or mRNAs
}
\description{Select significant miRNAs or mRNAs that will be used for correlationg their expression.
}
\usage{
addSig(obj, dataset, FC=NA, logratio=foldchange2logratio(FC), slope=NA, pval=NA,
  adj.pval=NA, min.meanExp=NA, up=FALSE, dw=FALSE, manual=NULL) 
}
\arguments{
  \item{obj}{a \kbd{corObject}.}
  \item{dataset}{\kbd{"miRNA"} or \kbd{"mRNA"}}
  \item{FC}{minimum FoldChange (in absolute value).}
  \item{logratio}{minimum logratio (in absolute value).}
  \item{slope}{for longitudinal-regression analysis, minimum absolute slope}
  \item{pval}{maximum uncorrected \emph{p} value.}
  \item{adj.pval}{maximum corrected \emph{p} value.}
  \item{min.meanExp}{minimum mean expression.}
  \item{up}{Select only upregulated items (\kbd{TRUE} or \kbd{FALSE}).}
  \item{dw}{Select only downregulated items (\kbd{TRUE} or \kbd{FALSE}).}
  \item{manual}{character vector with miRNA or mRNA names.}
}
\value{
\item{corObject}{The same \kbd{corObject} with the slots \kbd{sig.miRNA} or \kbd{sig.mRNA} including the names of the miRNAs or mRNAs that will be used for correlation.}
}
\seealso{
\code{\link{addDiffexp}}
}
\examples{
data(data.obj)

## select the significant miRNAs and mRNAs
data.obj<-addSig(data.obj, "miRNA", adj.pval=0.05)
data.obj<-addSig(data.obj, "mRNA", adj.pval=0.05, FC=1.5)
}
