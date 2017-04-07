\name{selSubsetExprs}
\alias{selSubsetExprs}

\title{
Show differentially expressed miRNAs or mRNAs
}
\description{
Show differentially expressed miRNAs or mRNAs, according to your criteria}
\usage{
selSubsetExprs(obj, dataset, FC = NA, logratio = foldchange2logratio(FC), slope = NA , pval = NA, adj.pval = NA, min.meanExp = NA,
 up = FALSE, dw = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{
a \code{corObject}
}
  \item{dataset}{
 \kbd{"miRNA"} or \kbd{"mRNA"}.
}
  \item{FC}{
minimum absolute FoldChange cutoff
}
  \item{logratio}{
minimum absolute logratio cutoff
}
  \item{slope}{
for longitudinal-regression analysis, minimum absolute slope
}
  \item{pval}{
\emph{p} value cutoff
}
  \item{adj.pval}{
adjusted \emph{p} value cutoff
}
  \item{min.meanExp}{
minimum mean expression cutoff
}
  \item{up}{
TRUE or FALSE. Select only upregulated miRNAs or mRNAs.
}
  \item{dw}{
TRUE or FALSE. Select only upregulated miRNAs or mRNAs.
}
}

\value{
A data.frame with the selected differentially expressed miRNAs or mRNAs and their characteristics.
}

\examples{
data(data.obj)
selSubsetExprs(data.obj, "miRNA", adj.pval=0.05, FC=1.5)
selSubsetExprs(data.obj, "mRNA", adj.pval=0.05, up=TRUE)

}
\keyword{selection}
\keyword{subset}
\keyword{expression}
