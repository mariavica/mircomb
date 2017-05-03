\name{selSubsetExprs}
\alias{selSubsetExprs}

\title{
Select differentially expressed miRNAs or mRNAs
}
\description{
Select differentially expressed miRNAs or mRNAs, according to your criteria}
\usage{
selSubsetExprs(obj, dataset, FC = NA, logratio = foldchange2logratio(FC), slope = NA , pval = NA, adj.pval = NA, min.meanExp = NA,
 up = FALSE, dw = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{
a \kbd{corObject}
}
  \item{dataset}{
 \kbd{"miRNA"} or \kbd{"mRNA"}.
}
  \item{FC}{
minimum absolute FoldChange of selected miRNAs or mRNAs
}
  \item{logratio}{
minimum absolute logratio cutoff of selected miRNAs or mRNAs
}
  \item{slope}{
for longitudinal-regression analysis, minimum absolute slope  of selected miRNAs or mRNAs
}
  \item{pval}{
maximum \emph{p} value of selected miRNAs or mRNAs
}
  \item{adj.pval}{
maximum adjusted \emph{p} value  of selected miRNAs or mRNAs
}
  \item{min.meanExp}{
minimum mean expression of selected miRNAs or mRNAs
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

\seealso{
\code{\link{writeCsv}, \code{\link{writeExcel}}, \code{\link{writeSif}}
}

\examples{
data(data.obj)
selSubsetExprs(data.obj, "miRNA", adj.pval=0.05, FC=1.5)
selSubsetExprs(data.obj, "mRNA", adj.pval=0.05, up=TRUE)

}
\keyword{selection}
\keyword{subset}
\keyword{expression}
