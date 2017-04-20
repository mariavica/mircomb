\name{addDiffexp}
\alias{addDiffexp}
\title{
Calculate differential expression}
\description{
Calculate miRNA or mRNA differential expression from a \kbd{corObject}
}
\usage{
addDiffexp(obj, dataset, classes, method.dif = "t.test", method.adj = "BH",
  var.t.test = FALSE, trend = FALSE)
}

\arguments{
  \item{obj}{
a \kbd{corObject}}
  \item{dataset}{
\kbd{"miRNA"} or \kbd{"mRNA"}
}
  \item{classes}{
column name of the \kbd{pheno.miRNA} or \kbd{pheno.mRNA} encoding the contrast to make. The column must contain \kbd{"0"} (reference) and \kbd{"1"} (case). Missing values are also allowed. More levels are accepted for \kbd{"anova"} method.}
  \item{method.dif}{method used for differential expression: \kbd{"t.test"}, \kbd{"wilcoxon"}, \kbd{"limma"}, \kbd{"rankprod"}, \kbd{"anova"} or \kbd{"only.fc"}
}
  \item{method.adj}{
multiple testing correction method used (only for \kbd{method.dif} = \kbd{"t.test"}, \kbd{"wilcoxon"}, \kbd{"limma"} or \kbd{"anova"}. One of \kbd{"holm"}, \kbd{"hochberg"}, \kbd{"hommel"}, \kbd{"bonferroni"}, \kbd{"BH"}, \kbd{"BY"}, \kbd{"fdr"}, \kbd{"none"} (see \kbd{\link{p.adjust}}).}
  \item{var.t.test}{
TRUE or FALSE (default). Only for \kbd{method.dif} = \kbd{"t.test"}. If TRUE, considers equality of variances in a T-test. 
}
  \item{trend}{
if TRUE use "limma-trend" method (only for \kbd{method.dif} = \kbd{"t.test"}). Recommended for log-normalised counts)
}
}
\value{
a  \kbd{corObject} with a \kbd{diffexp.miRNA} or \kbd{diffexp.mRNA} slot added, where the rownames of the data.frame are the names of the miRNAs or mRNAs. This slot is a data.frame with the following columns:
	\itemize{
		\item{FC}{: foldchange}
		\item{logratio}{: logratio}
		\item{meanExp}{: mean value log2-expression for the probe across all samples.}
		\item{pval}{: \emph{p} values. In RankProd the minimum of both \emph{p} values is reported.}
		\item{adj.pval}{: \emph{p} values corrected for multiple testing. In the RankProd the pfp (estimated percentage of false positives, which are, in theory, equivalent to False Discovery Rate) are added.}
	}
* If method "anova" is selected, then the proportion of SS related to factor is showed instead of FC and logratio.
** If method "only.fc" is selected, then only FC and logratio are computed (useful if there is no enough samples to perform a statistical test).


}
%\references{
%\url{http://www.bioconductor.org/packages/release/bioc/html/limma.html}

%\url{http://www.bioconductor.org/packages/release/bioc/html/RankProd.html}

%\url{http://www.bioconductor.org/packages/release/bioc/html/voom.html}
%}
\seealso{
\code{\link{corObject-class}},\code{\link{package:limma}, \code{\link{package:RankProd}}, \code{\link{voom}}, \code{\link[stats]{t.test}}, \code{\link{wilcoxon}}, \code{\link{aov}}, \code{\link[stats]{p.adjust}}
}

\examples{
data(miRNA)
data(mRNA)
data(pheno.miRNA)
data(pheno.mRNA)

minimal<-new("corObject",dat.miRNA=miRNA,dat.mRNA=mRNA,
	pheno.miRNA=pheno.miRNA,pheno.mRNA=pheno.mRNA)

minimal.diffexp<-addDiffexp(minimal, "miRNA", classes="DvH",
 method.dif="limma")
head(minimal.diffexp@diffexp.miRNA)

minimal.diffexp<-addDiffexp(minimal.diffexp, "mRNA", classes="DvH",
 method.dif="limma")
head(minimal.diffexp@diffexp.mRNA)

}
\keyword{differential}
