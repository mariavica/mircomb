\name{addDiffexp}
\alias{addDiffexp}
\title{
Calculate differential expression}
\description{
Calculate differential expression from a corObject
}
\usage{
addDiffexp(obj, dataset, classes, method.dif = "t.test", method.adj = "BH",
 var.t.test = FALSE, trend = FALSE)
}

\arguments{
  \item{obj}{
corObject}
  \item{dataset}{
\kbd{"miRNA"} or \kbd{"mRNA"}
}
  \item{classes}{
column name of the column to compare in the phenodata slot. The column must contain \kbd{"0"} (reference) and \kbd{"1"} (case). Missing values are also allowed.}
  \item{method.dif}{
\kbd{"t.test"}, \kbd{"wilcoxon"}, \kbd{"limma"}, \kbd{"rankprod"}, \kbd{"anova"} or \kbd{"only.fc"}
}
  \item{method.adj}{
Multiple testing correction method used (only for \kbd{"t.test"}, \kbd{"wilcoxon"} or \kbd{"limma"}. One of \kbd{"holm"}, \kbd{"hochberg"}, \kbd{"hommel"}, \kbd{"bonferroni"}, \kbd{"BH"}, \kbd{"BY"}, \kbd{"fdr"}, \kbd{"none"}.}
  \item{var.t.test}{
TRUE or FALSE (default). If TRUE, considers equality of variances in a T-test. 
}
  \item{trend}{
if TRUE and method.dif="limma", use "limma-trend" method. (Recommended for log-normalised counts)
}
}
\value{
\item{corObject}{ with a diffexp.miRNA or diffexp.mRNA slot added. The diffexp is a data.frame with:
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

}
\references{
\url{http://www.bioconductor.org/packages/release/bioc/html/limma.html}

\url{http://www.bioconductor.org/packages/release/bioc/html/RankProd.html}

\url{http://www.bioconductor.org/packages/release/bioc/html/voom.html}
}
\seealso{
\code{package:limma}, \code{package:RankProd}
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

minimal.diffexp<-addSig(minimal.diffexp,"miRNA",pval=1)
minimal.diffexp<-addSig(minimal.diffexp,"mRNA",pval=1)

minimal<-addNet(minimal.diffexp)
minimal<-addFoldchanges(minimal)
head(minimal@net)

}
\keyword{differential}
