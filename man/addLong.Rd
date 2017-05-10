\name{addLong}
\alias{addLong}
\title{
Time analysis}
\description{
Time analysis of a \kbd{corObject}
}
\usage{
addLong(obj, dataset, classes, method.dif = "time.point", method.adj = "BH",
 var.t.test = FALSE, trend = FALSE)
}

\arguments{
  \item{obj}{
corObject}
  \item{dataset}{
\kbd{"miRNA"} or \kbd{"mRNA"}
}
  \item{classes}{
column name of the \kbd{pheno.miRNA} or \kbd{pheno.mRNA} encoding codification of time variable. The column must contain \kbd{"0"} (t_init) and \kbd{"1"} (t_fin), or the specific times. Missing values are also allowed.}
  \item{method.dif}{
\kbd{"time.point"} or \kbd{"linear.regression"}
}
  \item{method.adj}{
Multiple testing correction method used. One of \kbd{"holm"}, \kbd{"hochberg"}, \kbd{"hommel"}, \kbd{"bonferroni"}, \kbd{"BH"}, \kbd{"BY"}, \kbd{"fdr"}, \kbd{"none"} (see \link{p.adjust}).}
  \item{var.t.test}{
\kbd{TRUE} or \kbd{FALSE} (default). If \kbd{TRUE}, considers equality of variances in a T-test. 
}
}

\details{
If \kbd{method.dif="time.point"}, both times are compared using a t-test, time variable should be encoded as kbd{"0"} (t_init) and \kbd{"1"} (t_fin).
}

\value{
a \kbd{corObject} with a \kbd{diffexp.miRNA} or \kbd{diffexp.mRNA} slot added, where the rownames of the data.frame are the names of the miRNAs or mRNAs. Depending on the selected method, the diffexp is a data.frame with the following columns:

If \kbd{method.dif="time.point"}:
	\itemize{
		\item{FC}{: foldchange between the two selected times.}
		\item{logratio}{: logratio between the two selected times.}
		\item{meanExp}{: mean value log2-expression for the probe across all samples.}
		\item{pval}{: \emph{p} values. In RankProd the minimum of both \emph{p} values is reported.}
		\item{adj.pval}{: \emph{p} values corrected for multiple testing. In the RankProd the pfp values (estimated percentage of false positives, which are, in theory, equivalent to False Discovery Rate) are added.}
	}

If \kbd{method.dif="linear.regression"}:
	\itemize{
		\item{slope}{: slope of the linear regression.}
		\item{meanExp}{: mean value log2-expression for the probe across all samples.}
		\item{pval}{: \emph{p} values associated to the slope of the linear regression.}
		\item{adj.pval}{: \emph{p} values corrected for multiple testing.}
	}



}
%\references{
%\url{http://www.bioconductor.org/packages/release/bioc/html/limma.html}
%\url{http://www.bioconductor.org/packages/release/bioc/html/RankProd.html}
%\url{http://www.bioconductor.org/packages/release/bioc/html/voom.html}
%}
\seealso{
 \code{\link{addDiffexp}}
}
\examples{
################################
### data special preparation ###
################################

data(miRNA)
data(mRNA)

### simulated phenotypical: 3 samples in 4 time-points, one treatment.

pheno.miRNA<-data.frame(sample=rep(1:3, each=4),
	time=rep(c(1,2,5,9),3),
	time_alt=rep(c(0,NA,NA,1),3))
rownames(pheno.miRNA)<-paste(pheno.miRNA$sample,pheno.miRNA$time,sep="_")

pheno.mRNA<-pheno.miRNA	#same sample distribution

## modify sample names from miRNA and mRNA dataset
colnames(miRNA)<-rownames(pheno.miRNA)
colnames(mRNA)<-rownames(pheno.mRNA)

############################
### start of the example ###
############################

minimal<-new("corObject",dat.miRNA=miRNA,dat.mRNA=mRNA,
	pheno.miRNA=pheno.miRNA,pheno.mRNA=pheno.mRNA)

#### comparing times 1 (t_init) and 9 (t_fin):
minimal.diffexp<-addLong(minimal, "miRNA", classes="time_alt",
 method.dif="time.point")
head(minimal.diffexp@diffexp.miRNA)

#### treated as linear regression:
minimal.diffexp<-addLong(minimal, "miRNA", classes="time",
 method.dif="linear.regression")
head(minimal.diffexp@diffexp.miRNA)


}
\keyword{differential}
