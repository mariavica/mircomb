\name{selOutliers}
\alias{selOutliers}
\title{
select outliers based on PCA analysis
}
\description{
select outliers from miRNA or mRNA samples based on PCA analysis

}
\usage{
selOutliers(obj, subset, method = "aq.plot", delete = FALSE, add.pheno = TRUE, n.dim = 2)
}
\arguments{
  \item{obj}{
a \kbd{corObject}
}
  \item{subset}{
\kbd{"miRNA"} or \kbd{"mRNA"}}
  \item{method}{
method used to select the outliers \kbd{"aq.plot"}
}
  \item{delete}{
TRUE or FALSE. If TRUE, outlier samples are removed.
}
  \item{add.pheno}{
TRUE or FALSE. If TRUE, \kbd{"is.outlier"} column is added to pheno slot and then a PCA plot highlighting outlier samples is produced.
}
  \item{n.dim}{
number of components of the PCA to use.
}
}
\details{
This is an implementation of mvoutlier package. Check the origninal source for more information.
}
\value{
\itemize{
	\item If delete=FALSE and add.pheno=FALSE, character vector with the names of the outlier samples.
	\item If delete=FALSE and add.pheno=TRUE, a \kbd{corObject} with a column called is.outlier added to the pheno slot indicating if a sample is outlier or not.
	\item If delete=TRUE, a \kbd{corObject} without the outlier samples in the data and pheno slot.
}
}

\seealso{
\code{\link{mvoutlier}, \code{\link{plotPca}}
}
\examples{
data(data.obj)
#selOutliers(data.obj, "miRNA", add.pheno=FALSE)
#data.obj.out <- selOutliers(data.obj, "miRNA", add.pheno=TRUE)
#head(data.obj@pheno.miRNA)
#head(data.obj.out@pheno.miRNA)
}
\keyword{outlier}
\keyword{plot}
\keyword{pca}
