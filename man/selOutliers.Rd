\name{selOutliers}
\alias{selOutliers}
\title{
select outliers from PCA
}
\description{
select outliers from PCA
}
\usage{
selOutliers(obj, subset, method = "aq.plot", delete = FALSE, add.pheno = TRUE, n.dim = 2)
}
\arguments{
  \item{obj}{
a corObject.
}
  \item{subset}{
"miRNA" or "mRNA".
}
  \item{method}{
"aq.plot"
}
  \item{delete}{
TRUE or FALSE. If TRUE, outlier samples are removed.
}
  \item{add.pheno}{
TRUE or FALSE. If TRUE, "is.outlier" column is added to pheno slot and then a pcaPlot highlighting outlier samples is produced.
}
  \item{n.dim}{
number of components of the PCA to use.
}
}
\details{


%%  ~~ If necessary, more details than the description above ~~
}
\value{
\itemize{
	\item If delete=FALSE and add.pheno=FALSE, character vector with the names of the outlier samples.
	\item If delete=FALSE and add.pheno=TRUE, a corOjbect with a column called is.outlier added to the pheno slot indicating if a sample is outlier or not.
	\item If delete=TRUE, a corObject without the outlier samples in the data and pheno slot.
}
}
\references{
https://cran.r-project.org/web/packages/mvoutlier/index.html
}

\seealso{
\code{mvoutlier}, \code{\link{plotPca}}
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
