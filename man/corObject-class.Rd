\name{corObject-class}
\Rdversion{1.1}
\docType{class}
\alias{corObject-class}

\title{Class \code{"corObject"}}
\description{
Class object for storing all the information for a miRNA-mRNA correlation analysis.
%%  ~~ A concise (1-5 lines) description of what the class is. ~~
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("corObject", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{dat.miRNA}:}{Object of class \code{"matrix"}. Contains the miRNA expression (rows for miRNAs and columns for samples).}
    \item{\code{dat.mRNA}:}{Object of class \code{"matrix"}. Contains the mRNA expression (rows for mRNAs and columns for samples).}
    \item{\code{pheno.miRNA}:}{Object of class \code{"data.frame"}. Rows for samples and columns for phenotypical information.}
    \item{\code{pheno.mRNA}:}{Object of class \code{"data.frame"}. Rows for samples and columns for phenotypical information.}
    \item{\code{cor}:}{Object of class \code{"matrix"}. Rows for miRNAs and columns for mRNAs.}
    \item{\code{pval}:}{Object of class \code{"matrix"}. Rows for miRNAs and columns for mRNAs.}
    \item{\code{net}:}{Object of class \code{"data.frame"}. Rows for unique miRNA:mRNA pairs and columns for their corresponding information (at least: miRNA name, mRNA name, coefficient of correlation and \emph{p} value).}
    \item{\code{diffexp.miRNA}:}{Object of class \code{"data.frame"}. Rows for miRNAs and columns for their corresponding infomation (usually FC, logratio, mean expression, \emph{p} value and corrected \emph{p} value).}
    \item{\code{diffexp.mRNA}:}{Object of class \code{"data.frame"}. Rows for mRNAs and columns for their corresponding infomation (usually FC, logratio, mean expression, \emph{p} value and corrected \emph{p} value).}
    \item{\code{sig.miRNA}:}{Object of class \code{"vector"}. Vector specifying the miRNAs that are used for correlation.}
    \item{\code{sig.mRNA}:}{Object of class \code{"vector"}. Vector specifying the mRNAs that are used for correlation.}
    \item{\code{GO.results}:}{Object of class \code{"list"}. It contains the results of a GO analysis.}
    \item{\code{info}:}{Object of class \code{"list"}. It contains the information of the tests and functions used.}
  }
}


\section{Methods}{
No methods defined with class "corObject" in the signature.
}
\details{
MiRNAs should be named according to miRBase 17, and mRNAs should be named according to HUGO gene symbol nomenclature.
}


\examples{
## minimal corObject:

data(miRNA)
data(mRNA)
data(pheno.miRNA)
data(pheno.mRNA)

minimal<-new("corObject",dat.miRNA=miRNA,dat.mRNA=mRNA,
	pheno.miRNA=pheno.miRNA,pheno.mRNA=pheno.mRNA)

str(minimal)

## corObject with more slots:
data(data.obj)
str(data.obj)
}
\keyword{classes}
