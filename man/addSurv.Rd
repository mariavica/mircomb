\name{addSurv}
\alias{addSurv}
\title{
Calculate hazard ratio -- STILL IN TEST!!! -- }
\description{
Calculate hazard ratio of miRNAs or mRNAs from a corObject -- STILL IN TEST!!! --
}
\usage{
addSurv(obj, dataset, time, event, adjusting = NULL)
}

\arguments{
  \item{obj}{a corObject.
corObject}
  \item{dataset}{
\kbd{"miRNA"} or \kbd{"mRNA"}.
}
  \item{time}{colname of the time column.
}
  \item{event}{colname of the event column.
}
  \item{adjusting}{(optional) colname of the adjusting variables.
}


}
\value{
\item{corObject}{ with a diffexp.miRNA or diffexp.mRNA slot added. The diffexp is a data.frame with:
	\itemize{
		\item{coef}{: risk associated to the expression of the miRNA/mRNA}
		\item{pval}{: \emph{p} values}
		\item{adj.pval}{: \emph{p} values corrected for multiple testing (FDR)}
	}

	}

}
\references{
\url{https://cran.r-project.org/web/packages/survival/index.html}

}
\seealso{
\code{package:survival}
}
\examples{
data(miRNA)
data(mRNA)
data(pheno.miRNA)
data(pheno.mRNA)

minimal<-new("corObject",dat.miRNA=miRNA,dat.mRNA=mRNA,
	pheno.miRNA=pheno.miRNA,pheno.mRNA=pheno.mRNA)

minimal@pheno.miRNA$time<-runif(nrow(minimal@pheno.miRNA),1,20)
minimal@pheno.miRNA$event<-rbinom(nrow(minimal@pheno.miRNA),1,0.5)
minimal@pheno.miRNA$a1<-rnorm(nrow(minimal@pheno.miRNA),1,0.5)
minimal@pheno.miRNA$a2<-rnorm(nrow(minimal@pheno.miRNA),1,0.5)


#plotSurv(minimal, "miRNA", "hsa-let-7c", "time", "event")

minimal.diffexp<-addSurv(minimal, "miRNA", "time", "event",c("a1","a2"))
head(minimal.diffexp@diffexp.miRNA)
minimal.diffexp<-addSurv(minimal, "miRNA", "time", "event",c("a1"))
head(minimal.diffexp@diffexp.miRNA)
minimal.diffexp<-addSurv(minimal, "miRNA", "time", "event")
head(minimal.diffexp@diffexp.miRNA)




}
\keyword{survival}
