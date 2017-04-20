\name{addSurv}
\alias{addSurv}
\title{
Calculate hazard ratio}
\description{
Calculate hazard ratio of miRNAs or mRNAs from a \kbd{corObject} using a Cox proportional hazards regression.
}
\usage{
addSurv(obj, dataset, time, event, adjusting = NULL)
}

\arguments{
  \item{obj}{a \kbd{corObject}}
  \item{dataset}{
\kbd{"miRNA"} or \kbd{"mRNA"}.
}
  \item{time}{colname of the time column in the \kbd{pheno.miRNA} or \kbd{pheno.mRNA} slot
}
  \item{event}{colname of the event column in the \kbd{pheno.miRNA} or \kbd{pheno.mRNA} slot
}
  \item{adjusting}{(optional) colname of the adjusting variables that will be used
}


}
\value{
a \kbd{corObject} with a \kbd{diffexp.miRNA} or \kbd{diffexp.mRNA} slot added. The slot is a data.frame in which row names are the names of the miRNAs or mRNAs and has the following columns:
	\itemize{
		\item{coef}{: risk associated to the expression of the miRNA/mRNA}
		\item{pval}{: \emph{p} values}
		\item{adj.pval}{: \emph{p} values corrected for multiple testing (FDR)}
	}



}
\seealso{
\code{\link{package:survival}}, \code{\link{coxph}}
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
