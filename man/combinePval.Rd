\name{combinePval}
\alias{combinePval}
\title{
Combine \emph{p} values
}
\description{
Combine two \emph{p} values into one.}
\usage{
combinePval(obj, pval.1 = "pval", pval.2 = "pval.database", method="stouffer",
 w=c(1,1))
}
\arguments{
  \item{obj}{\kbd{corObject} with a \kbd{cytofile} slot defined. It must have at least two columns with \emph{p} values to combine.}
  \item{pval.1}{column name (from the \kbd{cytofile} slot) of the first \emph{p} value. By default: the \emph{p} value of the correlation.}
  \item{pval.2}{column name (from the \kbd{cytofile} slot) of the second \emph{p} value. By default: the \emph{p} value from MicroCosm database.}
  \item{method}{\kbd{"stouffer"} or \kbd{"fisher"}}
  \item{w}{numeric vector of length two indicating the respective weights that will be applied to Stouffer combination. By default: no weighting.}
}

\details{
Stouffer and Weighted Stouffer (Lipták) combination is computed according to:

\deqn{p_{comb} = 1- \Phi \left(  \frac{1}{\sqrt{w_1^2+w_2^2}} \left( w_1 \left( \Phi^{-1} \left( 1- p_{1} \right) \right) + w_2 \left( \Phi^{-1} \left( 1- p_{2} \right) \right) \right) \right)}
where,
\deqn{\Phi(x) = \int_{-\infty}^{x} \frac{1}{\sqrt{2\pi}}e^{\frac{z^2}{2}} \textrm{d}z}

Fisher combination is computed according to:

\deqn{
t = -2 \left( \ln p_{1} + \ln p_{2} \right) \sim \chi_{4}^{2}
}


}


\value{
\item{corObject}{a \kbd{corObject} in which a column containing the combined \emph{p} values has been added to the \kbd{cytofile} slot.}
}
\references{
For more information about the combination methods, see: 

Zaykin D.V.. Optimally weighted Z-test is a powerful method for combining probabilities
in meta analysis. Journal of Evolutionary Biology, 2011.

Gade G., Porzelius C., Fälth M., Brase J.C., Wuttig D., Kuner R., Binder H., Sültmann
H., and Beissbarth T. Graph based fusion of miRNA and mRNA expression data
improves clinical outcome prediction of prostate cancer. BMC Bioinformatics, 12(488),
2011.
}


\examples{
data(data.obj)

## add column pval.database
data(microCosm_v5_18)
data.obj<-addDatabase(data.obj,"microCosm_v5_18_numeric")

## combine the two p-values
data.obj<-combinePval(data.obj, pval.1="pval", pval.2="pval.database")

head(data.obj@net)
}
\keyword{p.value}
\keyword{Type II error}
\keyword{meta-analysis}
