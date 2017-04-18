\name{addGlmnet}
\alias{addGlmnet}

\title{
Add elastic net estimates -- STILL IN TEST!!! --
}
\description{
Add elastic net estimates to a corObject -- STILL IN TEST!!! --
}
\usage{
addGlmnet(obj, response = "mRNAs", alpha = 0.5, upper.limit = 0, cluster = "manual", plot = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{
a  \kbd{corObject}}
  \item{response}{
\kbd{miRNA}, \kbd{mRNA} or \kbd{mRNAs}}
  \item{alpha}{
%%     ~~Describe \code{alpha} here~~
}
  \item{upper.limit}{
%%     ~~Describe \code{upper.limit} here~~
}
  \item{cluster}{
%%     ~~Describe \code{cluster} here~~
}
  \item{plot}{
%%     ~~Describe \code{cluster} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
https://cran.r-project.org/web/packages/glmnet/index.html
}

\seealso{
\code{\link{corObject-class}}, \code{\link{package:glmnet}}
}

\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.
}

\keyword{correlation}
\keyword{glmnet}
\keyword{elasticnet}
