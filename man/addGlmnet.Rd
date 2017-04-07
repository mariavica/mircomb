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
%%     ~~Describe \code{obj} here~~
}
  \item{response}{
%%     ~~Describe \code{response} here~~
}
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
\code{package:survival}
}

\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (obj, response = "mRNAs", alpha = 0.5, upper.limit = 0) 
{
    x <- t(obj@dat.miRNA[obj@sig.miRNA, ])
    y <- t(obj@dat.mRNA[obj@sig.mRNA, ])
    if (response == "mRNAs") {
        res <- cv.glmnet(x, y, alpha = alpha, family = "mgaussian", 
            upper.limit = upper.limit)
        fullcoefs <- matrix(NA, ncol = ncol(y), nrow = ncol(x))
        for (i in 1:ncol(fullcoefs)) {
            fullcoefs[coef(res, s = "lambda.min")[[i]]@i[-1], 
                i] <- as.vector(coef(res, s = "lambda.min")[[i]]@x)[-1]
            print(i)
        }
        rownames(fullcoefs) <- colnames(x)
        colnames(fullcoefs) <- colnames(y)
        obj@cor <- fullcoefs
    }
    if (response == "mRNA") {
        fullcoefs <- matrix(NA, ncol = ncol(y), nrow = ncol(x))
        for (i in 1:ncol(y)) {
            res <- cv.glmnet(x, y[, i], alpha = alpha, family = "gaussian", 
                upper.limit = upper.limit)
            fullcoefs[, i] <- coef(res, s = "lambda.min")[-1]
        }
        rownames(fullcoefs) <- colnames(x)
        colnames(fullcoefs) <- colnames(y)
        obj@cor <- fullcoefs
    }
    if (response == "miRNA") {
        fullcoefs <- matrix(NA, ncol = ncol(y), nrow = ncol(x))
        for (i in 1:ncol(x)) {
            res <- cv.glmnet(y, x[, i], alpha = alpha, family = "gaussian", 
                upper.limit = upper.limit)
            fullcoefs[i, ] <- coef(res, s = "lambda.min")[-1]
        }
        rownames(fullcoefs) <- colnames(x)
        colnames(fullcoefs) <- colnames(y)
        obj@cor <- fullcoefs
    }
    return(obj)
  }
}

\keyword{correlation}
\keyword{glmnet}
\keyword{elasticnet}
