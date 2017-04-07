\name{summary.corObject}
\alias{summary}
\alias{summary.corObject}
\title{
Brief report of a \kbd{corObject}
}
\description{Tests if a \kbd{corObject} is valid and prints the information about the tests performed.
}
\usage{
\method{summary}{corObject}(object, ...)
}
\arguments{
  \item{object}{a \kbd{corObject}.}
	\item{...}{other.}
}


\value{
Simple version of printing of the report.
}

\seealso{
\code{\link{mkReport}}
}
\examples{
data(data.obj)
summary(data.obj)
}
