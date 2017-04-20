\name{mkReport}
\alias{mkReport}
\title{
Creates a pdf report}
\description{
Creates a pdf report summarizing the contents of the \code{corObject}}
\usage{
mkReport(obj, file, title = "Default \\\\texttt{miRComb} output")
}
\arguments{
  \item{obj}{
a \code{corObject}
}
  \item{file}{
name of the file, for example "myExampleReport"
}
  \item{title}{
Title of the report
}
}
\details{
Documents \kbd{myExampleReport.tex} and \kbd{myExampleReport.pdf} will be created on the working directory.

See PLoS ONE publication for more details: Vila-Casadesús et al., "MiRComb: an R package to analyse miRNA-mRNA interactions. Examples across five digestive cancers".  PLoS ONE, 2016.
}

\note{
This function only works in Linux computers, with LaTeX and texlive already configured.

Some known problems and solutions:
\itemize{
\item If this happens:
\code{!  ==> Fatal error occurred, no output PDF file produced!}

Try:
\code{sudo apt-get install texlive-recommended-fonts}
or just
\code{sudo apt-get install texlive-full}

\item If you have problems with xcolor package:
\code{xcolor.sty not found}
\code{sudo apt-get install latex-xcolor}

\item If you have problems with tikz package:
\code{tikz.sty not found}
\code{sudo apt-get install pgf}

}

}

\examples{
### do not run

#data(data.obj)
#mkReport(data.obj,"myExampleReport")

### documents myExampleReport.tex and myExampleReport.pdf will be created
}
\keyword{report}
\keyword{LaTeX}
\keyword{pdf}
