\name{addDatabase}
\alias{addDatabase}
\title{
Intersect correlations with an external database.}
\description{
For each miRNA-mRNA pair, add if this pair has been predicted as miRNA-target according to the desired external databases.
}
\usage{
addDatabase(obj, database = "microCosm_v5_18_numeric", pval.ref=1, dat.sum=1)
}
\arguments{
  \item{obj}{
	a \kbd{corObject} with a \kbd{cytofile} slot already defined.
}
  \item{database}{
	\kbd{"microCosm_v5_18_numeric"}, or a character vector including: \kbd{"microCosm_v5_18"}, \kbd{"targetScan_v6.2_18"} or a character with the name of the \kbd{data.frame} containing the database.
}
  \item{pval.ref}{
	only for \kbd{"microCosm_v5_18_numeric"}: \emph{p} value to set if no information is given.		
}
  \item{dat.sum}{only if \kbd{"microCosm_v5_18_numeric"} is not selected. For future purposes, the minimun number of concurrences across databases that determine that a miRNA-mRNA pair is bioinformatically predicted. By default: 1.
}
}
\details{If the database is a customised \kbd{data.frame}, the row names must follow the format "miRNA_name:mRNA_name" (check \kbd{head(microCosm_v5_18)} or \kbd{head(targetScan_v6.2_18)})
}
\value{
\item{corObject}{The same \kbd{corObject} with new columns added. The name and content of the column is dependent on the database selected:
\itemize{
	\item \kbd{microCosm_v5_18_numeric}: \kbd{pval.database}: \emph{p} value of the microCosm database. If there is no \emph{p} value in the database, then the value \kbd{pval.ref} is assigned
	\item A character vector including the names of the databases: the following columns are created:
	\itemize{
		\item For each database: column called \kbd{"dat.xxx"}, where {"xxx"} is the name of the database. The values are "0" (no target) or "1" (target).
		\item A column called \kbd{"dat.sum"} which sums, for each row, the values of all \kbd{"dat.xxx"} columns.
	}
	
	%\item \kbd{targetScan_v6.2_18}: \kbd{targetscan}: "0" (not target) or "1" (target), according to TargetScan database.
	%\item custom database: \kbd{custom.dat}: "0" (not target) or "1" (target), according to the custom database.
	%\item "PITA": "pita.score": score of the PITA database.
}
}
}
\references{
\url{http://www.targetscan.org/}

\url{http://www.ebi.ac.uk/enright-srv/microcosm/htdocs/targets/v5/}
%\url{http://genie.weizmann.ac.il/pubs/mir07/mir07_data.html}
}

\seealso{
\code{data(microCosm_v5_18)}, \code{data(targetScan_v6.2_18)}
}
\examples{
## load databases
data(microCosm_v5_18)
data(targetScan_v6.2_18)

## load corObject example
data(data.obj)

## numeric example
data.obj.datadded.numeric <- addDatabase(data.obj, "microCosm_v5_18_numeric")
head(data.obj.datadded.numeric@net)

## non-numeric example, multiple databases
data.obj.datadded.multiple <- addDatabase(data.obj,
 c("microCosm_v5_18", "targetScan_v6.2_18"))
head(data.obj.datadded.multiple@net)
}
\keyword{miRNA}
\keyword{database}
\keyword{add}
\keyword{net}
