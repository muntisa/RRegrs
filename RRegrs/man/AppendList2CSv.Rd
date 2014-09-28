\name{AppendList2CSv}
\alias{AppendList2CSv}
\title{
Append list items to a CSV file.
}
\description{
AppendList2CSv appends the items of a list to a CSV file using write.table function.
}
\usage{
AppendList2CSv(l,csvFile)
}
\arguments{
 \item{l}{the list object to be writen into a CSV file.}
 \item{csvFile}{a string to define the CSV file name to be append the items of a list.}
}
\value{
 New info into CSV files.
}
\examples{
\dontrun{
AppendList2CSv(iris,"Results.csv")
}
}
\author{
Cristian R. Munteanu
}
