\name{AppendList2txt}
\alias{AppendList2txt}
\title{
Append list items to a CSV file.
}
\description{
AppendList2txt appends the items of a list to a text file using write.table function.
}
\usage{
AppendList2txt(l,csvFile)
}
\arguments{
 \item{l}{the list object to be writen into a text file.}
 \item{csvFile}{a string to define the text file name whre the items of a list will be appended.}
}
\details{
 AppendList2txt appends data using space separation, without quotes, column names or row names. 
}
\value{
 New info into TXT files.
}
\examples{
\dontrun{
AppendList2txt(iris,"Results.csv")
}
}
\author{
Cristian R. Munteanu
}
