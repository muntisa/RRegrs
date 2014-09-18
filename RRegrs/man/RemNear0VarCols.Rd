\name{RemNear0VarCols}
\alias{RemNear0VarCols}
\title{
Produce a data frame reduced by data predictors that have unique values or have very few unique values.
}
\description{
RemNear0VarCols filters out data predictors that have unique values or have very few unique values based on the nearZeroVar function by caret package.
}
\usage{
RemNear0VarCols(ds, fDet=FALSE, outFile="ds3.No0Var.csv")
}
\arguments{
 \item{ds}{the object to be filtered out, a numeric vecror, matrix or data frame. The first column corresponds to the dependent variable, and the remaining to independent variables.}
 \item{fDet}{a logical value (TRUE or FALSE). If TRUE the corrected data frame is included at the details folder.}
 \item{outFile}{a character sring naming a file (could include path). The default the a csv file "ds3.No0Var.csv".}
}
\value{
 A data frame containing the filtered data based on the nearZeroVar function by caret package.
}
\examples{
\dontrun{
RemNear0VarCols(iris,'filtered_iris.csv')
}
}
\author{
Cristian R. Munteanu
}
