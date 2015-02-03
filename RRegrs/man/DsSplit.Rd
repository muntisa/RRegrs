\name{DsSplit}
\alias{DsSplit}
\title{
Produce a list with training and test data frames for a dataset.
}
\description{
DsSplit split one data frame (a dataset) into training and test data frames using createDataPartition function by caret package.
}
\usage{
DsSplit(ds,trainFrac=3/4,fDet=FALSE,PathDataSet="",iSeed)
}
\arguments{
 \item{ds}{the object to be filtered out, a numeric vecror, matrix or data frame. The first column corresponds to the dependent variable, and the remaining to independent variables.}
 \item{trainFrac}{the fraction of the training set from the entire dataset; the default values is 3/4: 75\% of the cases are using for training and 25\% for the test.}
 \item{fDet}{a logical value (TRUE or FALSE). If TRUE the generated training and test data frames are included at the details folder.}
 \item{PathDataSet}{a character string naming the path for the working/results folder. }
 \item{iSeed}{the value of the seed function used to split dataset.}
}
\details{
 If fDet is TRUE, this functions generates two CSV files for training and test dataset: ds.Train.csv and ds.Test.csv.
}
\value{
 A list of two data frames for training and test sets based on the createDataPartition function by caret package.
}
\examples{
\dontrun{
DsSplit(iris,3/4,TRUE,"DataResults",1)
}
}
\author{
Georgia Tsiliki, Cristian R. Munteanu
}
