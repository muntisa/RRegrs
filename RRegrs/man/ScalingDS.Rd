\name{ScalingDS}
\alias{ScalingDS}
\title{
Produce a scaled or centered data frame.
}
\description{
ScalingDS produces a transformed data frame (centering, scaling, etc.).
}
\usage{
ScalingDS(ds,s=1,c=1,fDet=FALSE,outFileName="ds4.scaled.csv")
}
\arguments{
 \item{ds}{a matrix or data frame. The first column corresponds to the dependent variable, and the remaining to independent variables.}
 \item{s}{type of scaling, default value is 1(normalization). Please see details.}
 \item{c}{indicator variable to define whether the first column of the data should be considered for scaling (c=1) or not (c=2).}
 \item{fDet}{a logical value (TRUE or FALSE). If TRUE the corrected data frame is included at the details folder.}
 \item{outFileName}{a character sring naming a file (could include path). The default the a csv file "ds4.scaled.csv".}
}
\details{
 ScalingDS offers normilization options (s=1), standardization (s=2), the user can define tailor-made scaling by setting s=3. 
}
\value{
 A data frame containing the sacled data based on the scale function by caret package.
}
\examples{
\dontrun{
ScalingDS(iris,'scaled_iris.csv')
}
}
\author{
Georgia Tsiliki, Cristian R. Munteanu
}
