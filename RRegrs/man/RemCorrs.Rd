\name{RemCorrs}
\alias{RemCorrs}
\title{
Produce a data frame reduced by data predictors that are correlated.
}
\description{
RemCorrs filters out data predictors that are correlated based on the cor function by corrplot package.
}
\usage{
RemCorrs(ds,fDet,cutoff,outFile)
}
\arguments{
 \item{ds}{the object to be filtered out, a numeric vecror, matrix or data frame. The first column corresponds to the dependent variable, and the remaining to independent variables.}
 \item{fDet}{a logical value (TRUE or FALSE). If TRUE the corrected data frame is included at the details folder.}
 \item{cutoff}{the correlation cut off (ex: 0.9).}
 \item{outFile}{a character string naming an output file (could include path). }
}
\details{
 This functions generates several output plots if fDet = TRUE: [outFile].corrMAT.csv and [outFile].corrs.png as correlation matrix and plot before removal the corralated features, [outFile].corrMAT4Selected.csv and [outFile].afterRemCorr.png as correlation matrix and plot after removing the correlated features by using cutoff value and outfile as the filtered dataset.
}
\value{
 A data frame containing the filtered data based on the cor function by corrplot package.
}
\examples{
\dontrun{
RemCorrs(iris,TRUE,0.9,'filtered_iris.csv')
}
}
\author{
Cristian R. Munteanu
}
