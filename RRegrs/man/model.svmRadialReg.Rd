
\name{model.svmRadialReg}
\alias{model.svmRadialReg}
\docType{model.svmRadialReg}

\title{
Specifications needed for Support Vector machines wit Radial basis function kernel regression model 
}

\description{
A list of parameters needed for SVM  radial kernel regression model
}

\usage{data("model.svmRadialReg")}

\format{A list}

\details{All parameters needed for svm regression model in R caret package and specification for parameter space search, specifically:
\item{label}{a character string specifying the model}
\item{library}{a character string 'kernlab' to specify the library used}
\item{type}{a charcater vector c('Regression') for SVMRFEreg function to specify the problem studied}
\item{parameters}{a 3x3 data frame of characters specifying the parameters used, namely sigma, C and epsilon, their class and label}
\item{grid}{a grid function to specify the grid of the parameter spece for sigma; arguments are x (...), y (..) and len (length of C)}
\item{loop}{NULL} 
\item{fit}{a function to specify fitting model}
\item{predict}{a function to specify prediction model}
\tem{prob}{prediction probabilities}
\item{predictors}{a function to generate predictions}
\item{tags}{a character vector including all possible models}
\item{levels}{a function for S4 VIRTUAL class objects from kernlab R library}
\item{sort}{a sort function used for parameters sigma, C, epsilon}
}

\examples{
data(model.svmRadialReg)

svmRadialReg
}
\author{Jose A. Seoane, Carlos Fernandez-Lozano}
