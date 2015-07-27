\name{glmnetModel}
\alias{glmnetModel}
\docType{glmnetModel}

\title{
Specifications needed for glmnet model 
}

\description{
A list of parameters needed for glmnet model
}

\usage{data("glmnetModel")}

\format{A list}

\details{All parameters needed for glmnet function in R caret package and specification for parameter space search, specifically:
\item{label}{a character string 'glmnet' to specify the model}
\item{library}{a character string 'glmnet' to specify the library used}
\item{type}{a charcater vector c('Regression', 'Classification') for NNreg function to specify the problem studied}
\item{parameters}{a 2x3 data frame of characters specifying the parameters used, namely alpha and lambda, their class and label}
\item{grid}{a grid function to specify the grid of the parameter spece for alpha and lambda; arguments are x (number of alpha values), y (number of lambda values) and len (length of sequence for both x,y)}
\item{loop}{a loop function over grid parameter to specify unique entries for grid; a list is produced with items loop (maximum parameter values) and submodels (remaining unique lambda values)} 
\item{fit}{a function to specify fitting model}
\item{predict}{a function to specify prediction model}
\tem{prob}{prediction probabilities}
\item{predictors}{a function to generate predictions}
\item{varImp}{ a funtion calculate importance of variables}
\item{levels}{}
\item{tags}{a character vector including all possible models}
\item{sort}{a sort function used for parameters alpha, lambda}
}

\examples{
data(glmnetModel)
}
\author{Jose A. Seoane, Carlos Fernandez-Lozano}
