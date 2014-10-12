\name{SVMRFEreg}
\alias{SVMRFEreg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fitting Recursive Feature Elimination - Random Forest Models.
}
\description{
SVMRFEreg fits RFE SVM regression models and returns resampling based performance measure using the train function by caret package. 
}
\usage{
SVMRFEreg(my.datf.train,my.datf.test,sCV,iSplit1,
        fDet=F,outFile="",cs=c(1,5,15,50),eps=c(0.01,0.1,0.3))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{my.datf.train}{the training data set; an object where samples are in rows and features are in columns. The first column should be a numeric or factor vector containing the outcome for each sample.}
  \item{my.datf.test}{the test data set.}
  \item{sCV}{A string or a character vector specifying which resampling method to use. See details below.}
  \item{iSplit1}{a number indicating from which splitting of the data, the above train and test sets are derived. The default value is 1.}
  \item{fDet}{A logical value for saving model statistics; the default value is FALSE. See below for details.}
  \item{outFile}{A string specifying the output file (could include path) for details (fDet=TRUE).}
\item{cs}{A numeric vector specifying the possible values of the cost function for the SVM mimization.}
\item{eps}{A numeric vector specifying the possible values epsilon function eps-svr insensitive-loss function (more details in kernlab package)}
  }
        \details{
              RMSE is the summary metric used to select the optimal model.

              RFE SVM uses the RFE function of caret and the eps-svr function in ksvm from the kernlab package modeling to obtaing the best svm model with the best feature set.              
              To control the computational nuances of the train function, trainControl is used; number of folds or resampling iterations is set to 10, and the number of completed set of folds is set to 10 (for repeated k-fold cross-validation). 
              
              sCV can take the following values: boot, boot632, cv, repeatedcv, LOOCV, LGOCV (for repeated training/test splits), none (only fits one model to the entire training set), oob (only for random forest, bagged trees, bagged earth, bagged flexible discriminant analysis, or conditional tree forest models), "adaptive_cv", "adaptive_boot" or "adaptive_LGOCV".
              
              If fDet=TRUE, the following output is produced: a CSV file with detailed statistics about the regression model (Regression method, splitting number, cross-validation type, Training set summary, Test set summary, Fitting summary, List of predictors, Training predictors, Test predictors, resampling statistics, features importance, residuals of the fitted model, assessment of applicability domain (leverage analysis, Cook`s distances, points influence)), 5-12 plots for fitting statistics as a PDF file for each splitting and cross-validation method (Training Yobs-Ypred, Test Yobs-Ypred, Feature Importance, Fitted vs. Residuals for Fitted Model, Leverage for Fitted Model, Cook`s Distance for Fitted Model, 6 standard fitting plots using plot function with cutoff.Cook). 
              
              }
 \value{
    A list is returned containing:
      \item{stat.values}{model`s statistics}
                         \item{model}{the full svmrfe model, i.e. a list of class train}
      }
                         \examples{
                         \dontrun{
                         fDet <- FALSE   
                         iSeed <- i                 
                         
                         # the fraction of training set from the entire dataset;
                         trainFrac <- 0.75
                         
                         # dataset folder for input and output files
                         PathDataSet <- 'DataResults' 
                         
                         # upload data set
                         ds <- read.csv(ds.Housing,header=T)
                         
                         # split the data into training and test sets
                         dsList  <- DsSplit(ds,trainFrac,fDet,PathDataSet,iSeed) 
                         ds.train<- dsList$train
                         ds.test <- dsList$test
                         
                         # types of cross-validation methods
                         CVtypes <- c('repeatedcv','LOOCV')
                         
                         outLM<- 'svmRFEoutput.csv'
                         
                         svmrfe.fit <- SVMRFEreg(ds.train,ds.test,CVtypes[1],iSplit=1,fDet=F,outFile=outLM)
                         }
                         }
                         \author{
                         Jose A. Seoane, Carlos Fernandez-Lozano
                         }