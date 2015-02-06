\name{RRegrs}
\alias{RRegrs}
\title{
Fitting a set of ten fully utilized Predictive Models.
}
\description{
The current tool is a collection of regression tools from R that could be used to search the best regression models for any dataset.
}
\usage{
RRegrs(DataFileName="ds.House.csv",
PathDataSet="DataResults", noCores=1, ResAvgs="RRegsResAvgs.csv", ResBySplits="RRegrsResAllSplits.csv", ResBest="RRegrsResBest.csv", fDet="T", fFilters="F", fScaling="T", fRemNear0Var="T", fRemCorr="T",fLM="T", fGLM="T", fPLS="T", fLASSO="T", fSVRM="T", fNN="T", fRF="T", fRFRFE="T", fSVMRFE="T", fENET="T", RFE_SVM_C="1;5;10;15;20", RFE_SVM_epsilon="0.01;0.1;0.3", cutoff=0.9, iScaling=1, iScalCol=1, trainFrac=0.75, iSplitTimes=10, noYrand=100, CVtypes="repeatedcv;LOOCV", NoNAValFile="ds.NoNA.csv",  No0NearVarFile="ds.No0Var.csv", ScaledFile="ds.scaled.csv", NoCorrFile="ds.scaled.NoCorrs.csv",
lmFile="LM.details.csv", glmFile="GLM.details.csv", plsFile="PLS.details.csv",lassoFile="Lasso.details.csv", svrmFile="SVMRadial.details.csv", nnFile="NN.details.csv", rfFile="RF.details.csv", rfrfeFile="RFRFE.details.csv", svmrfeFile="SVMRFE.details.csv",enetFile="ENET.details.csv")
}
\arguments{
 \item{DataFileName}{ a character string naming the data set file; the data set should be an object where samples are in rows and features are in columns. The first column should be a numeric or factor vector containing the outcome for each sample. Default value is the "ds.House.csv" data.}
 \item{PathDataSet}{ a character string naming the directory where the data set file is stored. Default value is the "DataResults" directory.}
 \item{noCores}{The number of cores used for calculations. Default value is 1.}
 \item{ResAvgs}{a character string naming the file where average  values of all models are stored. Default value is the "RRegsResAvgs.csv" file.}
 \item{ResBySplits}{a character string naming the file where statistics for all splits and all models are stored. Default value is the "RRegrsResAllSplits.csv" file.}
 \item{ResBest}{a character string naming the file where statistics for the best model are stored. Default value is the "RRegrsResBest.csv" file.}
 \item{fDet}{logical. If TRUE all details for all splits and models are stored.}
 \item{fFilters}{If run custom flters. Not implemented yet.}
 \item{fScaling}{logical. If scaling will be performed. Default value is TRUE.}
 \item{fRemNear0Var}{logical. If run Removal of near zero variance columns. Default value is TRUE.}
 \item{fRemCorr}{logical. If correlated features will be removed based on 'cutoff' value. Default value is TRUE.}
 \item{fLM}{logical. If Linear model will be applied to the data. Default value is TRUE.}
 \item{fGLM}{logical. If Generalized Linear model will be applied to the data. Default value is TRUE.}
 \item{fPLS}{logical. If Partial Leasts Square regression model will be applied to the data. Default value is TRUE.}
 \item{fLASSO}{logical. If Lasso regression model will be applied to the data. Default value is TRUE.}
 \item{fSVRM}{logical. If Support vector machine regression model (radial basis function kernel) will be applied to the data. Default value is TRUE.}
 \item{fNN}{logical. If Neural network model will be applied to the data. Default value is TRUE.}
 \item{fRF}{logical. If Random Forest model will be applied to the data. Default value is TRUE.}
 \item{fRFRFE}{logical. If random Forest with Random Feature Elimination model will be applied to the data. Default value is TRUE.}
 \item{fSVMRFE}{logical. If Support vector machines with Random Feature elimination model (radial basis function kernel) will be applied to the data. Default value is TRUE.}
 \item{fENET}{logical. If Elastic net model will be applied to the data. Default value is TRUE.}
 \item{RFE_SVM_C}{Support vector machines cost parameter. Default values are c(1,5,10,15,20).}
 \item{RFE_SVM_epsilon}{Support vector machines epsilon parameter. Default values are c(0.01,0.1,0.3).}
 \item{cutoff}{=0.9}{Cutoff for correlated features (default = 0.9).}
 \item{iScaling}{Type of scaling: 1 = normalization; 2 = standardization; 3 = other; any other: no scaling}
 \item{iScalCol}{Scaling columns: 1 = including dependent variable; 2: only all the features}
 \item{trainFrac}{ Fraction of training set from the entire dataset (default = 0.75); the rest of dataset is the test set.}
 \item{iSplitTimes}{Number of splittings the dataset into train and test (default = 10)}
 \item{noYrand}{Number of Y-Randomization (default = 100)}
 \item{CVtypes}{a character string indicating which cross validation method will be used. Default value is c("repeatedcv","LOOCV")}
 \item{NoNAValFile}{a character string naming the file where the data set without NA values will be stored. Default value is the "ds.NoNA.csv" file.}
 \item{No0NearVarFile}{a character string naming the file where the data set without without near zero variance values will be stored. Default value is the "ds.No0Var.csv" file.}
 \item{ScaledFile}{a character string naming the file where the scaled data set will be stored. Default value is the "ds.scaled.csv" file.}
 \item{NoCorrFile}{a character string naming the file where the data set without correlated independent variables will be stored. Default value is the "ds.scaled.NoCorrs.csv" file.}
 \item{lmFile}{a character string naming the file where LM model details will be stored. Default value is the "LM.details.csv" file.}
 \item{glmFile}{a character string naming the file where GLM model details will be stored. Default value is the "GLM.details.csv" file.}
 \item{plsFile}{a character string naming the file where PLS model details will be stored. Default value is the "PLS.details.csv" file.}
 \item{lassoFile}{a character string naming the file where LASSO model details will be stored. Default value is the "Lasso.details.csv" file.}  
 \item{svrmFile}{a character string naming the file where SVMR model details will be stored. Default value is the "SVMRadial.details.csv" file.}
 \item{nnFile}{a character string naming the file where NN model details will be stored. Default value is the "NN.details.csv" file.}
 \item{rfFile}{a character string naming the file where RF model details will be stored. Default value is the "RF.details.csv" file.}
 \item{rfrfeFile}{a character string naming the file where RFRFE model details will be stored. Default value is the "RFRFE.details.csv" file.}
 \item{svmrfeFile}{a character string naming the file where SVMRRFE model details will be stored. Default value is the "SVMRFE.details.csv" file.}
 \item{enetFile}{a character string naming the file where ENET model details will be stored. Default value is the "ENET.details.csv" file.}
}
\details{
RMSE is the summary metric used to select the optimal model.

To control the computational nuances of the train function, trainControl is used; number of folds or resampling iterations is set to 10, and the number of completed set of folds is set to 10 (for repeated k-fold cross-validation).
    
sCV can take the following values: boot, boot632, cv, repeatedcv, LOOCV, LGOCV (for repeated training/test splits), none (only fits one model to the entire training set), oob (only for random forest, bagged trees, bagged earth, bagged flexible discriminant analysis, or conditional tree forest models), "adaptive_cv", "adaptive_boot" or "adaptive_LGOCV".

If fDet=TRUE, the following output is produced: a CSV file with detailed statistics about the regression model (Regression method, splitting number, cross-validation type, Training set summary, Test set summary, Fitting summary, List of predictors, Training predictors, Test predictors, resampling statistics, features'importance, residuals of the fitted model, assessment of applicability domain (leverage analysis, Cook's distances, points influence)), 5-12 plots for fitting statistics as a PDF file for each splitting and cross-validation method (Training Yobs-Ypred, Test Yobs-Ypred, Feature Importance, Fitted vs. Residuals for Fitted Model, Leverage for Fitted Model, Cook's Distance for Fitted Model, 6 standard fitting plots using plot function with cutoff.Cook). 

}
\value{
 A list is returned containing:
\item{stat.values}{model's statistics}
\item{model}{the full stepwise-selected model, i.e. a list of class train}
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

outGLM<- 'GLMoutput.csv'

GLM.fit <- GLMreg(ds.train,ds.test,CVtypes[1],iSplit=1,fDet=F,outFile=outGLM)
}
}
\author{
Georgia Tsiliki, Cristian R. Munteanu
}
