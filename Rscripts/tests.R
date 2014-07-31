# ======================================================================
# RRegrs - R Regression Models
#
# Best regression models for one dataset using R methods
# Developed as tool for nano-toxicity QSAR models
# NTUA and UM groups
# contact: Cristian R Munteanu | BiGCaT - UM    | muntisa@gmail.com
#          Georgia Tsiliki     | ChemEng - NTUA | g_tsiliki@hotmail.com
# ======================================================================
# Main input file: CSV
# ----------------------------------------------------------------------
# Variable names
# -----------------------------------------------------------------------
# DataSet = dataset including header, Y and X features values
# Yname, Xnames = names of dependent and features
# NoX = no of features
# NoCases = no of cases
# YData = values of Y
# XData = values of the features
# -----------------------------------------------------------------------
# Note:
# - Input file = CSV format and extension:
#    - there is only one input file with the original dataset:
#      Dependent Variable, Feature 1, ..., Feature n
#    - the names of the variables are located in the first row
# - the main script is importing separate file functions
# - the functions are reading datafiles in specific folders
# and it is writing the output as variables or datafiles
# - it generates several dataset files for each step
# - it generates PNG plots for the correlation matrix
# before and after the correlation removal
# - all the input and output files are placed into the same folder


#
#Methods included are: pls, lasso, lm, lmAIC, svmLinear, svmRadial, neural network (single hidden layer),
#som, rbf.
#
#For those which do not include feature selection, rfe() function has been used (stands for recursive feature selection).
#
#REsampling: 10fold CV (no repeats), boostrap for rfe using the default 25 samples
#
#adjsuted R2 and RMSE functions are included at the top.
#

r2.adj.funct<- function(y,y.new,num.pred){#y==y, y.new=predicted, num.pred=number of idependent variables (predictors)
  y.mean<- mean(y)
  x.in<- sum((y-y.new)^2)/sum((y-y.mean)^2)
  x.in<- 1-x.in #r squared
  
  x.in<- (1-x.in)*((length(y)-1)/(length(y)-num.pred-1))
  x.in<- 1 - x.in 
  return(x.in)
}

r2.adj.lm.funct<- function(y,y.new,num.pred){#y==y, y.new=predicted, num.pred=number of idependent variables (predictors)
  #only for lm with intercept
  #y.mean<- mean(y)
  #x.in<- sum((y-y.new)^2)/sum((y-y.mean)^2)
  #x.in<- 1-x.in #r squared
  x.in<- cor(y,y.new)^2
  x.in<- (1-x.in)*((length(y)-1)/(length(y)-num.pred-1))
  x.in<- 1 - x.in 
  return(x.in)
}


rmse.funct<- function(y,y.new){
  return(sqrt(mean((y.new - y)^2)))
}


#==========================================================================================
# (1) Load dataset and parameters
#     (these parameters will be read from an input file! TO BE IMPLEMENTED at the end)
#==========================================================================================
# (1.1) PARAMETERS
# -----------------------------------------------------------------------
# Option to run any step
# -----------------------------------------------------------------------
fDet = TRUE          # flag to calculate and print details for all the functions
fFilters=TRUE        # flag to apply filters                          (2)
fScaling=TRUE        # flag for dataset Scaling                       (3)
fRemNear0Var=TRUE    # flag for Removal of near zero variance columns (4)
fRemCorr=TRUE        # flag for Removal of correlated columns         (5)
fFeatureSel=TRUE     # flag for selection before the regression       (7)
cutoff=0.7           # cut off for corralations
# ----------------------------------------------------------------------------------------
iScaling = 1 # 1 = standardization; 2 = normalization, 3 = other; any other: no scaling
iScalCol = 1 # 1: including dependent variable in scaling; 2: only all features; etc.
# ----------------------------------------------------------------------------------------
trainFrac  = 3/4 # the fraction of training set from the entire dataset
#                # 1 - trainFrac = the rest of dataset, the test set
# ---------------------------------------------------------------------------------------------------
# Files
# ---------------------------------------------------------------------------------------------------
PathDataSet    = "DataResults"            # dataset folder
DataFileName   = "ds.csv"                 # input step 1 = ds original file name
No0NearVarFile = "ds3.No0Var.csv"         # output step 3 = ds without zero near vars
ScaledFile     = "ds4.scaled.csv"         # output step 4 = scaled ds file name (in the same folder)
NoCorrFile     = "ds5.scaled.NoCorrs.csv" # output step 5 = dataset after correction removal
ResFile        = "RRegresRezults.txt"     # the common output file! 

glmFile        = "glm.res.txt"

# Generate path + file name = original dataset
inFile <- file.path(PathDataSet, DataFileName)
# Set the file with results (append data using: sink(outRegrFile, append = TRUE) !!!)
outRegrFile <- file.path(PathDataSet,ResFile) # the same folder as the input 

# -----------------------------------
# (1.2) Load the ORIGINAL DATASET
# -----------------------------------
# (it can contain errors, correlations, near zero variance columns)
ds.dat0 <- read.csv(inFile,header=T)                              # original dataset frame

# resolving the text to number errors for future calculations
ds.indx<- colnames(ds.dat0)[2:dim(ds.dat0)[2]]                    # FEATURE names (no dependent variable)
ds.dat1<- ds.dat0[1:dim(ds.dat0)[1],2:dim(ds.dat0)[2]]            # dataset as columns
ds.dat1<- apply(ds.dat1,1,function(x)as.numeric(as.character(x))) # dataset as row vectors to be used with caret!!!

# dependent variable
net.c<- ds.dat0[,1]
net.c<- as.numeric(as.character(net.c)) # values
# full ds frame with training and test
ds<- as.data.frame(cbind(net.c,t(ds.dat1)))


# -----------------------------------------------------------------------
# (3) Remove near zero variance columns
# -----------------------------------------------------------------------
if (fRemNear0Var==TRUE) {
  print("-> [3] Removal of near zero variance columns ...")
  outFile <- file.path(PathDataSet,No0NearVarFile) # the same folder as input  
  
  # get the ds without near zero cols 
  source("s3.RemNearZeroVar.R")                    # add function
  ds <- RemNear0VarCols(ds,fDet,outFile)           # inputs: ds, flag for details, output file
}

# -----------------------------------------------------------------------
# (6) Dataset split: Training and Test sets
# -----------------------------------------------------------------------
print("-> [6] Splitting dataset in Training and Test sets ...")

source("s6.DsSplit.R")  # add external function
iSeed=1                 # to reapeat the ds splitting, different values of seed will be used
dsList <- DsSplit(ds,trainFrac,fDet,PathDataSet,iSeed) # return a list with 2 datasets = dsList$train, dsList$test
# get train and test from the resulted list
my.datf.train <- dsList$train
my.datf.test  <- dsList$test

my.datf.train
my.datf.test


####################
# TO TEST
####################

my.stats<- list() # create empty result
# specify CV parameters
ctrl<- trainControl(method = 'repeatedcv', number = 10,repeats = 10,
                    summaryFunction = defaultSummary)

# Training the model using only training set
set.seed(2)
attach(my.datf.train)
lm.fit<- train(net.c~.,data=my.datf.train,
               method = 'glmStepAIC', tuneLength = 10, trControl = ctrl,
               metric = 'RMSE')

# Training RESULTS
#------------------------------
# get statistics for training
RMSE = lm.fit$results[,2]
Rsquared = lm.fit$results[,3]
RMSE_SD = lm.fit$results[,4]
Rsquared_SD = lm.fit$results[,5]

# RMSE & R^2, for train/ test respectively
lm.train.res<- getTrainPerf(lm.fit)
lm.test.res <- postResample(predict(lm.fit,my.datf.test),my.datf.test[,1])

# Adj R2: y obs, y predicted, No of predictors
ds.full = rbind(my.datf.train,my.datf.test)
adjR2_both = r2.adj.funct(ds.full[,1], predict(lm.fit,ds.full), length(predictors(lm.fit)))
adjR2_train= r2.adj.funct(my.datf.train[,1],predict(lm.fit,my.datf.train),length(predictors(lm.fit)))
adjR2_test = r2.adj.funct(my.datf.test[,1], predict(lm.fit,my.datf.test), length(predictors(lm.fit)))
  
if (fDet==TRUE) {
  # write RESULTS
  sink(outFile)
  print(summary(my.datf.train))
  print(summary(my.datf.test))
  print(lm.fit)
  print(predictors(lm.fit))
  print(lm.train.res)
  print(lm.test.res)
  sink()
  #file.show(outFile)
}
my.stats = list("Regrrs"="GLM","Step"=1,"tr.RMSE.10CV"= RMSE,"tr.R2.10CV" = Rsquared,
                "tr.RMSESD.10CV" = RMSE_SD,"tr.R2SD.10CV" = Rsquared_SD,
                "ts.RMSE" = (lm.test.res["RMSE"][[1]]),"ts.R2" = (lm.test.res["Rsquared"][[1]]),
                "both.adjR2" = adjR2_both, "tr.adjR2" = adjR2_train, "ts.adjR2" = adjR2_test)
# default values: "Regrrs"="GLM","Step"=1 (or for one single use of funtion)
print(data.frame(my.stats))