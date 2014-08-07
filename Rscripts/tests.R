# TESTING ANY CODE

library(caret)

r2.adj.funct<- function(y,y.new,num.pred){#y==y, y.new=predicted, num.pred=number of idependent variables (predictors)
  y.mean<- mean(y)
  x.in<- sum((y-y.new)^2)/sum((y-y.mean)^2)
  x.in<- 1-x.in #r squared
  
  x.in<- (1-x.in)*((length(y)-1)/(length(y)-num.pred-1))
  x.in<- 1 - x.in 
  return(x.in)
}
#----------------------------------------------------------------------------------------------------------------------
r2.adj.lm.funct<- function(y,y.new,num.pred){ #y==y, y.new=predicted, num.pred=number of idependent variables (predictors)
  #only for lm with intercept
  #y.mean<- mean(y)
  #x.in<- sum((y-y.new)^2)/sum((y-y.mean)^2)
  #x.in<- 1-x.in #r squared
  x.in<- cor(y,y.new)^2
  x.in<- (1-x.in)*((length(y)-1)/(length(y)-num.pred-1))
  x.in<- 1 - x.in 
  return(x.in)
}
#----------------------------------------------------------------------------------------------------------------------
rmse.funct<- function(y,y.new){               #y==y, y.new=predicted
  return(sqrt(mean((y.new - y)^2)))
}
#----------------------------------------------------------------------------------------------------------------------
r2.funct<- function(y,y.new){                 #y==y, y.new=predicted
  y.mean<- mean(y)
  x.in<- sum((y-y.new)^2)/sum((y-y.mean)^2)
  x.in<- 1-x.in #r squared
  return(x.in)
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

sCV = "repeatedcv"
iSplit=1

attach(my.datf.train)   # make available the names of variables from training dataset
RegrMethod <- "GLM.AIC" # type of regression

# Define the CV conditions
ctrl<- trainControl(method = sCV, number = 10,repeats = 10,
                    summaryFunction = defaultSummary)

# Train the model using only training set
set.seed(iSplit)
lm.fit<- train(net.c~.,data=my.datf.train,
               method = 'glmStepAIC', tuneLength = 10, trControl = ctrl,
               metric = 'RMSE')

#------------------------------
# Training RESULTS
#------------------------------
RMSE.tr  <- lm.fit$results[,2]
R2.tr    <- lm.fit$results[,3]
if (sCV == "repeatedcv"){ # if 10-fold CV
  RMSEsd.tr <- lm.fit$results[,4]
  R2sd.tr   <- lm.fit$results[,5]
}
if (sCV == "LOOCV"){ # if LOOCV
  RMSEsd.tr <- 0 # formulas will be added later!
  R2sd.tr   <- 0 # formulas will be added later!
}

#------------------------------------------------
# RMSE & R^2, for train/test respectively
#------------------------------------------------
lm.train.res <- getTrainPerf(lm.fit)
lm.test.res  <- postResample(predict(lm.fit,my.datf.test),my.datf.test[,1])

#------------------------------------------------
# Adj R2, Pearson correlation
#------------------------------------------------
pred.tr     <- predict(lm.fit,my.datf.train) # predicted Y
pred.ts     <- predict(lm.fit,my.datf.test)  # predicted Y
noFeats.fit <- length(predictors(lm.fit))    # no. of features from the fitted model
Feats.fit   <- paste(predictors(lm.fit),collapse="+") # string with the features included in the fitted model

ds.full     <- rbind(my.datf.train,my.datf.test)
pred.both   <- predict(lm.fit,ds.full)       # predicted Y
adjR2.tr    <- r2.adj.funct(my.datf.train[,1],pred.tr,noFeats.fit)
adjR2.ts    <- r2.adj.funct(my.datf.test[,1],pred.ts,noFeats.fit)
corP.ts     <- cor(my.datf.test[,1],pred.ts)

adjR2.both  <- r2.adj.funct(ds.full[,1],pred.both,noFeats.fit)
RMSE.both   <- rmse.funct(ds.full[,1],pred.both)
r2.both     <- r2.funct(ds.full[,1],pred.both)

# Generate the output list with statistics for each cross-validation type
# --------------------------------------------------------------------
# There are 3 lists for Data set info, 10-fold CV and LOOCV
# These 3 lists will be merged into one in order obtain the function output including the header (statistics names)

my.stats <- list("RegrMeth"     = RegrMethod,
                 "Split No"     = as.numeric(iSplit),     # from function param
                 "CVtype"       = sCV,                    # from function param
                 "NoModelFeats" = as.numeric(noFeats.fit),
                 "ModelFeats"   = Feats.fit,
                 "adjR2.tr"  = as.numeric(adjR2.tr),
                 "RMSE.tr"   = as.numeric(RMSE.tr),
                 "R2.tr"     = as.numeric(R2.tr),
                 "RMSEsd.tr" = as.numeric(RMSEsd.tr),
                 "R2sd.tr"   = as.numeric(R2sd.tr),
                 "adjR2.ts"= as.numeric(adjR2.ts),
                 "RMSE.ts" = as.numeric((lm.test.res["RMSE"][[1]])),
                 "R2.ts"   = as.numeric((lm.test.res["Rsquared"][[1]])),
                 "corP.ts" = as.numeric(corP.ts),
                 "adjR2.both" = as.numeric(adjR2.both),
                 "RMSE.both"  = as.numeric(RMSE.both),
                 "R2.both"    = as.numeric(r2.both))

#---------------------------------------------------------------------
# Write to file DETAILS for GLM for each cross-validation method
#---------------------------------------------------------------------
if (fDet==TRUE) {   # if flag for details if true, print details about any resut
  write("RRegr package | eNanoMapper", file = outFile)
  write.table(paste("Regression method: ", RegrMethod), file = outFile,append = TRUE, sep = " ",col.names = FALSE,quote = FALSE)
  write.table(paste("Split no.: ", iSplit), file = outFile,append = TRUE, sep = " ",col.names = FALSE,quote = FALSE)
  write.table(paste("CV type: ", sCV), file = outFile,append = TRUE, sep = " ",col.names = FALSE,quote = FALSE)
  write.table("Training Set Summary: ", file = outFile,append = TRUE, sep = " ",col.names = FALSE,quote = FALSE)
  write.table(summary(my.datf.train), file = outFile,append = TRUE, sep = " ",col.names = TRUE,quote = FALSE)
  write.table("Test Set Summary: ", file = outFile,append = TRUE, sep = " ",col.names = FALSE,quote = FALSE)
  write.table(summary(my.datf.test), file = outFile,append = TRUE, sep = " ",col.names = TRUE,quote = FALSE)   
  
  write.table("Fitting Summary: ", file = outFile,append = TRUE, sep = " ",col.names = FALSE,quote = FALSE)
  write.table(data.frame(summary(lm.fit)$coefficients), file = outFile,append = TRUE, sep = " ",col.names = TRUE,quote = FALSE)
  
  write.table("Predictors: ", file = outFile,append = TRUE, sep = " ",col.names = FALSE,quote = FALSE)
  write.table(predictors(lm.fit), file = outFile,append = TRUE, sep = " ",col.names = TRUE,quote = FALSE)
  
  write.table("Trainig Results: ", file = outFile,append = TRUE, sep = " ",col.names = FALSE,quote = FALSE)
  write.table(predictors(lm.train.res), file = outFile,append = TRUE, sep = " ",col.names = TRUE,quote = FALSE)
  write.table("Test Results: ", file = outFile,append = TRUE, sep = " ",col.names = FALSE,quote = FALSE)
  write.table(predictors(lm.test.res), file = outFile,append = TRUE, sep = " ",col.names = TRUE,quote = FALSE)
  
  write.table("Full Statistics: ", file = outFile,append = TRUE, sep = " ",col.names = FALSE,quote = FALSE)
  write.table(my.stats, file = outFile,append = TRUE, sep = " ",col.names = TRUE,quote = FALSE)
}

#####################################################
print(data.frame(my.stats))
print(length(my.stats))
