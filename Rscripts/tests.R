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
#--------------------------------------------------------------------
# Write a LIST to CSV file
#--------------------------------------------------------------------
AppendList2CSv <- function(l,csvFile) {
  out_file <- file(csvFile, open="a")  #creates a file in append mode 
  for (i in seq_along(l)){ 
    write.table(names(l)[i], file=out_file, sep=",", dec=".", 
                quote=FALSE, col.names=FALSE, row.names=FALSE)  #writes the name of the list elements ("A", "B", etc) 
    write.table(l[[i]], file=out_file, sep=",", dec=".", quote=FALSE, 
                col.names=NA, row.names=TRUE)  #writes the data.frames 
  } 
  close(out_file)  #close connection to file.csv 
}  



fDet         = TRUE  # flag to calculate and print details for all the functions
fFilters     = TRUE  # flag to apply filters                          (2)
fScaling     = TRUE  # flag for dataset Scaling                       (3)
fRemNear0Var = TRUE  # flag for Removal of near zero variance columns (4)
fRemCorr     = TRUE  # flag for Removal of correlated columns         (5)
fFeatureSel  = FALSE  # flag for wrapper methods for feature selection (7)

cutoff       = 0.9   # cut off for correlated features
fGLM         = FALSE  # flag to run GLM (8.2)
fLASSO       = FALSE  # flag to run LASSO (8.4) 
fSVLM        = FALSE # flat to run svmRadial.RMSE (8.6)
fNN          = TRUE  # flat to run NN (8.6)

# ----------------------------------------------------------------------------------------
iScaling = 1 # 1 = normalization; 2 = standardization, 3 = other; any other: no scaling
iScalCol = 1 # 1 = including dependent variable in scaling; 2: only all features; etc.
# ----------------------------------------------------------------------------------------
trainFrac  = 3/4 # the fraction of training set from the entire dataset; trainFrac = the rest of dataset, the test set
iSplitTimes = 2 # default is 10; time to split the data in train and test (steps 6-11); report each step + average

CVtypes    <- c("repeatedcv","LOOCV")             # cross-validation types: 10-CV and LOOCV

# -------------------------------------------------------------------------------------------------------
# Files
# -------------------------------------------------------------------------------------------------------
PathDataSet    = "DataResults"            # dataset folder for input and output files
DataFileName   = "ds.csv"                 # input step 1 = ds original file name
No0NearVarFile = "ds3.No0Var.csv"         # output step 3 = ds without zero near vars
ScaledFile     = "ds4.scaled.csv"         # output step 4 = scaled ds file name (in the same folder)
NoCorrFile     = "ds5.scaled.NoCorrs.csv" # output step 5 = dataset after correction removal

ResAvgs        = "RRegsRes.csv"           # the main output file with averaged statistics for each regression method
ResBySplits    = "RRegrsResBySplit.csv"   # the output file with statistics for each split and the averaged values
glmFile        = "8.2.GLM.details.txt"    # GLM output file for details

# Generate path + file name = original dataset
inFile <- file.path(PathDataSet, DataFileName)

cat("======================================================================
    RRegrs - R Regression Models
    Get the best regression models for one dataset using R caret methods
    NTUA and UM groups, enanomapper.net
    
    Contacts:
    Cristian R Munteanu - muntisa [at] gmail [dot] com
    Georgia Tsiliki - g_tsiliki [at] hotmail [dot] com
    ======================================================================\n")

# -----------------------------------
# (1.2) Load the ORIGINAL DATASET
# -----------------------------------
cat("-> [1] Loading original dataset ...\n")
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

#========================================================
# (2) FILTERS
#     (it will be implemented in the future versions)
#========================================================
# 2.1 Outlier removal
# 2.2 Custom filter (percentage threshold)
# 2.3 Processing of missing values - use of preProcess();
#     caret employs knnImpute algorithm to impute values from a neighborhood of k
cat("-> [2] Filtering dataset ... No filter!\n")

# -----------------------------------------------------------------------
# (3) Remove near zero variance columns
# -----------------------------------------------------------------------
if (fRemNear0Var==TRUE) {
  cat("-> [3] Removal of near zero variance columns ...\n")
  outFile <- file.path(PathDataSet,No0NearVarFile) # the same folder as input  
  
  # get the ds without near zero cols 
  source("s3.RemNearZeroVar.R")                    # add function         
  ds <- cbind("net.c" = ds[,1],RemNear0VarCols(ds[,2:dim(ds)[2]],fDet,outFile))
  # use df without Y (predicted values), reconstruct the ds
  # inputs: ds, flag for details, output file
}

# -----------------------------------------------------------------------
# (4) Scaling dataset: normalization (default), standardization, other
# -----------------------------------------------------------------------
if (fScaling==TRUE) {
  cat("-> [4] Scaling original dataset ...\n")
  outFile <- file.path(PathDataSet,ScaledFile)       # the same folder as input
  
  # run fuction for scaling input dataset file
  source("s4.ScalingDataSet.R")                      # add function
  ds <- ScalingDS(ds,iScaling,iScalCol,fDet,outFile)
  # use df without Y (predicted values), reconstruct the ds
  # inputs: ds, type of scaling, flag for details, starting column, output file
}

# -----------------------------------------------------------------------
# (5) Remove correlated features
# -----------------------------------------------------------------------

library(corrplot) #corrplot: the library to compute correlation matrix.
library(caret)
DataSet <- ds[,2:dim(ds)[2]]               # input dataset
DataSetFiltered.scale <- ds[,2:dim(ds)[2]] # default results without any modification

# calculate the correlation matrix for the entire file!
# !!! NEED TO BE CORRECTED to avoid dependent variable (first column) but to report it!
corrMat <- cor(DataSet)                                             # get corralation matrix


highlyCor <- findCorrelation(corrMat, cutoff) # find corralated columns

# Apply correlation filter with the cutoff
# by removing all the variable correlated with more than cutoff
DataSetFiltered.scale <- DataSet[,-highlyCor]

if (fDet==TRUE) {
  corrMat <- cor(DataSetFiltered.scale)
  # plot again the rest of correlations after removing the correlated columns
  #corrplot(corrMat, order = "hclust")
  

  # correlation matrix for the rest of the columns after removal
  CorrMatFile2 <- paste(outFile,".corrMAT4Selected.csv",sep='')
  # write correlation matrix as output file
  write.csv(corrMat, CorrMatFile2, row.names=FALSE, quote=FALSE)
  # write the new dataset without the correlated features
  write.csv(DataSetFiltered.scale, outFile, row.names=FALSE, quote=FALSE)
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

sCV = "repeatedcv"
iSplit=1

####################
# TO TEST
####################

library(caret)

net.c = my.datf.train[,1]   # make available the names of variables from training dataset
RegrMethod <- "pls" # type of regression

# Define the CV conditions
ctrl<- trainControl(method = sCV, number = 10,repeats = 10,
                    summaryFunction = defaultSummary)

# Train the model using only training set
set.seed(iSplit)

pls.fit<- train(net.c~.,data=my.datf.train,
                method = 'pls', tuneLength = 10, trControl = ctrl,
                metric = 'RMSE',
                tuneGrid=expand.grid(.ncomp=c(1:(dim(my.datf.train)[2]-1))))
#------------------------------
# Training RESULTS
#------------------------------
RMSE.tr  <- pls.fit$results[,2]
R2.tr    <- pls.fit$results[,3]
if (sCV == "repeatedcv"){ # if 10-fold CV
  RMSEsd.tr <- pls.fit$results[,4]
  R2sd.tr   <- pls.fit$results[,5]
}
if (sCV == "LOOCV"){ # if LOOCV
  RMSEsd.tr <- 0 # formulas will be added later!
  R2sd.tr   <- 0 # formulas will be added later!
}

#------------------------------------------------
# RMSE & R^2, for train/test respectively
#------------------------------------------------
lm.train.res <- getTrainPerf(pls.fit)
lm.test.res  <- postResample(predict(pls.fit,my.datf.test),my.datf.test[,1])

#------------------------------------------------
# Adj R2, Pearson correlation
#------------------------------------------------
pred.tr     <- predict(pls.fit,my.datf.train) # predicted Y
pred.ts     <- predict(pls.fit,my.datf.test)  # predicted Y
noFeats.fit <- length(predictors(pls.fit))    # no. of features from the fitted model
Feats.fit   <- paste(predictors(pls.fit),collapse="+") # string with the features included in the fitted model

ds.full     <- rbind(my.datf.train,my.datf.test)
pred.both   <- predict(pls.fit,ds.full)       # predicted Y
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
                 "CVtype"       = "no CV",                # from function param
                 "NoModelFeats" = as.numeric(noFeats.fit),
                 "ModelFeats"   = Feats.fit,
                 "adjR2.tr"  = as.numeric(adjR2.tr),
                 
                 "RMSE.tr"   = RMSE.tr,  # these 4 lines correspond to the min of RMSE.tr !!!
                 "R2.tr"     = R2.tr,  
                 "RMSEsd.tr" = RMSEsd.tr,
                 "R2sd.tr"   = R2sd.tr,
                 
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
  
  
  write.table("Predictors: ", file = outFile,append = TRUE, sep = " ",col.names = FALSE,quote = FALSE)
  write.table(predictors(pls.fit), file = outFile,append = TRUE, sep = " ",col.names = TRUE,quote = FALSE)
  
  write.table("Trainig Results: ", file = outFile,append = TRUE, sep = " ",col.names = FALSE,quote = FALSE)
  write.table(predictors(lm.train.res), file = outFile,append = TRUE, sep = " ",col.names = TRUE,quote = FALSE)
  write.table("Test Results: ", file = outFile,append = TRUE, sep = " ",col.names = FALSE,quote = FALSE)
  write.table(predictors(lm.test.res), file = outFile,append = TRUE, sep = " ",col.names = TRUE,quote = FALSE)
  
  write.table("Full Statistics: ", file = outFile,append = TRUE, sep = " ",col.names = FALSE,quote = FALSE)
  write.table(my.stats, file = outFile,append = TRUE, sep = " ",col.names = TRUE,quote = FALSE)
}


#####################################################
dfRes = data.frame(my.stats)
print("*** RESULTS ***")
print(dfRes)
print(length(my.stats))
# print(dfRes[which.min(dfRes$RMSE.tr),])
