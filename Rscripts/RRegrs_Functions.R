# ======================================================================
# RRegrs functions based on caret package
# ----------------------------------------------------------------------
# eNanoMapper: enanomapper.net
# contact: Cristian R Munteanu | BiGCaT - UM    | muntisa@gmail.com
#          Georgia Tsiliki     | ChemEng - NTUA | g_tsiliki@hotmail.com


# library(caret)

#======================================================================================================================
# General functions
#======================================================================================================================
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
    # writes the name of the list elements ("A", "B", etc.)
    write.table(names(l)[i],file=out_file,sep=",",dec=".",quote=F,col.names=F,row.names=F) 
    write.table(l[[i]],     file=out_file,sep=",",dec=".",quote=F,col.names=NA,row.names=T) #writes the data.frames 
  } 
  close(out_file)  #close connection to file.csv 
}  
#--------------------------------------------------------------------
# Write a LIST to TXT file
#--------------------------------------------------------------------
AppendList2txt <- function(l,csvFile) {
  out_file <- file(csvFile, open="a")  #creates a file in append mode 
  for (i in seq_along(l)){ 
    #writes the name of the list elements ("A", "B", etc) 
    write.table(names(l)[i],file=out_file,sep=" ",dec=".",quote=F,col.names=F, row.names=F)
    write.table(l[[i]],     file=out_file,sep=" ",dec=".",quote=F,col.names=NA,row.names=T) #writes the data.frames 
  } 
  close(out_file)  #close connection to file.csv 
}  

# ************************************
# RRegrs Specific functions
# ************************************

#=======================================================================
# Removal of near zero variance columns (Step 3)
#=======================================================================
# inputs:
# - ds = dataset frame
# - fDet = flag for detais (TRUE/FALSE)
# - outFileName = new file name  (it could include the path)

# output = ds.Rem0NearVar  (ds without columns with near zero variance)
# if datails = TRUE, output the new ds as a file
# ------------------------------------------

RemNear0VarCols <- function(ds,fDet=FALSE,outFile="ds3.No0Var.csv") {
  # default parameters are no details, with a CSV file name
  library(caret)
  ds.Rem0NearVar <- ds               # default output without any modification
  ds.var <- nearZeroVar(ds)          # get the near zero columns
  if (!length(ds.var) == FALSE) {    # remove the columns only if nearZeroVar identified; if no columns to remove, ds will be the same
    ds.Rem0NearVar <- ds[,-(ds.var)] # get only the columns without this problem
    if (fDet == TRUE) {              # write as details the corrected ds file
      write.csv(ds.Rem0NearVar, outFile,row.names=F, quote=F)
    }
  }
  return(as.data.frame(ds.Rem0NearVar)) # return the new data frame without near zero variance
}

#===================================================================================
# Scaling dataset (Step 4)
#===================================================================================
# s = { 1,2,3 } - type of scaling: 1 = normalization, 2 = standardization, 3 = other
# c = the number of column into the dataset to start scaling
# - if c = 1: included the dependent variable
# - if c = 2: only the features will be scaled

# fDet = if details need to be printed    (TRUE/FALSE)
# inFileName  = file name      (it could include the path)
# outFileName = new file name  (it could include the path)
#-----------------------------------------------------------------------------------

ScalingDS <- function(ds,s=1,c=1,fDet=FALSE,outFileName="ds4.scaled.csv") {
  # Default scaling = NORMALIZATION !
  # if we use preProcess() from caret, errors are generating if threre are zero variance columns!
  
  # DEFAULT scaled dataset = original
  # if other s diffent of 1,2,3 is used, no scaling!
  DataSet.scaled <- ds
  
  # if NORMALIZATION
  if (s==1) {
    # Scale all the features (from column c; column 1 is the predictor output)
    DataSet.scaled <- ((ds-min(ds))/(max(ds)-min(ds))) # normalize all the columns
  }
  # if STADARDIZATION
  if (s==2) {
    # Scale all the features (from column c; column 1 is the predictor output)
    DataSet.scaled <- scale(ds[c:ncol(ds)],center=TRUE,scale=TRUE)  
  }
  
  # if other scaling
  if (s==3) {
    # Scale all the features (from feature 2 bacause feature 1 is the predictor output)
    # TO ADD THE CODE ! 
  }
  
  # if DETAILS
  if (fDet ==TRUE) {
    # write the result into a separated file
    write.csv(DataSet.scaled, outFileName,row.names=F, quote=F)  
  }
  return (as.data.frame(DataSet.scaled)) # return the scaled data frame
}

# =======================================================================
# Remove the correlated columns (Step 5)
# =======================================================================
# ds = dataset frame
# fDet = flag fro details (TRUE/FALSE)
# cutoff = correlation cut off (ex: 0.9)
# outFileName = new file name  (it could include the path)

# Generates 5 file results:
# - returns a dataset without the correlated columns  (1 file)
# - generate initial correlation matrix 
#   and the one after removing the correlated features (2 files)
# - plots for the before and after correlation removal (2 files)
# ------------------------------------------------------------------------

# another version of this function should be implemented using
# pairwise test between i and j descriptors- if(r2>=0.9){remove the j descriptor}
# using findCorelations() from caret

RemCorrs <- function(ds,fDet,cutoff,outFile) {
  
  library(corrplot) #corrplot: the library to compute correlation matrix.
  library(caret)
  
  DataSet <- ds               # input dataset
  DataSetFiltered.scale <- ds # default results without any modification
  
  # calculate the correlation matrix for the entire file!
  # !!! NEED TO BE CORRECTED to avoid dependent variable (first column) but to report it!
  corrMat <- cor(DataSet)                                             # get corralation matrix
  
  if (fDet==TRUE) {
    CorrMatFile <- paste(outFile,".corrMAT.csv",sep='')
    # write correlation matrix as output file
    write.csv(corrMat, CorrMatFile, row.names=F, quote=F)
    
    # Plot the matrix, clustering features by correlation index
    # corrplot(corrMat, order = "hclust")
    
    # plot the correlatio plot before correlation removal
    CorrPlotFile <-  paste(outFile,".corrs.png",sep='')
    png(height=1200, width=1200, pointsize=25, file=CorrPlotFile)
    col1 <-rainbow(100, s = 1, v = 1, start = 0, end = 0.9, alpha = 1)
    corrplot(corrMat,tl.cex=3,title="Initial feature correlation matrix",
             method="circle",is.corr=FALSE,type="full",
             cl.lim=c(-1,1),cl.cex=2,addgrid.col="red",
             addshade="positive",col=col1,
             addCoef.col = rgb(0,0,0, alpha = 0.6), mar=c(0,0,1,0), diag= FALSE) 
    dev.off()
  }
  
  highlyCor <- findCorrelation(corrMat, cutoff) # find corralated columns
  
  # Apply correlation filter with the cutoff
  # by removing all the variable correlated with more than cutoff
  DataSetFiltered.scale <- DataSet[,-highlyCor]
  
  if (fDet==TRUE) {
    corrMat <- cor(DataSetFiltered.scale)
    # plot again the rest of correlations after removing the correlated columns
    #corrplot(corrMat, order = "hclust")
    
    # plot the correlatio plot before correlation removal
    CorrPlotFile2 =  paste(outFile,".afterRemCorr.png",sep='')
    png(height=1200, width=1200, pointsize=25, file=CorrPlotFile2)
    col1 <-rainbow(100, s = 1, v = 1, start = 0, end = 0.9, alpha = 1)
    corrplot(corrMat,tl.cex=3,title="Correlation matrix after removing correlated features",
             method="circle",is.corr=FALSE,type="full",
             cl.lim=c(-1,1),cl.cex=2,addgrid.col="red",
             addshade="positive",col=col1,
             addCoef.col = rgb(0,0,0, alpha = 0.6), mar=c(0,0,1,0), diag= FALSE) 
    dev.off()
    # correlation matrix for the rest of the columns after removal
    CorrMatFile2 <- paste(outFile,".corrMAT4Selected.csv",sep='')
    # write correlation matrix as output file
    write.csv(corrMat, CorrMatFile2, row.names=F, quote=F)
    # write the new dataset without the correlated features
    write.csv(DataSetFiltered.scale, outFile, row.names=F, quote=F)
  }
  
  return(as.data.frame(DataSetFiltered.scale))
}

# ======================================================================
# Dataset spliting in Training and Test (Step 6)
# ======================================================================
# Inputs
# - ds = frame dataset object
# - fDet = flag for detais (TRUE/FALSE)
# - PathDataSet = pathway for results

# Output = training and test datasets (to be used for regressions in other functions)
# if datails = TRUE, output files will be created

# ------------------------------------------
DsSplit <- function(ds,trainFrac=3/4,fDet=FALSE,PathDataSet="",iSeed) {
  my.datf<- ds
  # create TRAIN and TEST sets to build a model
  set.seed(iSeed)
  inTrain <- createDataPartition(1:dim(my.datf)[1],p = trainFrac,list = FALSE,groups=2)
  # groups==2 forces to NOT partition
  # based on quantiles of numeric values
  my.datf.train<- my.datf[inTrain,]            # TRAIN dataset frame         
  my.datf.test <- my.datf[-inTrain,]           # TEST dataset frame
  
  if (fDet == TRUE) {
    # write the TRAIN and TEST set files
    # the index of each row will in the dataset will not be saved (row.names=F)
    outTrain <- file.path(PathDataSet,"ds.Train.csv") # the same folder as the input
    write.csv(my.datf.train,outTrain,row.names=FALSE)
    outTest <- file.path(PathDataSet,"ds.Test.csv") # the same folder as the input
    write.csv(my.datf.test,outTest,row.names=FALSE) 
  }
  MyList<- list("train"=my.datf.train, "test"=my.datf.test) 
  return(MyList)  # return a list with training and test datasets
}

# *************************************
# REGRESSION METHODS
# *************************************

#====================================================================================================
# 8.1. Basic LM
#====================================================================================================
LMreg <- function(my.datf.train,my.datf.test,sCV,iSplit=1,fDet=F,outFile="") {
  # ------------------------------------------
  # only Linux or Mac: parallel calculations
  # ------------------------------------------
  #library(doMC)
  #registerDoMC(cores = 2) # CPU cores
  # ------------------------------------------
  
  net.c = my.datf.train[,1]   # make available the names of variables from training dataset
  RegrMethod <- "lm" # type of regression
  
  # Define the CV conditions
  ctrl<- trainControl(method=sCV, number=10,repeats=10,
                      summaryFunction=defaultSummary)
  
  # Train the model using only training set
  set.seed(iSplit)
  lm.fit<- train(net.c~.,data=my.datf.train,
                 method='lm', tuneLength = 10,trControl=ctrl,
                 metric='RMSE')
  
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
  pred.tr     <- predict(lm.fit,my.datf.train) # predicted Y for training
  pred.ts     <- predict(lm.fit,my.datf.test)  # predicted Y for test
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
  if (fDet==T) {   # if flag for details if T, print details about any resut
    write("RRegr package | eNanoMapper", file=outFile, append=T)
    write.table(paste("Regression method: ", RegrMethod), file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(paste("Split no.: ", iSplit), file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(paste("CV type: ", sCV),      file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table("Training Set Summary: ",     file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(summary(my.datf.train),       file=outFile,append=T,sep=",",col.names=T,quote=F)
    write.table("Test Set Summary: ",         file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(summary(my.datf.test),        file=outFile,append=T,sep=",",col.names=T,quote=F)   
    
    write.table("Fitting Summary: ",                      file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(data.frame(summary(lm.fit)$coefficients), file=outFile,append=T,sep=",",col.names=T,quote=F)
    
    write.table("Predictors: ",     file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(predictors(lm.fit), file=outFile,append=T,sep=",",col.names=T,quote=F)
    
    write.table("Trainig Results: ",     file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(predictors(lm.train.res),file=outFile,append=T,sep=",",col.names=T,quote=F)
    write.table("Test Results: ",        file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(predictors(lm.test.res), file=outFile,append=T,sep=",",col.names=T,quote=F)
    
    write.table("Full Statistics: ", file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(my.stats,            file=outFile,append=T,sep=",",col.names=T,quote=F)
    
    # Variable Importance (max top 20)
    FeatImp <- varImp(lm.fit, scale = F)
    components = length(FeatImp)  # default plot all feature importance
    if (length(FeatImp)>20){     # if the number of features is greater than 20, use only 20
      components = 20
    }
    # Append feature importance to output details
    AppendList2CSv(FeatImp,outFile)
    
    
    fitModel <- lm.fit$finalModel
    
    # =============================================================================
    # Assessment of Applicability Domain (plot leverage)
    # =============================================================================
    
    # Residuals
    resids <- residuals(fitModel) # residuals
    write.table("Residuals of the fitted model: ",file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(data.frame(resids), file=outFile,append=T,sep=",",col.names=T,row.names=T, quote=F) # write residuals
    
    # Leverage / Hat values
    hat.fit <- hatvalues(fitModel)          # hat values
    hat.fit.df <- as.data.frame(hat.fit)    # hat data frame
    hat.mean <- mean(hat.fit)               # mean hat values
    hat.fit.df$warn <- ifelse(hat.fit.df[, 'hat.fit']>3*hat.mean, 'x3',ifelse(hat.fit.df[, 'hat.fit']>2*hat.mean, 'x2', '-' ))
    
    write.table("Leverage output: ",file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(paste("Mean of hat values: ", hat.mean), file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table("Leverage / Hat values with warnings (X3 & X2 = values 3 & 2 times than hat mean): ",file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(hat.fit.df, file=outFile,append=T,sep=",",col.names=T,row.names=T, quote=F) # write hat values and the levels X3, X2 (of hat mean)
    
    #THRESHOLD values: 3m/n, where m is the number of parameters, and n number of observations
    thresh.lever<- (3*(dim(my.datf.train)[2]-1))/dim(my.datf.train)[1] # leverage thresh
    hat.problems<- data.frame(hat.fit[hat.fit>thresh.lever]) # points with high leverage
    
    write.table(paste("Leverage Threshold: ", thresh.lever), file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table("Points with leverage > threshold: ",file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(hat.problems, file=outFile,append=T,sep=",",col.names=T,row.names=T, quote=F)
    
    # Cook's distance
    cook.dists<- cooks.distance(fitModel)
    cutoff.Cook <- 4/((nrow(my.datf.train)-length(fitModel$coefficients)-2)) # Cook's distance cutoff
    
    write.table("Cook's distances output: ",file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(paste("Cook's distance cutoff: ", cutoff.Cook), file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table("Cook's distances: ",file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(data.frame(cook.dists), file=outFile,append=T,sep=",",col.names=T,row.names=T, quote=F) # write residuals
    
    # Influence
    infl <- influence(fitModel)#produces several statistics of the kind
    
    write.table("Point influence output: ",file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table("Influences: ",file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(data.frame(infl), file=outFile,append=T,sep=",",col.names=T,row.names=T, quote=F) # write residuals
    
    # PDF with 12 plots
    # --------------------------------------------------------------
    pdf(file=paste(outFile,".",sCV,".","split",iSplit,".pdf",sep=""))
    # par(mfrow = c(3, 4)) # all plots into one page!
    
    plot(my.datf.train[,1],pred.tr,xlab="Yobs", ylab="Ypred", type="b", main="Train Yobs-Ypred") # plot 1
    plot(my.datf.test[,1], pred.ts,xlab="Yobs", ylab="Ypred", type="b", main="Test Yobs-Ypred")  # plot 2
    dotchart(as.matrix(FeatImp$importance),main="Feature Importance")                            # plot 3
    
    # Fitted vs Residuals - plot 4
    plot(fitted(fitModel),residuals(fitModel),
         main="Fitted vs. Residuals for Fitted Model",
         xlab="Fitted", ylab="Residuals")
    abline(h = 0, lty = 2)
    
    # Leverage plots - plot 5
    plot(hat.fit, type = "h",
         main="Leverage for Fitted Model",
         xlab="Index", ylab="Hat")
    abline(h = thresh.lever, lty = 2, col="red") # leverage thresh
    
    # Cook's distance - plot 6
    plot(cook.dists,
         main="Cook's Distance for Fitted Model",
         xlab="Index", ylab="Cook Distance")
    
    for (p in 1:6) {
      plot(fitModel, which=p, cook.levels=cutoff.Cook) # 6 standard fitting plots
    }
    
    # plot(FeatImp, top = components,main="Feature Importance") # ERROR !
    dev.off()
    # --------------------------------------------------------------

  }
  return(my.stats)  # return a list with statistics
}

#====================================================================================================
# 8.2- GLM stepwise regression - based on AIC (caret)
#====================================================================================================
# Inputs:
# - my.datf.train,my.datf.test = training and test dataset frames
# - sCV = type of cross-validation such as repeatedcv, LOOCV, etc.
# - iSplit = index of splitalse
# - fDet = flag for detais (True/F)
# - outFile = output file for GLM details
# Output:
# - list of statistics equal with the header introduced in the main script!
#   (tr = train, ts = test, both = tr+ts = full dataset)
# ---------------------------------------------------------------------------------------------------
GLMreg <- function(my.datf.train,my.datf.test,sCV,iSplit=1,fDet=F,outFile="") {
  # ------------------------------------------
  # only Linux or Mac: parallel calculations
  # ------------------------------------------
  #library(doMC)
  #registerDoMC(cores = 2) # CPU cores
  # ------------------------------------------
  library(caret)
  #attach(my.datf.train)    # make available the names of variables from training dataset
  net.c = my.datf.train[,1] # dependent variable is the first column in Training set
  RegrMethod <- "glmStepAIC" # type of regression
  
  # Define the CV conditions
  ctrl<- trainControl(method=sCV, number=10,repeats=10,
                      summaryFunction=defaultSummary)
  
  # Train the model using only training set
  set.seed(iSplit)
  glm.fit<- train(net.c~.,data=my.datf.train,
                  method='glmStepAIC', tuneLength=10, trControl=ctrl,
                  metric='RMSE')
  
  #------------------------------
  # Training RESULTS
  #------------------------------
  RMSE.tr  <- glm.fit$results[,2]
  R2.tr    <- glm.fit$results[,3]
  if (sCV == "repeatedcv"){ # if 10-fold CV
    RMSEsd.tr <- glm.fit$results[,4]
    R2sd.tr   <- glm.fit$results[,5]
  }
  if (sCV == "LOOCV"){ # if LOOCV
    RMSEsd.tr <- 0 # formulas will be added later!
    R2sd.tr   <- 0 # formulas will be added later!
  }
  
  #------------------------------------------------
  # RMSE & R^2, for train/test respectively
  #------------------------------------------------
  lm.train.res <- getTrainPerf(glm.fit)
  lm.test.res  <- postResample(predict(glm.fit,my.datf.test),my.datf.test[,1])
  
  #------------------------------------------------
  # Adj R2, Pearson correlation
  #------------------------------------------------
  pred.tr     <- predict(glm.fit,my.datf.train) # predicted Y
  pred.ts     <- predict(glm.fit,my.datf.test)  # predicted Y
  noFeats.fit <- length(predictors(glm.fit))    # no. of features from the fitted model
  Feats.fit   <- paste(predictors(glm.fit),collapse="+") # string with the features included in the fitted model
  
  ds.full     <- rbind(my.datf.train,my.datf.test)
  pred.both   <- predict(glm.fit,ds.full)       # predicted Y
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
  if (fDet==T) {   # if flag for details if T, print details about any resut
    write("RRegr package | eNanoMapper",                  file=outFile,append=T)
    write.table(paste("Regression method: ", RegrMethod), file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(paste("Split no.: ", iSplit), file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(paste("CV type: ", sCV),      file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table("Training Set Summary: ", file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(summary(my.datf.train),   file=outFile,append=T,sep=",",col.names=T,quote=F)
    write.table("Test Set Summary: ",  file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(summary(my.datf.test), file=outFile,append=T,sep=",",col.names=T,quote=F)   
    
    write.table("Fitting Summary: ", file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(data.frame(summary(glm.fit)$coefficients), file=outFile,append=T,sep=",",col.names=T,quote=F)
    
    write.table("Predictors: ",      file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(predictors(glm.fit), file=outFile,append=T,sep=",",col.names=T,quote=F)
    
    write.table("Trainig Results: ", file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(predictors(lm.train.res), file=outFile,append=T,sep=",",col.names=T,quote=F)
    write.table("Test Results: ", file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(predictors(lm.test.res), file=outFile,append=T,sep=",",col.names=T,quote=F)
    
    write.table("Full Statistics: ", file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(my.stats, file=outFile,append=T,sep=",",col.names=T,quote=F)
    
    # Variable Importance (max top 20)
    FeatImp <- varImp(glm.fit, scale = F)
    components = length(FeatImp)  # default plot all feature importance
    if (length(FeatImp)>20){     # if the number of features is greater than 20, use only 20
      components = 20
    }
    # Append feature importance to output details
    AppendList2CSv(FeatImp,outFile)
    
    fitModel <- glm.fit$finalModel
    
    # =============================================================================
    # Assessment of Applicability Domain (plot leverage)
    # =============================================================================
    
    # Residuals
    resids <- residuals(fitModel) # residuals
    write.table("Residuals of the fitted model: ",file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(data.frame(resids), file=outFile,append=T,sep=",",col.names=T,row.names=T, quote=F) # write residuals
    
    # Leverage / Hat values
    hat.fit <- hatvalues(fitModel)          # hat values
    hat.fit.df <- as.data.frame(hat.fit)    # hat data frame
    hat.mean <- mean(hat.fit)               # mean hat values
    hat.fit.df$warn <- ifelse(hat.fit.df[, 'hat.fit']>3*hat.mean, 'x3',ifelse(hat.fit.df[, 'hat.fit']>2*hat.mean, 'x2', '-' ))
    
    write.table("Leverage output: ",file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(paste("Mean of hat values: ", hat.mean), file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table("Leverage / Hat values with warnings (X3 & X2 = values 3 & 2 times than hat mean): ",file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(hat.fit.df, file=outFile,append=T,sep=",",col.names=T,row.names=T, quote=F) # write hat values and the levels X3, X2 (of hat mean)
    
    #THRESHOLD values: 3m/n, where m is the number of parameters, and n number of observations
    thresh.lever<- (3*(dim(my.datf.train)[2]-1))/dim(my.datf.train)[1] # leverage thresh
    hat.problems<- data.frame(hat.fit[hat.fit>thresh.lever]) # points with high leverage
    
    write.table(paste("Leverage Threshold: ", thresh.lever), file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table("Points with leverage > threshold: ",file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(hat.problems, file=outFile,append=T,sep=",",col.names=T,row.names=T, quote=F)
    
    # Cook's distance
    cook.dists<- cooks.distance(fitModel)
    cutoff.Cook <- 4/((nrow(my.datf.train)-length(fitModel$coefficients)-2)) # Cook's distance cutoff
    
    write.table("Cook's distances output: ",file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(paste("Cook's distance cutoff: ", cutoff.Cook), file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table("Cook's distances: ",file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(data.frame(cook.dists), file=outFile,append=T,sep=",",col.names=T,row.names=T, quote=F) # write residuals
    
    # Influence
    infl <- influence(fitModel)#produces several statistics of the kind
    
    write.table("Point influence output: ",file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table("Influences: ",file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(data.frame(infl), file=outFile,append=T,sep=",",col.names=T,row.names=T, quote=F) # write residuals
    
    # PDF with 12 plots
    # --------------------------------------------------------------
    pdf(file=paste(outFile,".",sCV,".","split",iSplit,".pdf",sep=""))
    # par(mfrow = c(3, 4)) # all plots into one page!
    
    plot(my.datf.train[,1],pred.tr,xlab="Yobs", ylab="Ypred", type="b", main="Train Yobs-Ypred") # plot 1
    plot(my.datf.test[,1], pred.ts,xlab="Yobs", ylab="Ypred", type="b", main="Test Yobs-Ypred")  # plot 2
    dotchart(as.matrix(FeatImp$importance),main="Feature Importance")                            # plot 3
    
    # Fitted vs Residuals - plot 4
    plot(fitted(fitModel),residuals(fitModel),
         main="Fitted vs. Residuals for Fitted Model",
         xlab="Fitted", ylab="Residuals")
    abline(h = 0, lty = 2)
    
    # Leverage plots - plot 5
    plot(hat.fit, type = "h",
         main="Leverage for Fitted Model",
         xlab="Index", ylab="Hat")
    abline(h = thresh.lever, lty = 2, col="red") # leverage thresh
    
    # Cook's distance - plot 6
    plot(cook.dists,
         main="Cook's Distance for Fitted Model",
         xlab="Index", ylab="Cook Distance")
    
    for (p in 1:6) {
      plot(fitModel, which=p, cook.levels=cutoff.Cook) # 6 standard fitting plots
    }
    
    # plot(FeatImp, top = components,main="Feature Importance") # ERROR !
    dev.off()
    # --------------------------------------------------------------
  }
  
  # my.stats.full <- c(my.stats.dsInfo,my.stats.10CV,my.stats.LOOCV)   # merge the CV results into one list that contains the names of each field!
  
  return(my.stats)  # return a list with statistics
}

#====================================================================================================
# 8.3. PLS regression (caret)
#====================================================================================================
PLSreg <- function(my.datf.train,my.datf.test,sCV,iSplit=1,fDet=F,outFile="") {
  # ------------------------------------------
  # only Linux or Mac: parallel calculations
  # ------------------------------------------
  #library(doMC)
  #registerDoMC(cores = 2) # CPU cores
  # ------------------------------------------
  library(caret)
  #attach(my.datf.train)   # make available the names of variables from training dataset
  net.c = my.datf.train[,1] # dependent variable is the first column in Training set
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
                   "CVtype"       = sCV,                    # from function param
                   "NoModelFeats" = as.numeric(noFeats.fit),
                   "ModelFeats"   = Feats.fit,
                   "adjR2.tr"  = as.numeric(adjR2.tr),
                   
                   "RMSE.tr"   = as.numeric(min(RMSE.tr)),  # these 4 lines correspond to the min of RMSE.tr !!!
                   "R2.tr"     = as.numeric(R2.tr[which.min(RMSE.tr)]),  
                   "RMSEsd.tr" = as.numeric(RMSEsd.tr[which.min(RMSE.tr)]),
                   "R2sd.tr"   = as.numeric(R2sd.tr[which.min(RMSE.tr)]),
                   
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
  if (fDet==T) {   # if flag for details if true, print details about any resut
    write("RRegr package | eNanoMapper", file=outFile,append=T)
    write.table(paste("Regression method: ", RegrMethod), file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(paste("Split no.: ", iSplit), file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(paste("CV type: ", sCV),  file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table("Training Set Summary: ", file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(summary(my.datf.train), file=outFile,append=T,sep=",",col.names=T,quote=F)
    write.table("Test Set Summary: ",   file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(summary(my.datf.test),  file=outFile,append=T,sep=",",col.names=T,quote=F)   
    
    
    write.table("Predictors: ", file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(predictors(pls.fit), file=outFile,append=T,sep=",",col.names=T,quote=F)
    
    write.table("Trainig Results: ",      file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(predictors(lm.train.res), file=outFile,append=T,sep=",",col.names=T,quote=F)
    write.table("Test Results: ",        file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(predictors(lm.test.res), file=outFile,append=T,sep=",",col.names=T,quote=F)
    
    write.table("Full Statistics: ", file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(my.stats, file=outFile,append=T,sep=",",col.names=T,quote=F)
    
    # Variable Importance (max top 20)
    FeatImp <- varImp(pls.fit, scale = F)
    components = length(FeatImp)  # default plot all feature importance
    if (length(FeatImp)>20){     # if the number of features is greater than 20, use only 20
      components = 20
    }
    # Append feature importance to output details
    AppendList2CSv(FeatImp,outFile)
    
    fitModel <- pls.fit$finalModel
    
    # =============================================================================
    # Assessment of Applicability Domain (plot leverage)
    # =============================================================================
    
    # Residuals
    resids <- residuals(fitModel) # residuals
    write.table("Residuals of the fitted model: ",file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(data.frame(resids), file=outFile,append=T,sep=",",col.names=T,row.names=T, quote=F) # write residuals
    
    # ADDED !
    predVals.pls.ad <- pred.ts
    Traind.pls= as.matrix(my.datf.train)
    Testd.pls = as.matrix(my.datf.test)
    Hat.train = diag(Traind.pls %*% solve(t(Traind.pls) %*%(Traind.pls), tol=1e-40)  %*% t(Traind.pls))
    Hat.test  = diag(Testd.pls  %*% solve(t(Traind.pls) %*%(Traind.pls), tol=1e-40)  %*% t(Testd.pls))  
    
    # Leverage / Hat values
    hat.fit <- Hat.test          # hat values
    hat.fit.df <- as.data.frame(hat.fit)    # hat data frame
    hat.mean <- mean(hat.fit)               # mean hat values
    hat.fit.df$warn <- ifelse(hat.fit.df[, 'hat.fit']>3*hat.mean, 'x3',ifelse(hat.fit.df[, 'hat.fit']>2*hat.mean, 'x2', '-' ))
    
    write.table("Leverage output: ",file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(paste("Mean of hat values: ", hat.mean), file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table("Leverage / Hat values with warnings (X3 & X2 = values 3 & 2 times than hat mean): ",file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(hat.fit.df, file=outFile,append=T,sep=",",col.names=T,row.names=T, quote=F) # write hat values and the levels X3, X2 (of hat mean)
    
    #THRESHOLD values: 3m/n, where m is the number of parameters, and n number of observations
    thresh.lever<- (3*(dim(my.datf.train)[2]-1))/dim(my.datf.train)[1] # leverage thresh
    hat.problems<- data.frame(hat.fit[hat.fit>thresh.lever]) # points with high leverage
    
    write.table(paste("Leverage Threshold: ", thresh.lever), file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table("Points with leverage > threshold: ",file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(hat.problems, file=outFile,append=T,sep=",",col.names=T,row.names=T, quote=F)
    
    # Cook's distance
    # Influence
    
    # PDF plots
    # --------------------------------------------------------------
    pdf(file=paste(outFile,".",sCV,".","split",iSplit,".pdf",sep=""))
    plot(my.datf.train[,1],pred.tr,xlab="Yobs", ylab="Ypred", type="b", main="Train Yobs-Ypred")
    plot(my.datf.test[,1], pred.ts,xlab="Yobs", ylab="Ypred", type="b", main="Test Yobs-Ypred")
    dotchart(as.matrix(FeatImp$importance),main="Feature Importance")
    
    # Fitted vs Residuals
    plot(fitted(fitModel),residuals(fitModel),
         main="Fitted vs. Residuals for Fitted Model",
         xlab="Fitted", ylab="Residuals")
    abline(h = 0, lty = 2)
    
    # Leverage plots
    plot(hat.fit, type = "h",
         main="Leverage for Fitted Model",
         xlab="Index", ylab="Hat")
    abline(h = thresh.lever, lty = 2, col="red") # leverage thresh
    
    dev.off()
    # --------------------------------------------------------------
  }
  
  return(my.stats)  # return a list with statistics
}

#====================================================================================================
# 8.4 Lasso Regression (caret)
#====================================================================================================
LASSOreg <- function(my.datf.train,my.datf.test,sCV,iSplit=1,fDet=F,outFile="") {
  # ------------------------------------------
  # only Linux or Mac: parallel calculations
  # ------------------------------------------
  #library(doMC)
  #registerDoMC(cores = 2) # CPU cores
  # ------------------------------------------
  library(caret)
  #attach(my.datf.train)   # make available the names of variables from training dataset
  net.c = my.datf.train[,1] # dependent variable is the first column in Training set
  RegrMethod <- "lasso.RMSE" # type of regression
  
  # Define the CV conditions
  ctrl<- trainControl(method = sCV, number = 10,repeats = 10,
                      summaryFunction = defaultSummary)
  
  # Train the model using only training set
  set.seed(iSplit)
  las.fit<- train(net.c~.,data=my.datf.train,
                  method='lasso', tuneLength = 10, trControl = ctrl,
                  metric='RMSE' ,tuneGrid=expand.grid(.fraction= seq(0.1,1,by=0.1)))
  
  #------------------------------
  # Training RESULTS
  #------------------------------
  RMSE.tr  <- las.fit$results[,2]
  R2.tr    <- las.fit$results[,3]
  if (sCV == "repeatedcv"){ # if 10-fold CV
    RMSEsd.tr <- las.fit$results[,4]
    R2sd.tr   <- las.fit$results[,5]
  }
  if (sCV == "LOOCV"){ # if LOOCV
    RMSEsd.tr <- 0 # formulas will be added later!
    R2sd.tr   <- 0 # formulas will be added later!
  }
  
  #------------------------------------------------
  # RMSE & R^2, for train/test respectively
  #------------------------------------------------
  lm.train.res <- getTrainPerf(las.fit)
  lm.test.res  <- postResample(predict(las.fit,my.datf.test),my.datf.test[,1])
  
  #------------------------------------------------
  # Adj R2, Pearson correlation
  #------------------------------------------------
  pred.tr     <- predict(las.fit,my.datf.train) # predicted Y
  pred.ts     <- predict(las.fit,my.datf.test)  # predicted Y
  noFeats.fit <- length(predictors(las.fit))    # no. of features from the fitted model
  Feats.fit   <- paste(predictors(las.fit),collapse="+") # string with the features included in the fitted model
  
  ds.full     <- rbind(my.datf.train,my.datf.test)
  pred.both   <- predict(las.fit,ds.full)       # predicted Y
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
                   
                   "RMSE.tr"   = as.numeric(min(RMSE.tr)),  # these 4 lines correspond to the min of RMSE.tr !!!
                   "R2.tr"     = as.numeric(R2.tr[which.min(RMSE.tr)]),  
                   "RMSEsd.tr" = as.numeric(RMSEsd.tr[which.min(RMSE.tr)]),
                   "R2sd.tr"   = as.numeric(R2sd.tr[which.min(RMSE.tr)]),
           
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
  if (fDet==T) {   # if flag for details if true, print details about any resut
    write("RRegr package | eNanoMapper", file=outFile,append=T)
    write.table(paste("Regression method: ", RegrMethod), file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(paste("Split no.: ", iSplit), file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(paste("CV type: ", sCV),      file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table("Training Set Summary: ", file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(summary(my.datf.train),   file=outFile,append=T,sep=",",col.names=T,quote=F)
    write.table("Test Set Summary: ",  file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(summary(my.datf.test), file=outFile,append=T,sep=",",col.names=T,quote=F)   
    
    
    write.table("Predictors: ",      file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(predictors(las.fit), file=outFile,append=T,sep=",",col.names=T,quote=F)
    
    write.table("Trainig Results: ",      file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    AppendList2CSv(predictors(lm.train.res),outFile)
    #write.table(predictors(lm.train.res), file=outFile,append=T,sep=",",col.names=T,quote=F)
    write.table("Test Results: ",        file=outFile,append=T,sep=",",col.names=F,quote=F)
    AppendList2CSv(predictors(lm.test.res),outFile)
    #write.table(predictors(lm.test.res), file=outFile,append=T,sep=",",col.names=T,quote=F)
    
    write.table("Full Statistics: ", file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(my.stats, file=outFile,append=T,sep=",",col.names=T,quote=F)
    
    # Variable Importance (max top 20)
    FeatImp <- varImp(las.fit, scale = F)
    components = length(FeatImp)  # default plot all feature importance
    if (length(FeatImp)>20){     # if the number of features is greater than 20, use only 20
      components = 20
    }
    # Append feature importance to output details
    AppendList2CSv(FeatImp,outFile)
    
    fitModel <- las.fit$finalModel
    
    # =============================================================================
    # Assessment of Applicability Domain (plot leverage)
    # =============================================================================
    
    # Residuals
    resids <- residuals(fitModel) # residuals
    write.table("Residuals of the fitted model: ",file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(data.frame(resids), file=outFile,append=T,sep=",",col.names=T,row.names=T, quote=F) # write residuals
    
    # ADDED !
    predVals.pls.ad <- pred.ts
    Traind.pls= as.matrix(my.datf.train)
    Testd.pls = as.matrix(my.datf.test)
    Hat.train = diag(Traind.pls %*% solve(t(Traind.pls) %*%(Traind.pls), tol=1e-40)  %*% t(Traind.pls))
    Hat.test  = diag(Testd.pls  %*% solve(t(Traind.pls) %*%(Traind.pls), tol=1e-40)  %*% t(Testd.pls))  
    
    # Leverage / Hat values
    hat.fit <- Hat.test          # hat values
    hat.fit.df <- as.data.frame(hat.fit)    # hat data frame
    hat.mean <- mean(hat.fit)               # mean hat values
    hat.fit.df$warn <- ifelse(hat.fit.df[, 'hat.fit']>3*hat.mean, 'x3',ifelse(hat.fit.df[, 'hat.fit']>2*hat.mean, 'x2', '-' ))
    
    write.table("Leverage output: ",file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(paste("Mean of hat values: ", hat.mean), file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table("Leverage / Hat values with warnings (X3 & X2 = values 3 & 2 times than hat mean): ",file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(hat.fit.df, file=outFile,append=T,sep=",",col.names=T,row.names=T, quote=F) # write hat values and the levels X3, X2 (of hat mean)
    
    #THRESHOLD values: 3m/n, where m is the number of parameters, and n number of observations
    thresh.lever<- (3*(dim(my.datf.train)[2]-1))/dim(my.datf.train)[1] # leverage thresh
    hat.problems<- data.frame(hat.fit[hat.fit>thresh.lever]) # points with high leverage
    
    write.table(paste("Leverage Threshold: ", thresh.lever), file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table("Points with leverage > threshold: ",file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(hat.problems, file=outFile,append=T,sep=",",col.names=T,row.names=T, quote=F)
    
    # Cook's distance
    # Influence
    
    # PDF plots
    # --------------------------------------------------------------
    pdf(file=paste(outFile,".",sCV,".","split",iSplit,".pdf",sep=""))
    plot(my.datf.train[,1],pred.tr,xlab="Yobs", ylab="Ypred", type="b", main="Train Yobs-Ypred")
    plot(my.datf.test[,1], pred.ts,xlab="Yobs", ylab="Ypred", type="b", main="Test Yobs-Ypred")
    dotchart(as.matrix(FeatImp$importance),main="Feature Importance")
    
    # Fitted vs Residuals
    plot(fitted(fitModel),residuals(fitModel),
         main="Fitted vs. Residuals for Fitted Model",
         xlab="Fitted", ylab="Residuals")
    abline(h = 0, lty = 2)
    
    # Leverage plots
    plot(hat.fit, type = "h",
         main="Leverage for Fitted Model",
         xlab="Index", ylab="Hat")
    abline(h = thresh.lever, lty = 2, col="red") # leverage thresh
    
    dev.off()
    # --------------------------------------------------------------
  }
  
  return(my.stats)  # return a list with statistics
}

#====================================================================================================
# 8.5. RBF network with the DDA algorithm regression (caret)
#====================================================================================================
RBF_DDAreg <- function(my.datf.train,my.datf.test,sCV,iSplit=1,fDet=F,outFile="") {
  # ------------------------------------------
  # only Linux or Mac: parallel calculations
  # ------------------------------------------
  #library(doMC)
  #registerDoMC(cores = 2) # CPU cores
  # ------------------------------------------
  library(caret)
  #attach(my.datf.train)   # make available the names of variables from training dataset
  net.c = my.datf.train[,1] # dependent variable is the first column in Training set
  RegrMethod <- "rbfDDA" # type of regression
  
  # Define the CV conditions
  ctrl<- trainControl(method=sCV,number=10,repeats=10,
                      summaryFunction=defaultSummary)
  
  # Train the model using only training set
  set.seed(iSplit)
  
  rbf.fit<- train(net.c~.,data=my.datf.train,
                  method='rbfDDA',trControl=ctrl,
                  tuneGrid=expand.grid(.negativeThreshold=seq(0,1,0.1)))
  
  #------------------------------
  # Training RESULTS
  #------------------------------
  RMSE.tr  <- rbf.fit$results[,2]
  R2.tr    <- rbf.fit$results[,3]
  if (sCV == "repeatedcv"){ # if 10-fold CV
    RMSEsd.tr <- rbf.fit$results[,4]
    R2sd.tr   <- rbf.fit$results[,5]
  }
  if (sCV == "LOOCV"){ # if LOOCV
    RMSEsd.tr <- 0 # formulas will be added later!
    R2sd.tr   <- 0 # formulas will be added later!
  }
  
  #------------------------------------------------
  # RMSE & R^2, for train/test respectively
  #------------------------------------------------
  lm.train.res <- getTrainPerf(rbf.fit)
  lm.test.res  <- postResample(predict(rbf.fit,my.datf.test),my.datf.test[,1])
  
  #------------------------------------------------
  # Adj R2, Pearson correlation
  #------------------------------------------------
  pred.tr     <- predict(rbf.fit,my.datf.train) # predicted Y
  pred.ts     <- predict(rbf.fit,my.datf.test)  # predicted Y
  noFeats.fit <- length(predictors(rbf.fit))    # no. of features from the fitted model
  Feats.fit   <- paste(predictors(rbf.fit),collapse="+") # string with the features included in the fitted model
  
  ds.full     <- rbind(my.datf.train,my.datf.test)
  pred.both   <- predict(rbf.fit,ds.full)       # predicted Y
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
                   
                   "RMSE.tr"   = as.numeric(min(RMSE.tr)),  # these 4 lines correspond to the min of RMSE.tr !!!
                   "R2.tr"     = as.numeric(R2.tr[which.min(RMSE.tr)]),  
                   "RMSEsd.tr" = as.numeric(RMSEsd.tr[which.min(RMSE.tr)]),
                   "R2sd.tr"   = as.numeric(R2sd.tr[which.min(RMSE.tr)]),
                   
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
  if (fDet==T) {   # if flag for details if true, print details about any resut
    write("RRegr package | eNanoMapper", file=outFile,append=T)
    write.table(paste("Regression method: ", RegrMethod), file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(paste("Split no.: ", iSplit),             file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(paste("CV type: ", sCV),  file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table("Training Set Summary: ", file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(summary(my.datf.train), file=outFile,append=T,sep=",",col.names=T,quote=F)
    write.table("Test Set Summary: ",   file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(summary(my.datf.test),  file=outFile,append=T,sep=",",col.names=T,quote=F)   
    
    
    write.table("Predictors: ",      file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(predictors(rbf.fit), file=outFile,append=T,sep=",",col.names=T,quote=F)
    
    write.table("Trainig Results: ",      file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(predictors(lm.train.res), file=outFile,append=T,sep=",",col.names=T,quote=F)
    write.table("Test Results: ",        file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(predictors(lm.test.res), file=outFile,append=T,sep=",",col.names=T,quote=F)
    
    write.table("Full Statistics: ", file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(my.stats, file=outFile,append=T,sep=",",col.names=T,row.names=F,quote=F)
    
    # Variable Importance (max top 20)
    FeatImp <- varImp(rbf.fit, scale = F)
    components = length(FeatImp)  # default plot all feature importance
    if (length(FeatImp)>20){     # if the number of features is greater than 20, use only 20
      components = 20
    }
    # Append feature importance to output details
    AppendList2CSv(FeatImp,outFile)
    
    fitModel <- rbf.fit$finalModel
    
    # =============================================================================
    # Assessment of Applicability Domain (plot leverage)
    # =============================================================================
    
    # Residuals
    resids <- residuals(fitModel) # residuals
    write.table("Residuals of the fitted model: ",file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(data.frame(resids), file=outFile,append=T,sep=",",col.names=T,row.names=T, quote=F) # write residuals
    
    # ADDED !
    predVals.pls.ad <- pred.ts
    Traind.pls= as.matrix(my.datf.train)
    Testd.pls = as.matrix(my.datf.test)
    Hat.train = diag(Traind.pls %*% solve(t(Traind.pls) %*%(Traind.pls), tol=1e-40)  %*% t(Traind.pls))
    Hat.test  = diag(Testd.pls  %*% solve(t(Traind.pls) %*%(Traind.pls), tol=1e-40)  %*% t(Testd.pls))  
    
    # Leverage / Hat values
    hat.fit <- Hat.test          # hat values
    hat.fit.df <- as.data.frame(hat.fit)    # hat data frame
    hat.mean <- mean(hat.fit)               # mean hat values
    hat.fit.df$warn <- ifelse(hat.fit.df[, 'hat.fit']>3*hat.mean, 'x3',ifelse(hat.fit.df[, 'hat.fit']>2*hat.mean, 'x2', '-' ))
    
    write.table("Leverage output: ",file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(paste("Mean of hat values: ", hat.mean), file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table("Leverage / Hat values with warnings (X3 & X2 = values 3 & 2 times than hat mean): ",file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(hat.fit.df, file=outFile,append=T,sep=",",col.names=T,row.names=T, quote=F) # write hat values and the levels X3, X2 (of hat mean)
    
    #THRESHOLD values: 3m/n, where m is the number of parameters, and n number of observations
    thresh.lever<- (3*(dim(my.datf.train)[2]-1))/dim(my.datf.train)[1] # leverage thresh
    hat.problems<- data.frame(hat.fit[hat.fit>thresh.lever]) # points with high leverage
    
    write.table(paste("Leverage Threshold: ", thresh.lever), file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table("Points with leverage > threshold: ",file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(hat.problems, file=outFile,append=T,sep=",",col.names=T,row.names=T, quote=F)
    
    # Cook's distance
    # Influence
    
    # PDF plots
    # --------------------------------------------------------------
    pdf(file=paste(outFile,".",sCV,".","split",iSplit,".pdf",sep=""))
    plot(my.datf.train[,1],pred.tr,xlab="Yobs", ylab="Ypred", type="b", main="Train Yobs-Ypred")
    plot(my.datf.test[,1], pred.ts,xlab="Yobs", ylab="Ypred", type="b", main="Test Yobs-Ypred")
    dotchart(as.matrix(FeatImp$importance),main="Feature Importance")
    
    # Fitted vs Residuals
    plot(fitted(fitModel),residuals(fitModel),
         main="Fitted vs. Residuals for Fitted Model",
         xlab="Fitted", ylab="Residuals")
    abline(h = 0, lty = 2)
    
    # Leverage plots
    plot(hat.fit, type = "h",
         main="Leverage for Fitted Model",
         xlab="Index", ylab="Hat")
    abline(h = thresh.lever, lty = 2, col="red") # leverage thresh
    
    dev.off()
    # --------------------------------------------------------------
  }
  
  return(my.stats)  # return a list with statistics
}

#====================================================================================================
# 8.6 SVM Radial Regression (caret)
#====================================================================================================
SVLMreg <- function(my.datf.train,my.datf.test,sCV,iSplit=1,fDet=F,outFile="") {
  # ------------------------------------------
  # only Linux or Mac: parallel calculations
  # ------------------------------------------
  #library(doMC)
  #registerDoMC(cores = 2) # CPU cores
  # ------------------------------------------
  library(caret)
  #attach(my.datf.train)   # make available the names of variables from training dataset
  net.c = my.datf.train[,1] # dependent variable is the first column in Training set
  RegrMethod <- "svmRadial.RMSE" # type of regression
  
  # Define the CV conditions
  ctrl<- trainControl(method=sCV,number=10,repeats=10,
                      summaryFunction=defaultSummary)

  # Train the model using only training set
  set.seed(iSplit)
  
  svmL.fit<- train(net.c~.,data=my.datf.train,
                   method='svmRadial',tuneLength=10,trControl=ctrl,
                   metric='RMSE',
                   tuneGrid=expand.grid(.sigma=seq(0,1,0.1),.C= c(1:10)))
  
  #------------------------------
  # Training RESULTS
  #------------------------------
  RMSE.tr  <- svmL.fit$results[,2]
  R2.tr    <- svmL.fit$results[,3]
  if (sCV == "repeatedcv"){ # if 10-fold CV
    RMSEsd.tr <- svmL.fit$results[,4]
    R2sd.tr   <- svmL.fit$results[,5]
  }
  if (sCV == "LOOCV"){ # if LOOCV
    RMSEsd.tr <- 0 # formulas will be added later!
    R2sd.tr   <- 0 # formulas will be added later!
  }
  
  #------------------------------------------------
  # RMSE & R^2, for train/test respectively
  #------------------------------------------------
  lm.train.res <- getTrainPerf(svmL.fit)
  lm.test.res  <- postResample(predict(svmL.fit,my.datf.test),my.datf.test[,1])
  
  #------------------------------------------------
  # Adj R2, Pearson correlation
  #------------------------------------------------
  pred.tr     <- predict(svmL.fit,my.datf.train) # predicted Y
  pred.ts     <- predict(svmL.fit,my.datf.test)  # predicted Y
  noFeats.fit <- length(predictors(svmL.fit))    # no. of features from the fitted model
  Feats.fit   <- paste(predictors(svmL.fit),collapse="+") # string with the features included in the fitted model
  
  ds.full     <- rbind(my.datf.train,my.datf.test)
  pred.both   <- predict(svmL.fit,ds.full)       # predicted Y
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
                   
                   "RMSE.tr"   = as.numeric(min(RMSE.tr)),  # these 4 lines correspond to the min of RMSE.tr !!!
                   "R2.tr"     = as.numeric(R2.tr[which.min(RMSE.tr)]),  
                   "RMSEsd.tr" = as.numeric(RMSEsd.tr[which.min(RMSE.tr)]),
                   "R2sd.tr"   = as.numeric(R2sd.tr[which.min(RMSE.tr)]),
                   
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
  if (fDet==T) {   # if flag for details if true, print details about any resut
    write("RRegr package | eNanoMapper", file=outFile,append=T)
    write.table(paste("Regression method: ", RegrMethod), file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(paste("Split no.: ", iSplit), file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(paste("CV type: ", sCV),      file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table("Training Set Summary: ", file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(summary(my.datf.train),   file=outFile,append=T,sep=",",col.names=T,quote=F)
    write.table("Test Set Summary: ",  file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(summary(my.datf.test), file=outFile,append=T,sep=",",col.names=T,quote=F)   
    
    
    write.table("Predictors: ",       file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(predictors(svmL.fit), file=outFile,append=T,sep=",",col.names=T,quote=F)
    
    write.table("Trainig Results: ",      file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(predictors(lm.train.res), file=outFile,append=T,sep=",",col.names=T,quote=F)
    write.table("Test Results: ",        file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(predictors(lm.test.res), file=outFile,append=T,sep=",",col.names=T,quote=F)
    
    write.table("Full Statistics: ", file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(my.stats, file=outFile,append=T,sep=",",col.names=T,quote=F)
    
    # Variable Importance (max top 20)
    FeatImp <- varImp(svmL.fit, scale = F)
    components = length(FeatImp)  # default plot all feature importance
    if (length(FeatImp)>20){     # if the number of features is greater than 20, use only 20
      components = 20
    }
    # Append feature importance to output details
    AppendList2CSv(FeatImp,outFile)
    
    fitModel <- svmL.fit$finalModel
    
    # =============================================================================
    # Assessment of Applicability Domain (plot leverage)
    # =============================================================================
    
    # Residuals
    resids <- residuals(fitModel) # residuals
    write.table("Residuals of the fitted model: ",file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(data.frame(resids), file=outFile,append=T,sep=",",col.names=T,row.names=T, quote=F) # write residuals
    
    # ADDED !
    predVals.pls.ad <- pred.ts
    Traind.pls= as.matrix(my.datf.train)
    Testd.pls = as.matrix(my.datf.test)
    Hat.train = diag(Traind.pls %*% solve(t(Traind.pls) %*%(Traind.pls), tol=1e-40)  %*% t(Traind.pls))
    Hat.test  = diag(Testd.pls  %*% solve(t(Traind.pls) %*%(Traind.pls), tol=1e-40)  %*% t(Testd.pls))  
    
    # Leverage / Hat values
    hat.fit <- Hat.test          # hat values
    hat.fit.df <- as.data.frame(hat.fit)    # hat data frame
    hat.mean <- mean(hat.fit)               # mean hat values
    hat.fit.df$warn <- ifelse(hat.fit.df[, 'hat.fit']>3*hat.mean, 'x3',ifelse(hat.fit.df[, 'hat.fit']>2*hat.mean, 'x2', '-' ))
    
    write.table("Leverage output: ",file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(paste("Mean of hat values: ", hat.mean), file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table("Leverage / Hat values with warnings (X3 & X2 = values 3 & 2 times than hat mean): ",file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(hat.fit.df, file=outFile,append=T,sep=",",col.names=T,row.names=T, quote=F) # write hat values and the levels X3, X2 (of hat mean)
    
    #THRESHOLD values: 3m/n, where m is the number of parameters, and n number of observations
    thresh.lever<- (3*(dim(my.datf.train)[2]-1))/dim(my.datf.train)[1] # leverage thresh
    hat.problems<- data.frame(hat.fit[hat.fit>thresh.lever]) # points with high leverage
    
    write.table(paste("Leverage Threshold: ", thresh.lever), file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table("Points with leverage > threshold: ",file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(hat.problems, file=outFile,append=T,sep=",",col.names=T,row.names=T, quote=F)
    
    # Cook's distance
    # Influence
    
    # PDF plots
    # --------------------------------------------------------------
    pdf(file=paste(outFile,".",sCV,".","split",iSplit,".pdf",sep=""))
    plot(my.datf.train[,1],pred.tr,xlab="Yobs", ylab="Ypred", type="b", main="Train Yobs-Ypred")
    plot(my.datf.test[,1], pred.ts,xlab="Yobs", ylab="Ypred", type="b", main="Test Yobs-Ypred")
    dotchart(as.matrix(FeatImp$importance),main="Feature Importance")
    
    # Fitted vs Residuals
    plot(fitted(fitModel),residuals(fitModel),
         main="Fitted vs. Residuals for Fitted Model",
         xlab="Fitted", ylab="Residuals")
    abline(h = 0, lty = 2)
    
    # Leverage plots
    plot(hat.fit, type = "h",
         main="Leverage for Fitted Model",
         xlab="Index", ylab="Hat")
    abline(h = thresh.lever, lty = 2, col="red") # leverage thresh
    
    dev.off()
    # --------------------------------------------------------------
  }
  
  return(my.stats)  # return a list with statistics
}

#====================================================================================================
# 8.8 Neural Network Regression (caret)
#====================================================================================================
NNreg <- function(my.datf.train,my.datf.test,sCV,iSplit=1,fDet=F,outFile="") {
  # ------------------------------------------
  # only Linux or Mac: parallel calculations
  # ------------------------------------------
  #library(doMC)
  #registerDoMC(cores = 2) # CPU cores
  # ------------------------------------------
  library(caret)
  #attach(my.datf.train)   # make available the names of variables from training dataset
  net.c = my.datf.train[,1] # dependent variable is the first column in Training set
  RegrMethod <- "nnet" # type of regression
  
  # Define the CV conditions
  ctrl<- trainControl(method = sCV, number = 10,repeats = 10,
                      summaryFunction = defaultSummary)
  
  # Train the model using only training set
  set.seed(iSplit)
  
  nn.fit<- train(net.c~.,data=my.datf.train,
                 method = 'nnet',trControl = ctrl,
                 linout=T, trace = F,MaxNWts=20000,
                 #Grid of tuning parameters to try:
                 tuneGrid=expand.grid(.size=c(1,5,10,15),.decay=c(0,0.001,0.1)))
  #Grid parameters are appearing at the print out of the model
  #size==#of units in hidden layer, decay==parameter of weight decay (default:0)
  
  #------------------------------
  # Training RESULTS
  #------------------------------
  RMSE.tr  <- nn.fit$results[,2]
  R2.tr    <- nn.fit$results[,3]
  if (sCV == "repeatedcv"){ # if 10-fold CV
    RMSEsd.tr <- nn.fit$results[,4]
    R2sd.tr   <- nn.fit$results[,5]
  }
  if (sCV == "LOOCV"){ # if LOOCV
    RMSEsd.tr <- 0 # formulas will be added later!
    R2sd.tr   <- 0 # formulas will be added later!
  }
  
  #------------------------------------------------
  # RMSE & R^2, for train/test respectively
  #------------------------------------------------
  lm.train.res <- getTrainPerf(nn.fit)
  lm.test.res  <- postResample(predict(nn.fit,my.datf.test),my.datf.test[,1])
  
  #------------------------------------------------
  # Adj R2, Pearson correlation
  #------------------------------------------------
  pred.tr     <- predict(nn.fit,my.datf.train) # predicted Y
  pred.ts     <- predict(nn.fit,my.datf.test)  # predicted Y
  noFeats.fit <- length(predictors(nn.fit))    # no. of features from the fitted model
  Feats.fit   <- paste(predictors(nn.fit),collapse="+") # string with the features included in the fitted model
  
  ds.full     <- rbind(my.datf.train,my.datf.test)
  pred.both   <- predict(nn.fit,ds.full)       # predicted Y
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
                   
                   "RMSE.tr"   = as.numeric(min(RMSE.tr)),  # these 4 lines correspond to the min of RMSE.tr !!!
                   "R2.tr"     = as.numeric(R2.tr[which.min(RMSE.tr)]),  
                   "RMSEsd.tr" = as.numeric(RMSEsd.tr[which.min(RMSE.tr)]),
                   "R2sd.tr"   = as.numeric(R2sd.tr[which.min(RMSE.tr)]),
                   
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
  if (fDet==T) {   # if flag for details if true, print details about any resut
    write("RRegr package | eNanoMapper",                  file=outFile,append=T)
    write.table(paste("Regression method: ", RegrMethod), file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(paste("Split no.: ", iSplit), file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(paste("CV type: ", sCV),      file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table("Training Set Summary: ", file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(summary(my.datf.train),   file=outFile,append=T,sep=",",col.names=T,quote=F)
    write.table("Test Set Summary: ",  file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(summary(my.datf.test), file=outFile,append=T,sep=",",col.names=T,quote=F)   
    
    
    write.table("Predictors: ",     file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(predictors(nn.fit), file=outFile,append=T,sep=",",col.names=T,quote=F)
    
    write.table("Trainig Results: ",      file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(predictors(lm.train.res), file=outFile,append=T,sep=",",col.names=T,quote=F)
    write.table("Test Results: ",        file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(predictors(lm.test.res), file=outFile,append=T,sep=",",col.names=T,quote=F)
    
    write.table("Full Statistics: ", file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(my.stats,            file=outFile,append=T,sep=",",col.names=T,quote=F)
    
    # Variable Importance (max top 20)
    FeatImp <- varImp(nn.fit, scale = F)
    components = length(FeatImp)  # default plot all feature importance
    if (length(FeatImp)>20){     # if the number of features is greater than 20, use only 20
      components = 20
    }
    # Append feature importance to output details
    AppendList2CSv(FeatImp,outFile)
    
    fitModel <- nn.fit$finalModel
    
    # =============================================================================
    # Assessment of Applicability Domain (plot leverage)
    # =============================================================================
    
    # Residuals
    resids <- residuals(fitModel) # residuals
    write.table("Residuals of the fitted model: ",file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(data.frame(resids), file=outFile,append=T,sep=",",col.names=T,row.names=T, quote=F) # write residuals
    
    # ADDED !
    predVals.pls.ad <- pred.ts
    Traind.pls= as.matrix(my.datf.train)
    Testd.pls = as.matrix(my.datf.test)
    Hat.train = diag(Traind.pls %*% solve(t(Traind.pls) %*%(Traind.pls), tol=1e-40)  %*% t(Traind.pls))
    Hat.test  = diag(Testd.pls  %*% solve(t(Traind.pls) %*%(Traind.pls), tol=1e-40)  %*% t(Testd.pls))  
    
    # Leverage / Hat values
    hat.fit <- Hat.test          # hat values
    hat.fit.df <- as.data.frame(hat.fit)    # hat data frame
    hat.mean <- mean(hat.fit)               # mean hat values
    hat.fit.df$warn <- ifelse(hat.fit.df[, 'hat.fit']>3*hat.mean, 'x3',ifelse(hat.fit.df[, 'hat.fit']>2*hat.mean, 'x2', '-' ))
    
    write.table("Leverage output: ",file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(paste("Mean of hat values: ", hat.mean), file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table("Leverage / Hat values with warnings (X3 & X2 = values 3 & 2 times than hat mean): ",file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(hat.fit.df, file=outFile,append=T,sep=",",col.names=T,row.names=T, quote=F) # write hat values and the levels X3, X2 (of hat mean)
    
    #THRESHOLD values: 3m/n, where m is the number of parameters, and n number of observations
    thresh.lever<- (3*(dim(my.datf.train)[2]-1))/dim(my.datf.train)[1] # leverage thresh
    hat.problems<- data.frame(hat.fit[hat.fit>thresh.lever]) # points with high leverage
    
    write.table(paste("Leverage Threshold: ", thresh.lever), file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table("Points with leverage > threshold: ",file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(hat.problems, file=outFile,append=T,sep=",",col.names=T,row.names=T, quote=F)
    
    # Cook's distance
    # Influence
    
    # PDF plots
    # --------------------------------------------------------------
    pdf(file=paste(outFile,".",sCV,".","split",iSplit,".pdf",sep=""))
    plot(my.datf.train[,1],pred.tr,xlab="Yobs", ylab="Ypred", type="b", main="Train Yobs-Ypred")
    plot(my.datf.test[,1], pred.ts,xlab="Yobs", ylab="Ypred", type="b", main="Test Yobs-Ypred")
    dotchart(as.matrix(FeatImp$importance),main="Feature Importance")
    
    # Fitted vs Residuals
    plot(fitted(fitModel),residuals(fitModel),
         main="Fitted vs. Residuals for Fitted Model",
         xlab="Fitted", ylab="Residuals")
    abline(h = 0, lty = 2)
    
    # Leverage plots
    plot(hat.fit, type = "h",
         main="Leverage for Fitted Model",
         xlab="Index", ylab="Hat")
    abline(h = thresh.lever, lty = 2, col="red") # leverage thresh
    
    dev.off()
    # --------------------------------------------------------------
  }
  return(my.stats)  # return a list with statistics
}


# **************************************
# WRAPPER METHODS
# **************************************

GLMregW <- function(my.datf,my.datf.train,my.datf.test,fDet=F,outFile="") { # with wrapper
  # to be implemented
  return("")  # return statistics
}

#====================================================================================================
# 8.3W. PLS regression with wrapper feature selection (caret)
#====================================================================================================
PLSregWSel <- function(my.datf.train,my.datf.test,sCV,iSplit=1,fDet=F,outFile="") {
  # ------------------------------------------
  # only Linux or Mac: parallel calculations
  # ------------------------------------------
  #library(doMC)
  #registerDoMC(cores = 2) # CPU cores
  # ------------------------------------------
  library(caret)
  
  net.c = my.datf.train[,1] # dependent variable is the first column in Training set
  RegrMethod <- "pls.WSel" # type of regression
  
  # Define the CV conditions
  ctrlw <- rfeControl(method = 'boot', number = 25,saveDetails=T)
  ctrl  <- trainControl(method = sCV, number = 10,repeats = 1,
                        summaryFunction = defaultSummary,savePred=T)
  
  subsetsx<- seq(2,dim(my.datf.train)[2]-1, by = 10)
  
  # Train the model using only training set
  set.seed(iSplit)
  
  pls.fit<- rfe(net.c~.,data=my.datf.train, 
                method = 'pls',
                rfeControl = ctrlw, trControl=ctrl, sizes=subsetsx, importance=T,
                metric = 'RMSE',
                tuneGrid=expand.grid(.ncomp=c(1:5)))
  
  #------------------------------
  # Training RESULTS
  #------------------------------
  
  pls.fit.best <- subset(pls.fit$results, pls.fit$results$Variables == 5) # best selected fit
  
  RMSE.tr  <- pls.fit.best$RMSE
  R2.tr    <- pls.fit.best$Rsquared
  if (sCV == "repeatedcv"){ # if 10-fold CV
    RMSEsd.tr <- pls.fit.best$RMSESD
    R2sd.tr   <- pls.fit.best$RsquaredSD
  }
  if (sCV == "LOOCV"){ # if LOOCV
    RMSEsd.tr <- 0 # formulas will be added later!
    R2sd.tr   <- 0 # formulas will be added later!
  }
  
  #------------------------------------------------
  # RMSE & R^2, for train/test respectively
  #------------------------------------------------
  lm.train.res <- pls.fit # ??
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
                   "CVtype"       = sCV,                    # from function param
                   "NoModelFeats" = as.numeric(noFeats.fit),
                   "ModelFeats"   = Feats.fit,
                   "adjR2.tr"     = as.numeric(adjR2.tr),
                   
                   "RMSE.tr"   = as.numeric(RMSE.tr),  # these 4 lines correspond to the min of RMSE.tr !!!
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
  if (fDet==T) {   # if flag for details if True, print details about any resut
    write("RRegr package | eNanoMapper",                  file=outFile,append=T)
    write.table(paste("Regression method: ", RegrMethod), file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(paste("Split no.: ", iSplit), file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(paste("CV type: ", sCV),      file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table("Training Set Summary: ", file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(summary(my.datf.train),   file=outFile,append=T,sep=",",col.names=T,quote=F)
    write.table("Test Set Summary: ", file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(summary(my.datf.test),file=outFile,append=T,sep=",",col.names=T,quote=F)   
    
    
    write.table("Predictors: ",      file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(predictors(pls.fit), file=outFile,append=T,sep=",",col.names=T,quote=F)
    
    write.table("Trainig Results: ",      file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(predictors(lm.train.res), file=outFile,append=T,sep=",",col.names=T,quote=F)
    write.table("Test Results: ",        file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(predictors(lm.test.res), file=outFile,append=T,sep=",",col.names=T,quote=F)
    
    write.table("NNet variable importance: ", file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    AppendList2txt(varImp(pls.fit),outFile)
    
    write.table("Full Statistics: ", file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(my.stats, file=outFile,append=T,sep=",",col.names=T,quote=F)
    
    # Variable Importance (max top 20)
    FeatImp <- varImp(pls.fit, scale = F)
    components = length(FeatImp)  # default plot all feature importance
    if (length(FeatImp)>20){     # if the number of features is greater than 20, use only 20
      components = 20
    }
    # Append feature importance to output details
    AppendList2CSv(FeatImp,outFile)
    
    # PDF plots
    pdf(file=paste(outFile,".",sCV,".","split",iSplit,".pdf",sep=""))
    par(mfrow = c(2, 2))
    plot(my.datf.train[,1],pred.tr,xlab="Yobs", ylab="Ypred", type="b", main="Train Yobs-Ypred")
    plot(my.datf.test[,1], pred.ts,xlab="Yobs", ylab="Ypred", type="b", main="Test Yobs-Ypred")
    dotchart(as.matrix(FeatImp$importance),main="Feature Importance")
    #plot(FeatImp, top = components,main="Feature Importance")
    dev.off()
  }
  
  return(my.stats)  # return a list with statistics
}

#=======================================================================
# Y-randomization for the best model (Step 12)
#=======================================================================
#    - 1 splitting, 1 CV type, best method
#    - best.R2.ts will be compared with Yrand.R2.ts
#    - returns ratios DiffsR2/bestR2
# --------------------------------------------------
Yrandom<- function(best.reg,best.R2.ts,noYrand,outFile){
  cat("-> Best model Y-Randomization ...\n")
  ds[,1] <- sample(ds[,1]) # randomize Y values for the entire dataset
  
  # splitting dataset in training and test
  #---------------------------------------
  Yrand.R2.ts <- NULL     # all values of R2 for each Y randomization
  for (i in 1:noYrand){
    iSeed=i               
    dsList  <- DsSplit(ds,trainFrac,fDet,PathDataSet,iSeed) # return a list with 2 datasets = dsList$train, dsList$test
    # get train and test from the resulted list
    ds.train<- dsList$train
    ds.test <- dsList$test
    
    # Run the caret function with the method from the best method
    #    for one training-test split only; no details, we need only R2 values
    if (best.reg=="lm") {
      my.stats.reg  <- LMreg(ds.train,ds.test,"repeatedcv",i,F,ResBestF) # run GLM for each CV and regr method
    }
    if (best.reg=="glmStepAIC") {
      my.stats.reg  <- GLMreg(ds.train,ds.test,"repeatedcv",i,F,ResBestF) # run GLM for each CV and regr method
    }
    if (best.reg=="pls") {
      my.stats.reg  <- PLSreg(ds.train,ds.test,"repeatedcv",i,F,ResBestF) # run SVLM Radial for each CV and regr method
    }
    if (best.reg=="lasso") {
      my.stats.reg  <- LASSOreg(ds.train,ds.test,"repeatedcv",i,F,ResBestF) # run SVLM Radial for each CV and regr method
    }
    if (best.reg=="rbfDDA") {  
      my.stats.reg  <- RBF_DDAreg(ds.train,ds.test,"repeatedcv",i,F,ResBestF) # run SVLM Radial for each CV and regr method
    }
    if (best.reg=="svmRadial") {  
      my.stats.reg  <- SVLMreg(ds.train,ds.test,"repeatedcv",i,F,ResBestF) # run SVLM Radial for each CV and regr method
    }
    if (best.reg=="nnet") {  
      my.stats.reg  <- NNreg(ds.train,ds.test,"repeatedcv",i,F,ResBestF) # run NNet for each CV and regr method
    } 
    
    Yrand.R2.ts <- c(Yrand.R2.ts,my.stats.reg$R2.ts) # adding test R2 value Y randomization 
  }
  
  # get histogram for differences between best R2 and the values for each Y randomization
  R2diffs     <- abs(Yrand.R2.ts - best.R2.ts) # absolute differences between R2 values (best model vs Y randomized results)
  R2diffsPerBestR2 <- R2diffs/best.R2.ts        # the same difference in percents
  
  pdf(file=paste(outFile,".Yrand.Hist.pdf",sep=""))    # save histogram if ratio diffs R2 into PDF for Y random
  Yrand.hist  <- hist(R2diffsPerBestR2)         # draw histogram of the ratio diffs/Best R2 for Y random
  dev.off()
  
  write.table("Y randomization test: ",file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
  write.table("=====================", file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
  write.table("Diffs R2 (Best Model - Y rand):",file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
  AppendList2CSv(R2diffs, outFile)
  write.table("Summary Difs:",file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
  AppendList2CSv(summary(R2diffs), outFile)
  write.table("Ratio Diffs R2 / Best R2 (Best Model - Y rand):",file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
  AppendList2CSv(R2diffsPerBestR2, outFile)
  write.table("Summary Difs %:",file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
  AppendList2CSv(summary(R2diffsPerBestR2),outFile)
  
  return(R2diffsPerBestR2) # return the ratio of diffs with the best R2 from the same 
}
#====================================================================================================
#  Basic RandomForest
# Jose A. Seoane seoane@stanford.edu
# Carlos Fernandez-Lozano carlos.fernandez@udc.es
#====================================================================================================
RFreg <- function(my.datf.train,my.datf.test,sCV,iSplit=1,fDet=F,outFile="") {
  # ------------------------------------------
  # only Linux or Mac: parallel calculations
  # ------------------------------------------
  #library(doMC)
  #registerDoMC(cores = 2) # CPU cores
  # ------------------------------------------
  
  
  net.c = my.datf.train[,1]   # make available the names of variables from training dataset
  RegrMethod <- "rf" # type of regression
  
  # Define the CV conditions
  ctrl<- trainControl(method=sCV, number=5,repeats=2,
                      summaryFunction=defaultSummary)
  
  # Train the model using only training set
  set.seed(iSplit)
  rf.fit<- train(net.c~.,data=my.datf.train,
                 method='rf', tuneLength = 10,trControl=ctrl,
                 metric='RMSE')
  
  #------------------------------
  # Training RESULTS
  #------------------------------
  RMSE.tr  <- rf.fit$results[rownames(rf.fit$bestTune),2]
  R2.tr    <- rf.fit$results[rownames(rf.fit$bestTune),3]
  if (sCV == "repeatedcv"){ # if 10-fold CV
    RMSEsd.tr <- rf.fit$results[rownames(rf.fit$bestTune),4]
    R2sd.tr   <- rf.fit$results[rownames(rf.fit$bestTune),5]
  }
  if (sCV == "LOOCV"){ # if LOOCV
    RMSEsd.tr <- 0 # formulas will be added later!  TODOOOOOOOOOOOOOO
    R2sd.tr   <- 0 # formulas will be added later!
  }
  
  #------------------------------------------------
  # RMSE & R^2, for train/test respectively
  #------------------------------------------------
  rf.train.res <- getTrainPerf(rf.fit)
  rf.test.res  <- postResample(predict(rf.fit,my.datf.test),my.datf.test[,1])
  
  #------------------------------------------------
  # Adj R2, Pearson correlation
  #------------------------------------------------
  pred.tr     <- predict(rf.fit,my.datf.train) # predicted Y for training
  pred.ts     <- predict(rf.fit,my.datf.test)  # predicted Y for test
  noFeats.fit <- length(predictors(rf.fit))    # no. of features from the fitted model
  Feats.fit   <- paste(predictors(rf.fit),collapse="+") # string with the features included in the fitted model
  
  ds.full     <- rbind(my.datf.train,my.datf.test)
  pred.both   <- predict(rf.fit,ds.full)       # predicted Y
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
                   "RMSE.ts" = as.numeric((rf.test.res["RMSE"])),
                   "R2.ts"   = as.numeric((rf.test.res["Rsquared"])),
                   "corP.ts" = as.numeric(corP.ts),
                   "adjR2.both" = as.numeric(adjR2.both),
                   "RMSE.both"  = as.numeric(RMSE.both),
                   "R2.both"    = as.numeric(r2.both))
  
  #---------------------------------------------------------------------
  # Write to file DETAILS for GLM for each cross-validation method
  #---------------------------------------------------------------------
  if (fDet==T) {   # if flag for details if T, print details about any resut
    write("RRegr package | eNanoMapper", file=outFile, append=T)
    write.table(paste("Regression method: ", RegrMethod), file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(paste("Split no.: ", iSplit), file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(paste("CV type: ", sCV),      file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table("Training Set Summary: ",     file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(summary(my.datf.train),       file=outFile,append=T,sep=",",col.names=T,quote=F)
    write.table("Test Set Summary: ",         file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(summary(my.datf.test),        file=outFile,append=T,sep=",",col.names=T,quote=F)   
    
    
    write.table("Predictors: ",     file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(predictors(rf.fit), file=outFile,append=T,sep=",",col.names=T,quote=F)
    
    write.table("Trainig Results: ",     file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(predictors(rf.train.res),file=outFile,append=T,sep=",",col.names=T,quote=F)
    write.table("Test Results: ",        file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(predictors(rf.test.res), file=outFile,append=T,sep=",",col.names=T,quote=F)
    
    write.table("Full Statistics: ", file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(my.stats,            file=outFile,append=T,sep=",",col.names=T,quote=F)
    
    # Variable Importance (max top 20)
    FeatImp <- importance(rf.fit$finalModel, scale = T)
    FeatImp = FeatImp[order(FeatImp,decreasing=T),]
    components = length(FeatImp)  # default plot all feature importance
    if (length(FeatImp)>20){     # if the number of features is greater than 20, use only 20
      components = 20
    }
    # Append feature importance to output details
    AppendList2CSv(FeatImp,outFile)
    
    
    fitModel <- rf.fit$finalModel
    
    # =============================================================================
    # Assessment of Applicability Domain (plot leverage)
    # =============================================================================
    
    # Residuals
    resids <- fitModel$predicted-fitModel$y # residuals
    write.table("Residuals of the fitted model: ",file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(data.frame(resids), file=outFile,append=T,sep=",",col.names=T,row.names=T, quote=F) # write residuals
    
    predVals.pls.ad <- pred.ts
    Traind.pls= as.matrix(my.datf.train)
    Testd.pls = as.matrix(my.datf.test)
    Hat.train = diag(Traind.pls %*% solve(t(Traind.pls) %*%(Traind.pls), tol=1e-40)  %*% t(Traind.pls))
    Hat.test  = diag(Testd.pls  %*% solve(t(Traind.pls) %*%(Traind.pls), tol=1e-40)  %*% t(Testd.pls))  
    
    
    
    # Leverage / Hat values
    hat.fit <- Hat.test          # hat values
    hat.fit.df <- as.data.frame(hat.fit)    # hat data frame
    hat.mean <- mean(hat.fit)               # mean hat values
    hat.fit.df$warn <- ifelse(hat.fit.df[, 'hat.fit']>3*hat.mean, 'x3',ifelse(hat.fit.df[, 'hat.fit']>2*hat.mean, 'x2', '-' ))
    
    write.table("Leverage output: ",file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(paste("Mean of hat values: ", hat.mean), file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table("Leverage / Hat values with warnings (X3 & X2 = values 3 & 2 times than hat mean): ",file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(hat.fit.df, file=outFile,append=T,sep=",",col.names=T,row.names=T, quote=F) # write hat values and the levels X3, X2 (of hat mean)
    
    #THRESHOLD values: 3m/n, where m is the number of parameters, and n number of observations
    thresh.lever<- (3*(dim(my.datf.train)[2]-1))/dim(my.datf.train)[1] # leverage thresh
    hat.problems<- data.frame(hat.fit[hat.fit>thresh.lever]) # points with high leverage
    
    write.table(paste("Leverage Threshold: ", thresh.lever), file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table("Points with leverage > threshold: ",file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(hat.problems, file=outFile,append=T,sep=",",col.names=T,row.names=T, quote=F)
    
    
    # PDF with 12 plots
    # --------------------------------------------------------------
    pdf(file=paste(outFile,".",sCV,".","split",iSplit,".pdf",sep=""))
    # par(mfrow = c(3, 4)) # all plots into one page!
    
    plot(my.datf.train[,1],pred.tr,xlab="Yobs", ylab="Ypred", type="b", main="Train Yobs-Ypred") # plot 1
    plot(my.datf.test[,1], pred.ts,xlab="Yobs", ylab="Ypred", type="b", main="Test Yobs-Ypred")  # plot 2
    dotchart(as.matrix(FeatImp),main="Feature Importance")                            # plot 3
    
    # Fitted vs Residuals - plot 4
    plot(fitModel$predicted,resids,
         main="Fitted vs. Residuals for Fitted Model",
         xlab="Fitted", ylab="Residuals")
    abline(h = 0, lty = 2)
    
    # Leverage plots - plot 5
    plot(hat.fit, type = "h",
         main="Leverage for Fitted Model",
         xlab="Index", ylab="Hat")
    abline(h = thresh.lever, lty = 2, col="red") # leverage thresh
    
    
    
    # plot(FeatImp, top = components,main="Feature Importance") # ERROR !
    dev.off()
    # --------------------------------------------------------------
    
  }
  return(my.stats)  # return a list with statistics
}


### HELPER FUNCTIONS FOR SVM-RFE
library(caret)
library(kernlab)
source("RFE-svm-reg-function.R")



#====================================================================================================
#  SVM-RFE
# Jose A. Seoane seoane@stanford.edu
# Carlos Fernandez-Lozano carlos.fernandez@udc.es
#====================================================================================================
SVMRFEreg <- function(my.datf.train,my.datf.test,sCV,iSplit=1,fDet=F,outFile="") {
  # ------------------------------------------
  # only Linux or Mac: parallel calculations
  # ------------------------------------------
  #library(doMC)
  #registerDoMC(cores = 2) # CPU cores
  # ------------------------------------------
  # parallel for windows:
  #library(doSNOW)
  #library(foreach)
  #cl<-makeCluster(5) #change the 2 to your number of CPU cores
  #registerDoSNOW(cl)
  
  
  net.c = my.datf.train[,1]   # make available the names of variables from training dataset
  RegrMethod <- "svmRFE" # type of regression
  
  # Define the CV conditions
  ctrl<- trainControl(method=sCV, number=3,repeats=1,
                      summaryFunction=defaultSummary,verboseIter = F)
  
  rfeCtr = rfeControl(functions=svmFuncsGradW,method="cv",number=5,repeats=1, saveDetails = T, verbose=T,rerank = T,allowParallel=T)
  sigma = sigest (as.matrix(my.datf.train[,-1]))[2]
  cs = c(0.0001,0.1,1,5,15,50)  # EXTERNAL PARAMETERS!!!
  cs = c(1,5,15,50)
  eps=c(0.01,0.1,0.3)   # EXTERNAL PARAMETERS!!!
  sizes = 2^(1:sqrt(ncol(my.datf.train)-1))
  tuneVars = expand.grid(.C=cs, .sigma=sigma, .epsilon=eps)    
  
  # Train the model using only training set
  set.seed(iSplit)
  
  rfesvm.fit = rfe(as.matrix(my.datf.train[,-1]),net.c,sizes = sizes,rfeControl=rfeCtr,prob.model =F,method=svmRadialReg,tuneGrid = tuneVars,trControl=ctrl ,allowParallel=T)
  
  # warning in some of the parameters is a extreme value
  if(rfesvm.fit$fit$bestTune$C %in% cs[c(1,length(cs))])
    warning("Best fitted value of C=",rfesvm.fit$fit$bestTune$C," is a extreme value in your possible c values. You may want to reset your C paramenter options", call. = FALSE)
  if(rfesvm.fit$fit$bestTune$epsilon %in% eps[c(1,length(cs))])
    warning("Best fitted value of eps=",rfesvm.fit$fit$bestTune$epsilon," is a extreme value in your possible epsilon values. You may want to reset your eps paramenter options", call. = FALSE)
  
  
  
  
  #------------------------------
  # Training RESULTS
  #------------------------------
  RMSE.tr  <- rfesvm.fit$results[rfesvm.fit$results$Variables ==rfesvm.fit$bestSubset,2]
  R2.tr    <- rfesvm.fit$results[rfesvm.fit$results$Variables ==rfesvm.fit$bestSubset,3]
  if (sCV == "repeatedcv"){ # if 10-fold CV
    RMSEsd.tr <- rfesvm.fit$results[rfesvm.fit$results$Variables ==rfesvm.fit$bestSubset,4]
    R2sd.tr   <- rfesvm.fit$results[rfesvm.fit$results$Variables ==rfesvm.fit$bestSubset,5]
  }
  if (sCV == "LOOCV"){ # if LOOCV
    RMSEsd.tr <- 0 # formulas will be added later!  TODOOOOOOOOOOOOOO
    R2sd.tr   <- 0 # formulas will be added later!
  }
  
  #------------------------------------------------
  # RMSE & R^2, for train/test respectively
  #------------------------------------------------
  rfesvm.train.res <- getTrainPerf(rfesvm.fit$fit)
  rfesvm.test.res  <- postResample(predict(rfesvm.fit,my.datf.test),my.datf.test[,1])
  
  #------------------------------------------------
  # Adj R2, Pearson correlation
  #------------------------------------------------
  pred.tr     <- predict(rfesvm.fit,my.datf.train) # predicted Y for training
  pred.ts     <- predict(rfesvm.fit,my.datf.test)  # predicted Y for test
  noFeats.fit <- length(predictors(rfesvm.fit))    # no. of features from the fitted model
  Feats.fit   <- paste(predictors(rfesvm.fit),collapse="+") # string with the features included in the fitted model
  
  ds.full     <- rbind(my.datf.train,my.datf.test)
  pred.both   <- predict(rfesvm.fit,ds.full)       # predicted Y
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
                   "RMSE.ts" = as.numeric((rf.test.res["RMSE"])),
                   "R2.ts"   = as.numeric((rf.test.res["Rsquared"])),
                   "corP.ts" = as.numeric(corP.ts),
                   "adjR2.both" = as.numeric(adjR2.both),
                   "RMSE.both"  = as.numeric(RMSE.both),
                   "R2.both"    = as.numeric(r2.both))
  
  #---------------------------------------------------------------------
  # Write to file DETAILS for GLM for each cross-validation method
  #---------------------------------------------------------------------
  if (fDet==T) {   # if flag for details if T, print details about any resut
    write("RRegr package | eNanoMapper", file=outFile, append=T)
    write.table(paste("Regression method: ", RegrMethod), file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(paste("Split no.: ", iSplit), file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(paste("CV type: ", sCV),      file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table("Training Set Summary: ",     file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(summary(my.datf.train),       file=outFile,append=T,sep=",",col.names=T,quote=F)
    write.table("Test Set Summary: ",         file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(summary(my.datf.test),        file=outFile,append=T,sep=",",col.names=T,quote=F)   
    
    
    write.table("Predictors: ",     file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(predictors(rfesvm.fit), file=outFile,append=T,sep=",",col.names=T,quote=F)
    
    write.table("Trainig Results: ",     file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(rfesvm.train.res,file=outFile,append=T,sep=",",col.names=T,quote=F)
    write.table("Test Results: ",        file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(rfesvm.test.res, file=outFile,append=T,sep=",",col.names=T,quote=F)
    
    write.table("Full Statistics: ", file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(my.stats,            file=outFile,append=T,sep=",",col.names=T,quote=F)
    
    # Variable Importance (max top 20)
    
    FeatImp <- svmFuncsGradW$rank(rfesvm.fit$fit,as.matrix(ds.full[,-1]),ds.full[,1])
    FeatImp = FeatImp[order(FeatImp[,1],decreasing=T),]
    components = length(FeatImp)  # default plot all feature importance
    if (length(FeatImp)>20){     # if the number of features is greater than 20, use only 20
      components = 20
    }
    # Append feature importance to output details
    AppendList2CSv(FeatImp,outFile)
    
    
    fitModel <- rfesvm.fit$fit
    
    # =============================================================================
    # Assessment of Applicability Domain (plot leverage)
    # =============================================================================
    
    # Residuals
    resids <- pred.both-ds.full[,1]    # residuals
    write.table("Residuals of the fitted model: ",file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(data.frame(resids), file=outFile,append=T,sep=",",col.names=T,row.names=T, quote=F) # write residuals
    
    predVals.pls.ad <- pred.ts
    Traind.pls= as.matrix(my.datf.train)
    Testd.pls = as.matrix(my.datf.test)
    Hat.train = diag(Traind.pls %*% solve(t(Traind.pls) %*%(Traind.pls), tol=1e-40)  %*% t(Traind.pls))
    Hat.test  = diag(Testd.pls  %*% solve(t(Traind.pls) %*%(Traind.pls), tol=1e-40)  %*% t(Testd.pls))  
    
    
    
    # Leverage / Hat values
    hat.fit <- Hat.test          # hat values
    hat.fit.df <- as.data.frame(hat.fit)    # hat data frame
    hat.mean <- mean(hat.fit)               # mean hat values
    hat.fit.df$warn <- ifelse(hat.fit.df[, 'hat.fit']>3*hat.mean, 'x3',ifelse(hat.fit.df[, 'hat.fit']>2*hat.mean, 'x2', '-' ))
    
    write.table("Leverage output: ",file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(paste("Mean of hat values: ", hat.mean), file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table("Leverage / Hat values with warnings (X3 & X2 = values 3 & 2 times than hat mean): ",file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(hat.fit.df, file=outFile,append=T,sep=",",col.names=T,row.names=T, quote=F) # write hat values and the levels X3, X2 (of hat mean)
    
    #THRESHOLD values: 3m/n, where m is the number of parameters, and n number of observations
    thresh.lever<- (3*(dim(my.datf.train)[2]-1))/dim(my.datf.train)[1] # leverage thresh
    hat.problems<- data.frame(hat.fit[hat.fit>thresh.lever]) # points with high leverage
    
    write.table(paste("Leverage Threshold: ", thresh.lever), file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table("Points with leverage > threshold: ",file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(hat.problems, file=outFile,append=T,sep=",",col.names=T,row.names=T, quote=F)
    
    
    # PDF with 12 plots
    # --------------------------------------------------------------
    pdf(file=paste(outFile,".",sCV,".","split",iSplit,".pdf",sep=""))
    # par(mfrow = c(3, 4)) # all plots into one page!
    
    plot(my.datf.train[,1],pred.tr,xlab="Yobs", ylab="Ypred", type="b", main="Train Yobs-Ypred") # plot 1
    plot(my.datf.test[,1], pred.ts,xlab="Yobs", ylab="Ypred", type="b", main="Test Yobs-Ypred")  # plot 2
    fi = as.matrix(FeatImp[,1])
    rownames(fi)=FeatImp[,2]
    dotchart(fi,main="Feature Importance")                            # plot 3
    
    # Fitted vs Residuals - plot 4
    plot(pred.both,resids,
         main="Fitted vs. Residuals for Fitted Model",
         xlab="Fitted", ylab="Residuals")
    abline(h = 0, lty = 2)
    
    # Leverage plots - plot 5
    plot(hat.fit, type = "h",
         main="Leverage for Fitted Model",
         xlab="Index", ylab="Hat")
    abline(h = thresh.lever, lty = 2, col="red") # leverage thresh
    
    
    
    # plot(FeatImp, top = components,main="Feature Importance") # ERROR !
    dev.off()
    # --------------------------------------------------------------
    
  }
  return(my.stats)  # return a list with statistics
}