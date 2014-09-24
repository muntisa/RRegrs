# ======================================================================
# RRegrs - R Regressions
# ======================================================================
# Get the best regression models for one dataset using R caret methods
# eNanoMapper.net
# ----------------------------------------------------------------------
# AUTHORS: 
# ----------------------------------------------------------------------
# Georgia Tsiliki | ChemEng - NTUA, Greece | g_tsiliki [at] hotmail [dot] com
# Cristian R Munteanu | RNASA-IMEDIR, University of A Coruña | muntisa [at] gmail [dot] com
# Jose A. Seoane | Stanford Cancer Institute | seoane [at] stanford [dot] edu
# Carlos Fernandez-Lozano | RNASA-IMEDIR, University of A Coruña | carlos.fernandez [at] udc [dot] es
# Haralambos Sarimveis | ChemEng - NTUA, Greece | hsarimv [at] central [dot] ntua [dot] gr
# Egon Willighagen | BiGCaT - Maastricht University | egon.willighagen [at] gmail [dot] com
# ----------------------------------------------------------------------

library(caret)
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
#----------------------------------------------------------------------------------------------------------------------

AppendList2CSv <- function(l,csvFile) {
  #--------------------------------------------------------------------
  # Write a LIST to CSV file
  #--------------------------------------------------------------------
  out_file <- file(csvFile, open="a")  #creates a file in append mode 
  for (i in seq_along(l)){ 
    # writes the name of the list elements ("A", "B", etc.)
    write.table(names(l)[i],file=out_file,sep=",",dec=".",quote=F,col.names=F,row.names=F) 
    write.table(l[[i]],     file=out_file,sep=",",dec=".",quote=F,col.names=NA,row.names=T) #writes the data.frames 
  } 
  close(out_file)  #close connection to file.csv 
}  
#----------------------------------------------------------------------------------------------------------------------

AppendList2txt <- function(l,csvFile) {
  #--------------------------------------------------------------------
  # Write a LIST to TXT file
  #--------------------------------------------------------------------
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

RemNear0VarCols <- function(ds,fDet=FALSE,outFile="ds3.No0Var.csv") {
  #================================================
  # Removal of near zero variance columns (Step 3)
  #================================================
  # inputs:
  # - ds = dataset frame
  # - fDet = flag for detais (TRUE/FALSE)
  # - outFileName = new file name  (it could include the path)
  
  # output = ds.Rem0NearVar  (ds without columns with near zero variance)
  # if datails = TRUE, output the new ds as a file
  # ------------------------------------------
  
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
#----------------------------------------------------------------------------------------------------------------------

ScalingDS <- function(ds,s=1,c=1,fDet=FALSE,outFileName="ds4.scaled.csv") {
  #===========================
  # Scaling dataset (Step 4)
  #===========================
  # s = { 1,2,3 } - type of scaling: 1 = normalization, 2 = standardization, 3 = other
  # c = the number of column into the dataset to start scaling
  # - if c = 1: included the dependent variable
  # - if c = 2: only the features will be scaled
  
  # fDet = if details need to be printed    (TRUE/FALSE)
  # inFileName  = file name      (it could include the path)
  # outFileName = new file name  (it could include the path)
  #-----------------------------------------------------------------------------------
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
#----------------------------------------------------------------------------------------------------------------------

RemCorrs <- function(ds,fDet,cutoff,outFile) {
  # ========================================
  # Remove the correlated columns (Step 5)
  # ========================================
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
#----------------------------------------------------------------------------------------------------------------------

DsSplit <- function(ds,trainFrac=3/4,fDet=FALSE,PathDataSet="",iSeed) {
  # ===============================================
  # Dataset spliting in Training and Test (Step 6)
  # ===============================================
  # Inputs
  # - ds = frame dataset object
  # - fDet = flag for detais (TRUE/FALSE)
  # - PathDataSet = pathway for results
  
  # Output = training and test datasets (to be used for regressions in other functions)
  # if datails = TRUE, output files will be created
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

LMreg <- function(my.datf.train,my.datf.test,sCV,iSplit=1,fDet=F,outFile="") {
  #==================
  # 8.1. Basic LM
  #==================
  
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


GLMreg <- function(my.datf.train,my.datf.test,sCV,iSplit=1,fDet=F,outFile="") {
  #======================================================
  # 8.2- GLM stepwise regression - based on AIC (caret)
  #======================================================
  # Inputs:
  # - my.datf.train,my.datf.test = training and test dataset frames
  # - sCV = type of cross-validation such as repeatedcv, LOOCV, etc.
  # - iSplit = index of splitalse
  # - fDet = flag for detais (True/F)
  # - outFile = output file for GLM details
  # Output:
  # - list of statistics equal with the header introduced in the main script!
  #   (tr = train, ts = test, both = tr+ts = full dataset)
  # ---------------------------------------------------------------------------------
  
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


PLSreg <- function(my.datf.train,my.datf.test,sCV,iSplit=1,fDet=F,outFile="") {
  #================================
  # 8.3. PLS regression (caret)
  #================================
  
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


LASSOreg <- function(my.datf.train,my.datf.test,sCV,iSplit=1,fDet=F,outFile="") {
  #================================
  # 8.4 Lasso Regression (caret)
  #================================
  
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


RBF_DDAreg <- function(my.datf.train,my.datf.test,sCV,iSplit=1,fDet=F,outFile="") {
  #============================================================
  # 8.5. RBF network with the DDA algorithm regression (caret)
  #============================================================
  
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


SVLMreg <- function(my.datf.train,my.datf.test,sCV,iSplit=1,fDet=F,outFile="") {
  #====================================
  # 8.6 SVM Radial Regression (caret)
  #====================================
  
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


NNreg <- function(my.datf.train,my.datf.test,sCV,iSplit=1,fDet=F,outFile="") {
  #========================================
  # 8.8 Neural Network Regression (caret)
  #========================================
  
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


PLSregWSel <- function(my.datf.train,my.datf.test,sCV,iSplit=1,fDet=F,outFile="") {
  #====================================================================================================
  # 8.3W. PLS regression with filter feature selection (caret)
  #====================================================================================================
  
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


Yrandom<- function(dss,trainFrac,best.reg,best.R2.ts,noYrand,ResBestF){
  #================================================
  # Y-randomization for the best model (Step 12)
  #================================================
  #    - 1 splitting, 1 CV type, best method
  #    - best.R2.ts will be compared with Yrand.R2.ts
  #    - returns ratios DiffsR2/bestR2
  # --------------------------------------------------
  
  cat("-> Best model Y-Randomization ...\n")
  dss[,1] <- sample(dss[,1]) # randomize Y values for the entire dataset
  
  # splitting dataset in training and test
  #---------------------------------------
  Yrand.R2.ts <- NULL     # all values of R2 for each Y randomization
  for (i in 1:noYrand){
    iSeed=i               
    dsList  <- DsSplit(dss,trainFrac,F,PathDataSet,iSeed) # return a list with 2 datasets = dsList$train, dsList$test
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
    if (best.reg=="rf") {  
      my.stats.reg  <- RFreg(ds.train,ds.test,"repeatedcv",i,T,ResBestF) # run NNet for each CV and regr method
    } 
    if (best.reg=="svmRFE") {  
      my.stats.reg  <- SVMRFEreg(ds.train,ds.test,"repeatedcv",i,T,ResBestF) # run NNet for each CV and regr method
    } 
    
    Yrand.R2.ts <- c(Yrand.R2.ts,my.stats.reg$R2.ts) # adding test R2 value Y randomization 
  }
  
  # get histogram for differences between best R2 and the values for each Y randomization
  R2diffs     <- abs(Yrand.R2.ts - best.R2.ts) # absolute differences between R2 values (best model vs Y randomized results)
  R2diffsPerBestR2 <- R2diffs/best.R2.ts        # the same difference in percents
  
  pdf(file=paste(ResBestF,".Yrand.Hist.pdf",sep=""))    # save histogram if ratio diffs R2 into PDF for Y random
  Yrand.hist  <- hist(R2diffsPerBestR2)         # draw histogram of the ratio diffs/Best R2 for Y random
  dev.off()
  
  write.table("Y randomization test: ",file=ResBestF,append=T,sep=",",col.names=F,row.names=F,quote=F)
  write.table("=====================", file=ResBestF,append=T,sep=",",col.names=F,row.names=F,quote=F)
  write.table("Diffs R2 (Best Model - Y rand):",file=ResBestF,append=T,sep=",",col.names=F,row.names=F,quote=F)
  AppendList2CSv(R2diffs, ResBestF)
  write.table("Summary Difs:",file=ResBestF,append=T,sep=",",col.names=F,row.names=F,quote=F)
  AppendList2CSv(summary(R2diffs), ResBestF)
  write.table("Ratio Diffs R2 / Best R2 (Best Model - Y rand):",file=ResBestF,append=T,sep=",",col.names=F,row.names=F,quote=F)
  AppendList2CSv(R2diffsPerBestR2, ResBestF)
  write.table("Summary Difs %:",file=ResBestF,append=T,sep=",",col.names=F,row.names=F,quote=F)
  AppendList2CSv(summary(R2diffsPerBestR2),ResBestF)
  
  return(R2diffsPerBestR2) # return the ratio of diffs with the best R2 from the same 
}

# -----------------------------------------------------------------------
# svm regression function
# -----------------------------------------------------------------------
# jseoane
# use:
# svmFuncsW: regular ranking using w
# svmFuncsLi15:  ranking from eq15 "SVM Feature Selection and Sample Regression for Chinese Medicine Research"
# svmFuncsLi17: ranking from eq17 "SVM Feature Selection and Sample Regression for Chinese Medicine Research"
# svmFuncsGradW: RAKOTOMAMONJY gradient w
# svmFuncsWCor: remove correlated features W. WARNING: computationally expensive

load("model.svmRadialReg.RData")

svmFuncsW = caretFuncs    ## regular ranking using w
svmFuncsW$fit=function(x,y,first,last,...,tuneGrid){
  #cat(param$sigma,"\n")
  library(kernlab)
  sigma = sigest(x)[2]  
  cs = tuneGrid$.C
  eps = tuneGrid$.epsilon
  tuneGrid = expand.grid(.C=cs,.sigma=sigma,.epsilon=eps)  
  train(x,y,...,tuneGrid=tuneGrid)
}
svmFuncsW$rank=function(object,x,y){
  alphas = alpha(object$finalModel)
  alpha.idxs = alphaindex(object$finalModel)
  y.sv = as.numeric(y[alpha.idxs])
  w = (y.sv * alphas) %*% xmatrix(object$finalModel)
  sig = ifelse(object$finalModel@fitted>y,yes=1,no=-1)
  
  avImp = t(w*w)
  out = data.frame(avImp)
  colnames(out) = "Overall"
  out = out[order(out$Overall, decreasing = TRUE), , drop = FALSE]
  out$var <- rownames(out)
  out
}
svmFuncsW$pred= function(object, x)
{
  tmp = predict(object, newdata=x)
  if(object$modelType == "Classification" &
       !is.null(object$modelInfo$prob))
  {
    out1 =cbind(data.frame(pred = tmp),
                as.data.frame(predict(object$finalModel, newdata=x, type = "prob")))
  } else out1 <- tmp
  out1
}

svmFuncsLi15 = svmFuncsW   # ranking from eq15 "SVM Feature Selection and Sample Regression for Chinese Medicine Research"
svmFuncsLi15$rank=function(object,x,y){
  alphas = alpha(object$finalModel)
  alpha.idxs = alphaindex(object$finalModel)
  y.sv = as.numeric(y[alpha.idxs])
  w = (y.sv * alphas) %*% xmatrix(object$finalModel)
  sig = ifelse(object$finalModel@fitted>y,yes=1,no=-1)
  avImp = t(w)* (t(x)%*%sig)  
  out = data.frame(avImp)
  colnames(out) = "Overall"
  out = out[order(out$Overall, decreasing = TRUE), , drop = FALSE]
  out$var <- rownames(out)
  out
}

svmFuncsLi17 = svmFuncsLi15 # ranking from eq17 "SVM Feature Selection and Sample Regression for Chinese Medicine Research"
svmFuncsLi17$rank=function(object,x,y){
  alphas = alpha(object$finalModel)
  alpha.idxs = alphaindex(object$finalModel)
  y.sv = as.numeric(y[alpha.idxs])
  w = (y.sv * alphas) %*% xmatrix(object$finalModel)
  sig = ifelse(object$finalModel@fitted>y,yes=1,no=-1)
  avImp = t(w * (colMeans(x[sig==1,])+colMeans(x[sig==-1,])))
  out = data.frame(avImp)
  colnames(out) = "Overall"
  out = out[order(out$Overall, decreasing = TRUE), , drop = FALSE]
  out$var <- rownames(out)
  out
}

# Based on the gradient of svm coefs
svmFuncsGradW = svmFuncsW
svmFuncsGradW$rank=function(object,x,y){ # RAKOTOMAMONJY gradient w  
  alphas = alpha(object$finalModel)#[[1]]
  alpha.idxs = alphaindex(object$finalModel)#[[1]]
  y.sv = y[alpha.idxs]  
  krnFun = kernelf(object$finalModel)
  kernel = kernelMatrix(krnFun,x)
  sigma = krnFun@kpar$sigma
  xmat = xmatrix(object$finalModel)[[1]]
  kerSV = kernel[alpha.idxs,alpha.idxs]  
  nSV = length(alpha.idxs)
  nfeat = dim(x)[2]
  avImp = numeric(nfeat)
  names(avImp) = colnames(x)
  
  for(i in 1:nfeat){
    deraux =  (  x[alpha.idxs,i] %*% t(as.matrix(rep(1,nSV))) ) -  (as.matrix(rep(1,nSV)) %*% t(x[alpha.idxs,i])     )    
    kernelDeriv1 = -(deraux * kerSV) / (sigma^2)
    kernelDeriv2 =  (deraux * kerSV) / (sigma^2)
    gradMarg1= -t(y.sv*alphas) %*% kernelDeriv1 %*%  (y.sv*alphas)
    gradMarg2= -t(y.sv*alphas) %*% kernelDeriv2 %*%  (y.sv*alphas)
    avImp[i] = gradMarg1^2 + gradMarg2^2
  }
  
  out = data.frame(avImp)
  colnames(out) = "Overall"
  out = out[order(out$Overall, decreasing = TRUE), , drop = FALSE]
  out$var <- rownames(out)
  out
}

# Remove correlated svm coefs
svmFuncsWCor = svmFuncsW  # correlated rfe from "A Greedy Correlation-incorporated SVM-based Algorithm for Gene Selection"
svmFuncsWCor$rank=function(object,x,y){
  alphas = alpha(object$finalModel)
  alpha.idxs = alphaindex(object$finalModel)
  y.sv = as.numeric(y[alpha.idxs])
  w = (y.sv * alphas) %*% xmatrix(object$finalModel)
  
  cormat = cor(x)
  
  diag(cormat)=0
  imp = data.frame(t(w*w))
  imp.s = order(imp,decreasing=T)
  picked = imp.s
  thr = 0.75  
  for(i in 1:length(imp.s)){    
    corvec = cormat[,picked[i]]
    oldpicked = picked
    id.out= which(abs(corvec)>thr)
    id.in = intersect(which(abs(corvec[])<=thr),which(corvec!=0))   # and distinct 0
    len.in = length(id.in)
    if(length(id.in)>0){  
      picked[(i+1):(i+len.in)]= intersect(oldpicked[(i+1):length(oldpicked)],id.in)
    }
    if(length(id.out)>0){    
      picked[(i+len.in+1):length(imp.s)]= intersect(oldpicked[(i+1):length(oldpicked)],id.out)     }    
    cormat[picked[i],]=0        
  }  
  out = imp
  colnames(out) = "Overall"
  out = out[picked, , drop=FALSE]
  out$var <- rownames(out)
  out  
}


RFreg <- function(my.datf.train,my.datf.test,sCV,iSplit=1,fDet=F,outFile="") {
  #======================================
  # Basic RandomForest
  #======================================
  
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


SVMRFEreg <- function(my.datf.train,my.datf.test,sCV,iSplit=1,fDet=F,outFile="") {
  #==========
  # SVM-RFE
  #==========
  
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
  library(kernlab)
  # source("RFE-svm-reg-function.R")
  
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

# MAIN FUNCTION #########################################################################################################

RRegrs<- function(paramFile) { # input = file with all parameters
  # ========================
  # MAIN function RRegrs
  # ========================
  
  #==========================================================================================
  # (1) Load dataset and parameters
  #     (these parameters will be read from an input file! TO BE IMPLEMENTED at the end)
  #==========================================================================================
  # (1.1) PARAMETERS
  #-------------------------
  # Read parameters from a file as data frame
  # -----------------------------------------
  Param.df <- read.csv(paramFile,header=T)
  
  fDet         = as.logical(Param.df[which(Param.df$RRegrs.Parameters=="fDet"),2])         # flag to calculate and print details for all the functions
  fFilters     = as.logical(Param.df[which(Param.df$RRegrs.Parameters=="fFilters"),2])     # flag to apply filters                          (2)
  fScaling     = as.logical(Param.df[which(Param.df$RRegrs.Parameters=="fScaling"),2])     # flag for dataset Scaling                       (3)
  fRemNear0Var = as.logical(Param.df[which(Param.df$RRegrs.Parameters=="fRemNear0Var"),2]) # flag for Removal of near zero variance columns (4)
  fRemCorr     = as.logical(Param.df[which(Param.df$RRegrs.Parameters=="fRemCorr"),2])     # flag for Removal of correlated columns         (5)
  fFeatureSel  = as.logical(Param.df[which(Param.df$RRegrs.Parameters=="fFeatureSel"),2])  # flag for wrapper methods for feature selection (7)
  
  cutoff       = as.numeric(as.character(Param.df[which(Param.df$RRegrs.Parameters=="cutoff"),2]))  # cut off for correlated features
  fLM          = as.logical(Param.df[which(Param.df$RRegrs.Parameters=="fLM"),2])     # flag to run LM            (8.1)
  fGLM         = as.logical(Param.df[which(Param.df$RRegrs.Parameters=="fGLM"),2])    # flag to run GLM           (8.2)
  fPLS         = as.logical(Param.df[which(Param.df$RRegrs.Parameters=="fPLS"),2])    # flag to run PLS           (8.3)
  fLASSO       = as.logical(Param.df[which(Param.df$RRegrs.Parameters=="fLASSO"),2])  # flag to run LASSO         (8.4)
  fRBFdda      = as.logical(Param.df[which(Param.df$RRegrs.Parameters=="fRBFdda"),2]) # flat to run RBF DDA       (8.5)
  fSVLM        = as.logical(Param.df[which(Param.df$RRegrs.Parameters=="fSVLM"),2])   # flat to run svmRadial     (8.6)
  fNN          = as.logical(Param.df[which(Param.df$RRegrs.Parameters=="fNN"),2])     # flat to run NN            (8.8)
  fRF          = as.logical(Param.df[which(Param.df$RRegrs.Parameters=="fRF"),2])      # flag to run RandomForest        (8.9)
  fSVMRFE      = as.logical(Param.df[which(Param.df$RRegrs.Parameters=="fSVMRFE"),2])  # flag to run SVM RFE      (8.10)
  
  # ----------------------------------------------------------------------------------------
  iScaling = as.numeric(as.character(Param.df[which(Param.df$RRegrs.Parameters=="iScaling"),2])) # 1 = normalization; 2 = standardization, 3 = other; any other: no scaling
  iScalCol = as.numeric(as.character(Param.df[which(Param.df$RRegrs.Parameters=="iScalCol"),2])) # 1 = including dependent variable in scaling; 2: only all features; etc.
  # ----------------------------------------------------------------------------------------
  trainFrac   = as.numeric(as.character(Param.df[which(Param.df$RRegrs.Parameters=="trainFrac"),2]))   # the fraction of training set from the entire dataset; trainFrac = the rest of dataset, the test set
  iSplitTimes = as.numeric(as.character(Param.df[which(Param.df$RRegrs.Parameters=="iSplitTimes"),2])) # default is 10; time to split the data in train and test (steps 6-11); report each step + average
  noYrand     = as.numeric(as.character(Param.df[which(Param.df$RRegrs.Parameters=="noYrand"),2]))     # number of Y randomization (default = 100)
  
  CVtypes = strsplit(as.character(Param.df[which(Param.df$RRegrs.Parameters=="CVtypes"),2]),";")[[1]] # types of cross-validation methods
  
  # -------------------------------------------------------------------------------------------------------
  # Files
  # -------------------------------------------------------------------------------------------------------
  PathDataSet    = as.character(Param.df[which(Param.df$RRegrs.Parameters=="PathDataSet"),2])    # dataset folder for input and output files
  DataFileName   = as.character(Param.df[which(Param.df$RRegrs.Parameters=="DataFileName"),2])   # input step 1 = ds original file name
  No0NearVarFile = as.character(Param.df[which(Param.df$RRegrs.Parameters=="No0NearVarFile"),2]) # output step 3 = ds without zero near vars
  ScaledFile     = as.character(Param.df[which(Param.df$RRegrs.Parameters=="ScaledFile"),2])     # output step 4 = scaled ds file name (in the same folder)
  NoCorrFile     = as.character(Param.df[which(Param.df$RRegrs.Parameters=="NoCorrFile"),2])     # output step 5 = dataset after correction removal
  
  ResAvgs        = as.character(Param.df[which(Param.df$RRegrs.Parameters=="ResAvgs"),2])     # the output file with averaged statistics for each regression method
  ResBySplits    = as.character(Param.df[which(Param.df$RRegrs.Parameters=="ResBySplits"),2]) # the output file with statistics for each split and the averaged values
  ResBest        = as.character(Param.df[which(Param.df$RRegrs.Parameters=="ResBest"),2])     # the output file with statistics for the best model
  
  lmFile         = as.character(Param.df[which(Param.df$RRegrs.Parameters=="lmFile"),2])     # LM output file for details
  glmFile        = as.character(Param.df[which(Param.df$RRegrs.Parameters=="glmFile"),2])    # GLM output file for details
  plsFile        = as.character(Param.df[which(Param.df$RRegrs.Parameters=="plsFile"),2])    # PLS output file for details
  lassoFile      = as.character(Param.df[which(Param.df$RRegrs.Parameters=="lassoFile"),2])  # Lasoo Radial output file for details
  rbfDDAFile     = as.character(Param.df[which(Param.df$RRegrs.Parameters=="rbfDDAFile"),2]) # RBF DDA output file for details
  svlmFile       = as.character(Param.df[which(Param.df$RRegrs.Parameters=="svlmFile"),2])   # SVM Radial output file for details
  nnFile         = as.character(Param.df[which(Param.df$RRegrs.Parameters=="nnFile"),2])     # NN Radial output file for details
  rfFile         = as.character(Param.df[which(Param.df$RRegrs.Parameters=="rfFile"),2])      # RF  output file for details
  svmrfeFile     = as.character(Param.df[which(Param.df$RRegrs.Parameters=="svmrfeFile"),2])  # svMRFE  output file for details
  
  # Generate path + file name = original dataset
  inFile <- file.path(PathDataSet, DataFileName)
  
  cat("======================================================================
RRegrs - R Regression Models
Get the best regression models for one dataset using R caret methods
eNanoMapper.net
      
AUTHORS:
Georgia Tsiliki | ChemEng - NTUA, Greece | g_tsiliki [at] hotmail [dot] com
Cristian R Munteanu | RNASA-IMEDIR, University of A Coruña | muntisa [at] gmail [dot] com
Jose A. Seoane | Stanford Cancer Institute | seoane [at] stanford [dot] edu
Carlos Fernandez-Lozano | RNASA-IMEDIR, University of A Coruña | carlos.fernandez [at] udc [dot] es
Haralambos Sarimveis | ChemEng - NTUA, Greece | hsarimv [at] central [dot] ntua [dot] gr
Egon Willighagen | BiGCaT - Maastricht University | egon.willighagen [at] gmail [dot] com
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
  if (fFilters==T) {
    # cat("-> [2] Filtering dataset ... \n")
  }
  
  # -----------------------------------------------------------------------
  # (3) Remove near zero variance columns
  # -----------------------------------------------------------------------
  if (fRemNear0Var==T) {
    cat("-> [3] Removal of near zero variance columns ...\n")
    outFile <- file.path(PathDataSet,No0NearVarFile) # the same folder as input  
    
    # get the ds without near zero cols 
    ds <- cbind("net.c" = ds[,1],RemNear0VarCols(ds[,2:dim(ds)[2]],fDet,outFile))
    # use df without Y (predicted values), reconstruct the ds
    # inputs: ds, flag for details, output file
  }
  
  # -----------------------------------------------------------------------
  # (4) Scaling dataset: normalization (default), standardization, other
  # -----------------------------------------------------------------------
  if (fScaling==T) {
    cat("-> [4] Scaling original dataset ...\n")
    outFile <- file.path(PathDataSet,ScaledFile)       # the same folder as input
    
    # run fuction for scaling input dataset file
    ds <- ScalingDS(ds,iScaling,iScalCol,fDet,outFile)
    # use df without Y (predicted values), reconstruct the ds
    # inputs: ds, type of scaling, flag for details, starting column, output file
  }
  
  # -----------------------------------------------------------------------
  # (5) Remove correlated features
  # -----------------------------------------------------------------------
  if (fRemCorr==T) {    
    cat("-> [5] Removing correlated features ...\n") 
    outFile <- file.path(PathDataSet,NoCorrFile)    # the same folder as the input
    
    # run function to remove the correlations between the features
    ds <- cbind("net.c" = ds[,1],RemCorrs(ds[,2:dim(ds)[2]],fDet,cutoff,outFile))
  }
  
  #=================================================================================================
  # Steps 6 - 11 will be repeated 10 times for reporting each result and average
  #                  (iSplitTimes = 10, default)
  #=================================================================================================
  
  # Initialize the list with the statistics results; the same HEADER as the function output
  dfRes <- list("RegrMeth"     = NULL,
                "Split No"     = NULL,    
                "CVtype"       = NULL,      
                "NoModelFeats" = NULL,
                "ModelFeats"   = NULL,
                "adjR2.tr"  = NULL,
                "RMSE.tr"   = NULL,
                "R2.tr"     = NULL,
                "RMSEsd.tr" = NULL,
                "R2sd.tr"   = NULL,
                "adjR2.ts"= NULL,
                "RMSE.ts" = NULL,
                "R2.ts"   = NULL,
                "corP.ts" = NULL,
                "adjR2.both" = NULL,
                "RMSE.both"  = NULL,
                "R2.both"    = NULL)
  
  #-------------------------------------------------------------------------------------------------
  for (i in 1:iSplitTimes) {                      # Step splitting number = i
    # -----------------------------------------------------------------------
    # (6) Dataset split: Training and Test sets
    # -----------------------------------------------------------------------
    cat("-> [6] Splitting dataset in Training and Test sets ...\n")
    cat(paste("--> Split No.",i,"from",iSplitTimes,"\n"))
    
    iSeed=i                 # to reapeat the ds splitting, different values of seed will be used
    dsList  <- DsSplit(ds,trainFrac,fDet,PathDataSet,iSeed) # return a list with 2 datasets = dsList$train, dsList$test
    # get train and test from the resulted list
    ds.train<- dsList$train
    ds.test <- dsList$test
    
    # Training data set info - to be added later!
    # ------------------------------------------------------------------------
    #nCases     <- dim(ds)[1]                     # number of total cases
    #nCases.tr  <- dim(ds.train)[1]               # number of training cases
    #nFeatures  <- dim(ds.train)[2] - 1           # number of input features (it could be different with the fitted model features!)
    #FeatureList<- paste(names(ds.train)[2:nFeatures], collapse="+") # list of input feature names, from second name, first name is the output variable
    #OutVar     <- names(ds.train)[1]             # name of predicted variable (dependent variable)
    
    # -----------------------------------------------------------------------
    # (7) Feature selection
    # -----------------------------------------------------------------------
    # Two types of regression funtions:
    # -> without wrapper methods (no feature selection excepts the built-in feature selection: Lasso, PLS)
    # -> with wrapper methods (all the functions without built-in feature selection)
    # ===>>>> if fFeatureSel = T, the wrapper versions will be used in step 8
    
    # Note: additional feature selection could be implemented in the future
    
    # -----------------------------------------------------------------------
    # (8) REGRESSION METHODS
    # -----------------------------------------------------------------------
    #
    cat("-> [8] Run Regressions ...\n")
    
    # --------------------------------------------
    # 8.1. Basic LM : default
    # --------------------------------------------
    if (fLM==T) {   # if LM was selected, run the method
      cat("-> [8.1] LM ...\n")
      outFile.LM <- file.path(PathDataSet,lmFile)   # the same folder as the input is used for the output
      
      # Both wrapper and nont-wrapper function are placed in the same external file s8.RegrrMethods.R
      if (fFeatureSel==F) {    # if there is no need of feature selection ->> use normal functions
        
        # For each type of CV do all the statistics
        # -----------------------------------------------------
        for (cv in 1:length(CVtypes)) { # there is no CV but it will be implemented in the future!!!
          my.stats.LM   <- LMreg(ds.train,ds.test,CVtypes[cv],i,fDet,outFile.LM) # run GLM for each CV and regr method
          
          #-------------------------------------------------------
          # Add output from GLM to the list of results
          #-------------------------------------------------------
          # List of results for each splitting, CV type & regression method
          dfRes = mapply(c, my.stats.LM, dfRes, SIMPLIFY=F)
          
        } # end CV types
      } 
      else    # if there is a need for previous feature selection ->> use wrapper functions
      {                     
        # run LM with wrapper method (TO BE IMPLEMENTED!)
      }
      
    } # end LM
    
    # -----------------------------------------------------------------------------------------
    # (8.2) GLM based on AIC regression - Generalized Linear Model with Stepwise Feature Selection
    # -----------------------------------------------------------------------------------------
    if (fGLM==T) {   # if GLM was selected, run the method
      cat("-> [8.2] GLM stepwise - based on AIC ...\n")
      outFile.GLM <- file.path(PathDataSet,glmFile)   # the same folder as the input is used for the output
      
      # Both wrapper and nont-wrapper function are placed in the same external file s8.RegrrMethods.R
      if (fFeatureSel==F) {    # if there is no need of feature selection ->> use normal functions
        
        # List with data set and method information
        #my.stats.dsInfo <- list("RegrMethod"= RegrMethod,           # regression method name
        #"NoCases"= as.numeric(nCases),      # no. of cases
        #"InNoVars"= as.numeric(nFeatures),  # no. of input features
        #"InFeatures"= FeatureList,          # list with feature names
        #"PredVar" = OutVar)                 # name of predicted variable (dependent variable)
        
        # For each type of CV do all the statistics
        # -----------------------------------------------------
        for (cv in 1:length(CVtypes)) {
          my.stats.GLM   <- GLMreg(ds.train,ds.test,CVtypes[cv],i,fDet,outFile.GLM) # run GLM for each CV and regr method
          #my.stats.split <- c(my.stats.dsInfo,my.stats.GLM) # merge the ds info with statistics results for each Cv & reg method
          
          #-------------------------------------------------------
          # Add output from GLM to the list of results
          #-------------------------------------------------------
          # List of results for each splitting, CV type & regression method
          dfRes = mapply(c, my.stats.GLM, dfRes, SIMPLIFY=F)
          
        } # end CV types
      } 
      else    # if there is a need for previous feature selection ->> use wrapper functions
      {                     
        # run GLM with wrapper method (TO BE IMPLEMENTED!)
      }
      
    } # end GLM
    
    # --------------------------------------------
    # 8.3. PLS
    # --------------------------------------------
    if (fPLS==T) {   # if LASSO was selected, run the method
      outFile.PLS <- file.path(PathDataSet,plsFile)   # the same folder as the input is used for the output
      
      # Both wrapper and nont-wrapper function are placed in the same external file s8.RegrrMethods.R
      
      if (fFeatureSel==F) {    # if there is no need of feature selection ->> use normal functions
        cat("-> [8.3] PLS ...\n")
        # For each type of CV do all the statistics
        # -----------------------------------------------------
        for (cv in 1:length(CVtypes)) {
          my.stats.PLS  <- PLSreg(ds.train,ds.test,CVtypes[cv],i,fDet,outFile.PLS) # run SVLM Radial for each CV and regr method
          #-------------------------------------------------------
          # Add output from GLM to the list of results
          #-------------------------------------------------------
          # List of results for each splitting, CV type & regression method
          dfRes = mapply(c, my.stats.PLS, dfRes, SIMPLIFY=F)
        } # end CV types
      } 
      else    # if there is a need for previous feature selection ->> use wrapper functions
      {                     
        cat("-> [8.3] PLS Wrapper Feature Selection ...\n")
        # For each type of CV do all the statistics
        # -----------------------------------------------------
        for (cv in 1:length(CVtypes)) {
          my.stats.PLS  <- PLSregWSel(ds.train,ds.test,CVtypes[cv],i,fDet,outFile.PLS) # run SVLM Radial for each CV and regr method
          #-------------------------------------------------------
          # Add output from GLM to the list of results
          #-------------------------------------------------------
          # List of results for each splitting, CV type & regression method
          dfRes = mapply(c, my.stats.PLS, dfRes, SIMPLIFY=F)
        } # end CV types 
      }
      
    } # end PLS with wrapper
    
    # --------------------------------------------
    # 8.4. LASSO regression
    # --------------------------------------------
    if (fLASSO==T) {   # if LASSO was selected, run the method
      cat("-> [8.4] LASSO ...\n")
      outFile.LASSO <- file.path(PathDataSet,lassoFile)   # the same folder as the input is used for the output
      
      # Both wrapper and nont-wrapper function are placed in the same external file s8.RegrrMethods.R
      if (fFeatureSel==F) {    # if there is no need of feature selection ->> use normal functions
        # For each type of CV do all the statistics
        # -----------------------------------------------------
        for (cv in 1:length(CVtypes)) {
          my.stats.LASSO  <- LASSOreg(ds.train,ds.test,CVtypes[cv],i,fDet,outFile.LASSO) # run SVLM Radial for each CV and regr method
          #-------------------------------------------------------
          # Add output from GLM to the list of results
          #-------------------------------------------------------
          # List of results for each splitting, CV type & regression method
          dfRes = mapply(c, my.stats.LASSO, dfRes, SIMPLIFY=F)
        } # end CV types
      } 
      else    # if there is a need for previous feature selection ->> use wrapper functions
      {                     
        # run Lasso with wrapper method (TO BE IMPLEMENTED!)
      }
      
    } # end Lasso
    
    # --------------------------------------------
    # 8.5. RBF network with the DDA algorithm regression (caret)
    # --------------------------------------------
    if (fRBFdda==T) {   # if SVM Radial was selected, run the method
      cat("-> [8.5] RBF network with the DDA ...\n")
      outFile.rbfDDA <- file.path(PathDataSet,rbfDDAFile)   # the same folder as the input is used for the output
      
      # Both wrapper and nont-wrapper function are placed in the same external file s8.RegrrMethods.R
      if (fFeatureSel==F) {    # if there is no need of feature selection ->> use normal functions
        # For each type of CV do all the statistics
        # -----------------------------------------------------
        for (cv in 1:length(CVtypes)) {
          my.stats.rbfDDA  <- RBF_DDAreg(ds.train,ds.test,CVtypes[cv],i,fDet,outFile.rbfDDA) # run SVLM Radial for each CV and regr method
          #-------------------------------------------------------
          # Add output from SVM Radial to the list of results
          #-------------------------------------------------------
          # List of results for each splitting, CV type & regression method
          dfRes = mapply(c, my.stats.rbfDDA, dfRes, SIMPLIFY=F)
        } # end CV types
      } 
      else    # if there is a need for previous feature selection ->> use wrapper functions
      {                     
        # run rbfDDA with wrapper method (TO BE IMPLEMENTED!)
      }
      
    } # end rbfDDA
    
    # --------------------------------------------
    # 8.6. SVM radial regression
    # --------------------------------------------
    if (fSVLM==T) {   # if SVM Radial was selected, run the method
      cat("-> [8.6] SVM radial ...\n")
      outFile.SVLM <- file.path(PathDataSet,svlmFile)   # the same folder as the input is used for the output
      
      # Both wrapper and nont-wrapper function are placed in the same external file s8.RegrrMethods.R
      if (fFeatureSel==F) {    # if there is no need of feature selection ->> use normal functions
        # For each type of CV do all the statistics
        # -----------------------------------------------------
        for (cv in 1:length(CVtypes)) {
          my.stats.SVLM  <- SVLMreg(ds.train,ds.test,CVtypes[cv],i,fDet,outFile.SVLM) # run SVLM Radial for each CV and regr method
          #-------------------------------------------------------
          # Add output from SVM Radial to the list of results
          #-------------------------------------------------------
          # List of results for each splitting, CV type & regression method
          dfRes = mapply(c, my.stats.SVLM, dfRes, SIMPLIFY=F)
        } # end CV types
      } 
      else    # if there is a need for previous feature selection ->> use wrapper functions
      {                     
        # run SVLM with wrapper method (TO BE IMPLEMENTED!)
      }
      
    } # end SVLM
    
    # --------------------------------------------
    # 8.7. SVM linear
    # --------------------------------------------
    
    # --------------------------------------------
    # 8.8. Neural Networks Regression
    # --------------------------------------------
    if (fNN==T) {   # if NNet was selected, run the method
      cat("-> [8.8] Neural Networks ...\n")
      outFile.NN <- file.path(PathDataSet,nnFile)   # the same folder as the input is used for the output
      
      # Both wrapper and nont-wrapper function are placed in the same external file s8.RegrrMethods.R
      if (fFeatureSel==F) {    # if there is no need of feature selection ->> use normal functions
        # For each type of CV do all the statistics
        # -----------------------------------------------------
        for (cv in 1:length(CVtypes)) {
          my.stats.NN  <- NNreg(ds.train,ds.test,CVtypes[cv],i,fDet,outFile.NN) # run NNet for each CV and regr method
          #-------------------------------------------------------
          # Add output from NNet to the list of results
          #-------------------------------------------------------
          # List of results for each splitting, CV type & regression method
          dfRes = mapply(c, my.stats.NN, dfRes, SIMPLIFY=F)
        } # end CV types
      } 
      else    # if there is a need for previous feature selection ->> use wrapper functions
      {                     
        # run NNet with wrapper method (TO BE IMPLEMENTED!)
      }
      
    } # end NNet
    # --------------------------------------------
    # 8.9. RF
    # --------------------------------------------
    if (fRF==T) {   # if NNet was selected, run the method
      cat("-> [8.9] Random Forest ...\n")
      outFile.RF <- file.path(PathDataSet,rfFile)   # the same folder as the input is used for the output
      0
      if (fFeatureSel==F) {    # if there is no need of feature selection ->> use normal functions
        # For each type of CV do all the statistics
        # -----------------------------------------------------
        for (cv in 1:length(CVtypes)) {
          my.stats.RF  <- RFreg(ds.train,ds.test,CVtypes[cv],i,fDet,outFile.RF) # run NNet for each CV and regr method
          #-------------------------------------------------------
          # Add output from NNet to the list of results
          #-------------------------------------------------------
          # List of results for each splitting, CV type & regression method
          dfRes = mapply(c, my.stats.RF, dfRes, SIMPLIFY=F)
        } # end CV types
      } 
      else    # if there is a need for previous feature selection ->> use wrapper functions
      {                     
        # run RFet with wrapper method (TO BE IMPLEMENTED!)
      }
      
    } # end RF
    # --------------------------------------------
    # 8.10. SVM RFE
    # --------------------------------------------
    if (fSVMRFE==T) {   # if NNet was selected, run the method
      cat("-> [8.10] SVM RFE ...\n")
      outFile.SVMRFE <- file.path(PathDataSet,svmrfeFile)   # the same folder as the input is used for the output
      
      # Both wrapper and nont-wrapper function are placed in the same external file s8.RegrrMethods.R
      if (fFeatureSel==F) {    # if there is no need of feature selection ->> use normal functions
        # For each type of CV do all the statistics
        # -----------------------------------------------------
        for (cv in 1:length(CVtypes)) {
          my.stats.SVMRFE  <- SVMRFEreg(ds.train,ds.test,CVtypes[cv],i,fDet,outFile.SVMRFE) # run SVM RFEet for each CV and regr method
          #-------------------------------------------------------
          # Add output from SVM RFE to the list of results
          #-------------------------------------------------------
          # List of results for each splitting, CV type & regression method
          dfRes = mapply(c, my.stats.SVMRFE, dfRes, SIMPLIFY=F)
        } # end CV types
      } 
      else    # if there is a need for previous feature selection ->> use wrapper functions
      {                     
        # run RFet with wrapper method (TO BE IMPLEMENTED!)
      }
      
    } # end SVM RFE
    # END OF REGRESSION Functions !!!
  }
  
  #------------------------------------------------------------------------------
  # 9. Results for all splittings (not ordered)
  #-------------------------------------------------------------------------------
  cat("[12] Results for all splitings ...\n")
  df.res <- data.frame(dfRes)
  print(df.res) # print all results as data frame
  
  # Writing the statistics into output files: one with detailed splits, other with only averages
  # File names includin paths for the statistics outputs (only averages and split detailed; +averages [to be implemented])
  ResBySplitsF <- file.path(PathDataSet,ResBySplits)   # the main output file with statistics for each split
  write.csv(df.res, file = ResBySplitsF)    # write statistics data frame into a CSV output file
  # file.show(ResBySplitsF)  # show the statistics file!
  
  #-------------------------------------------------------------------------------------
  # Averaged values of the results by each Regression Method & CV type
  #-------------------------------------------------------------------------------------
  cat("-> Averaged statistics ...\n")
  
  ResAvgsF <- file.path(PathDataSet,ResAvgs)           # the main output file with averaged statistics for each regression method
  library(data.table)
  dt.res  <- data.table(df.res) # convert data frame into data table (for sorting abd averaging)
  
  # MEANS for each Regression Method & CV type
  #--------------------------------------------------------------------------------------------------------------
  # means for all CV types, not only 10CV
  dt.mean <- dt.res[,list(adjR2.tr.Avg=mean(adjR2.tr),RMSE.tr.Avg=mean(RMSE.tr),R2.tr.Avg=mean(R2.tr),
                          RMSEsd.tr.Avg=mean(RMSEsd.tr),R2sd.tr.Avg=mean(R2sd.tr),adjR2.ts.Avg=mean(adjR2.ts),
                          RMSE.ts.Avg=mean(RMSE.ts),R2.ts.Avg=mean(R2.ts),corP.ts.Avg=mean(corP.ts),
                          adjR2.both.Avg=mean(adjR2.both),RMSE.both.Avg=mean(RMSE.both),
                          R2.both.Avg=mean(R2.both)),by="RegrMeth,CVtype"]
  
  dt.mean     <- dt.mean[dt.mean$CVtype=="repeatedcv",]   # keep only the 10CV results to be used to find the best model
  dt.mean.ord <- dt.mean[order(-rank(adjR2.ts.Avg))]      # descendent order the averages by adjR2.ts.Avg
  
  # Write averages descendent ordered by adjR2.ts.Avg
  #-------------------------------------------------------------------------------
  write.csv(data.frame(dt.mean.ord), file = ResAvgsF)    # write statistics data frame into a CSV output file
  
  #------------------------------------------------------------------------------
  # 10. Best model selection - detailed statistics
  #-------------------------------------------------------------------------------
  cat("-> Best model analysis ...\n")
  
  # ADD an algorithm to verifty similar adjR2 values
  # From the best ones (+/- 0.05 of adjR2), chose the one with less variables, after that the one with min RMSE!!!
  
  best.dt  <- dt.mean.ord[1] # the best model (adjR2.ts) should be the first value in the descendent ordered results
  # best.reg <- paste(best.dt$RegrMeth,collapse="") # best regrression method
  
  # New conditions
  # +/- 0.05 adjR2ts --> min(RMSE)
  best.adjR2.ts <- as.numeric(data.frame(best.dt)[,8]) # best adjR2.ts avgs
  
  
  # best model with adjR2.ts +/- 0.05 and min of RMSE for Avgs
  best.dt  <- dt.mean.ord[adjR2.ts.Avg %between% c(best.adjR2.ts-0.05,best.adjR2.ts+0.05)][RMSE.ts.Avg == min(RMSE.ts.Avg)]
  best.reg <- paste(best.dt$RegrMeth,collapse="") # best regrression method
  
  # best model non-averaged ? no. of features
  # best.method <- dt.res[CVtype == "repeatedcv"][RegrMeth == best.reg] # best modes corresponding with the avg best values
  # best.method.mean <- mean(as.numeric(data.frame(best.method)[,11])) # best adjR2.ts)
  # dt.res[CVtype == "repeatedcv"][RegrMeth == best.reg][NoModelFeats == min(NoModelFeats)][RMSE.ts == min(RMSE.ts)]
  
  #----------------------------------------------------------
  # 11. Best model detailed statistics 
  #----------------------------------------------------------
  # Write the best model statistics
  ResBestF <- file.path(PathDataSet,ResBest)
  write.table("Averaged values for all spits: ",file=ResBestF,append=T,sep=",",col.names=F,row.names=F,quote=F)
  # write.csv(data.frame(best.dt), file = ResBestF)    # write statistics data frame into a CSV output file
  write.table(data.frame(best.dt), file=ResBestF,append=T,sep=",",col.names=T,quote=F) # write statistics data frame into a CSV output file
  
  # Use the last split for dataset (ds.train & ds.test) ! (or chose other one?)
  
  # Run the caret function with the method from the best method, for one training-test split only
  # and append the details in the best model output file
  
  if (best.reg=="lm") {
    my.stats.reg  <- LMreg(ds.train,ds.test,"repeatedcv",i,T,ResBestF) # run GLM for each CV and regr method
  }
  if (best.reg=="glmStepAIC") {
    my.stats.reg  <- GLMreg(ds.train,ds.test,"repeatedcv",i,T,ResBestF) # run GLM for each CV and regr method
  }
  if (best.reg=="pls") {
    my.stats.reg  <- PLSreg(ds.train,ds.test,"repeatedcv",i,T,ResBestF) # run SVLM Radial for each CV and regr method
  }
  if (best.reg=="lasso") {
    my.stats.reg  <- LASSOreg(ds.train,ds.test,"repeatedcv",i,T,ResBestF) # run SVLM Radial for each CV and regr method
  }
  if (best.reg=="rbfDDA") {  
    my.stats.reg  <- RBF_DDAreg(ds.train,ds.test,"repeatedcv",i,T,ResBestF) # run SVLM Radial for each CV and regr method
  }
  if (best.reg=="svmRadial") {  
    my.stats.reg  <- SVLMreg(ds.train,ds.test,"repeatedcv",i,T,ResBestF) # run SVLM Radial for each CV and regr method
  }
  if (best.reg=="nnet") {  
    my.stats.reg  <- NNreg(ds.train,ds.test,"repeatedcv",i,T,ResBestF) # run NNet for each CV and regr method
  } 
  if (best.reg=="rf") {  
    my.stats.reg  <- RFreg(ds.train,ds.test,"repeatedcv",i,T,ResBestF) # run NNet for each CV and regr method
  } 
  if (best.reg=="svmRFE") {  
    my.stats.reg  <- SVMRFEreg(ds.train,ds.test,"repeatedcv",i,T,ResBestF) # run NNet for each CV and regr method
  } 
  
  #------------------------------------------------------------------------------
  # 12. Test best model with test dataset
  #                   (+ Y randomization 100 times, bootstaping)
  #-------------------------------------------------------------------------------
  # ratios Yrand R2 - Best model R2 / Best model R2
  R2Diff.Yrand <- Yrandom(ds,trainFrac,best.reg,my.stats.reg$R2.ts,noYrand,ResBestF) # mean value of ratio (deatails are printed to output file)
  
  #------------------------------------------------------------------------------
  # Assessment of Applicability Domain (plot leverage) was included as details
  # in lm and glm methods
  #-------------------------------------------------------------------------------
  
  #----------------------------
  # Indicate main result files
  #----------------------------
  cat("\n MAIN RESULT FILES:\n")
  cat("===================================\n")
  cat("Statistics for each data set splitting/method/CV type:", ResBySplits,"\n")
  cat("Averages for all data set splittings by method/CV type:",ResAvgsF,"\n")
  cat("Best model statistics:",ResBestF,"\n")
  cat("Best model plots:",paste(ResBestF,".repeatedcv.split",i,".pdf"),"\n")
  cat("Best model Y-randomization plot:",paste(ResBestF,".Yrand.Hist.pdf"),"\n")
  cat("\n* if you choose Details, additional CSV files will be create for each method.\n")
  
  return(ResBestF) # return the statistics for the best model (all the other info could be found in output files)
}
