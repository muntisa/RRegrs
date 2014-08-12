# ======================================================================
# Regression function based on caret package
# ----------------------------------------------------------------------
# contact: Cristian R Munteanu | BiGCaT - UM    | muntisa@gmail.com
#          Georgia Tsiliki     | ChemEng - NTUA | g_tsiliki@hotmail.com


library(caret)
#----------------------------------------------------------------------------------------------------------------------
# General functions
#----------------------------------------------------------------------------------------------------------------------
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

#--------------------------------------------------------------------
# Write a LIST to TXT file
#--------------------------------------------------------------------
AppendList2txt <- function(l,csvFile) {
  out_file <- file(csvFile, open="a")  #creates a file in append mode 
  for (i in seq_along(l)){ 
    write.table(names(l)[i], file=out_file, sep=" ", dec=".", 
                quote=FALSE, col.names=FALSE, row.names=FALSE)  #writes the name of the list elements ("A", "B", etc) 
    write.table(l[[i]], file=out_file, sep=" ", dec=".", quote=FALSE, 
                col.names=NA, row.names=TRUE)  #writes the data.frames 
  } 
  close(out_file)  #close connection to file.csv 
}  

#====================================================================================================
# 8.1. Basic LM no CV
#====================================================================================================
LMreg <- function(my.datf.train,my.datf.test,sCV,iSplit=1,fDet=FALSE,outFile="") {
  net.c = my.datf.train[,1]   # make available the names of variables from training dataset
  RegrMethod <- "lm" # type of regression
  
  # Define the CV conditions
  ctrl<- trainControl(method = sCV, number = 10,repeats = 10,
                      summaryFunction = defaultSummary)
  
  # Train the model using only training set
  set.seed(iSplit)
  X <- as.matrix(my.datf.train[,2:dim(my.datf.train)[2]])
  Y <- as.matrix(my.datf.train[,1])
  lm.fit<- lm(Y~X)
  lm.fit.sum<- summary(lm.fit)
  
  #------------------------------
  # Training RESULTS
  #------------------------------
  RMSE.tr  <- 0 # to be calculated
  R2.tr    <- lm.fit.sum$r.squared
  
  RMSEsd.tr <- 0 # formulas will be added later!
  R2sd.tr   <- 0 # formulas will be added later!
  
  #------------------------------------------------
  # RMSE & R^2, for train/test respectively
  #------------------------------------------------
  lm.train.res <- 0
  lm.test.res  <- 0
  
  #------------------------------------------------
  # Adj R2, Pearson correlation
  #------------------------------------------------
  pred.tr     <- 0 # predicted Y
  pred.ts     <- 0  # predicted Y
  noFeats.fit <- dim(lm.fit.sum$coefficients)[1] - 1    # no. of features from the fitted model
  Feats.fit   <- paste(names(my.datf.train)[-1],collapse = "+") # string with the features included in the fitted model (all the training features!)
  
  ds.full     <- rbind(my.datf.train,my.datf.test)
  pred.both   <- 0       # predicted Y
  adjR2.tr    <- lm.fit.sum$adj.r.squared
  adjR2.ts    <- 0
  corP.ts     <- 0
  
  adjR2.both  <- 0
  RMSE.both   <- 0
  r2.both     <- 0
  # zero values should be corrected!!!
  
  
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
                   "RMSE.ts" = 0,
                   "R2.ts"   = 0,
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
    
    # TO BE MODIFIED!
    #write.table("Fitting results: ", file = outFile,append = TRUE, sep = " ",col.names = FALSE,quote = FALSE)
    #AppendList2txt(lm.fit.sum,outFile)
    
    write.table("Full Statistics: ", file = outFile,append = TRUE, sep = " ",col.names = FALSE,quote = FALSE)
    write.table(my.stats, file = outFile,append = TRUE, sep = " ",col.names = TRUE,quote = FALSE)
  }
  
  return(my.stats)  # return a list with statistics
}

#====================================================================================================
# 8.2- GLM stepwise regression - based on AIC (caret)
#====================================================================================================
# Inputs:
# - my.datf.train,my.datf.test = training and test dataset frames
# - sCV = type of cross-validation such as repeatedcv, LOOCV, etc.
# - iSplit = index of split
# - fDet = flag for detais (TRUE/FALSE)
# - outFile = output file for GLM details
# Output:
# - list of statistics equal with the header introduced in the main script!
#   (tr = train, ts = test, both = tr+ts = full dataset)
# ---------------------------------------------------------------------------------------------------
GLMreg <- function(my.datf.train,my.datf.test,sCV,iSplit=1,fDet=FALSE,outFile="") {
  library(caret)
  #attach(my.datf.train)    # make available the names of variables from training dataset
  net.c = my.datf.train[,1] # dependent variable is the first column in Training set
  RegrMethod <- "glmStepAIC" # type of regression
  
  # Define the CV conditions
  ctrl<- trainControl(method = sCV, number = 10,repeats = 10,
                      summaryFunction = defaultSummary)
  
  # Train the model using only training set
  set.seed(iSplit)
  glm.fit<- train(net.c~.,data=my.datf.train,
                  method = 'glmStepAIC', tuneLength = 10, trControl = ctrl,
                  metric = 'RMSE')
  
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
    write.table(data.frame(summary(glm.fit)$coefficients), file = outFile,append = TRUE, sep = " ",col.names = TRUE,quote = FALSE)
    
    write.table("Predictors: ", file = outFile,append = TRUE, sep = " ",col.names = FALSE,quote = FALSE)
    write.table(predictors(glm.fit), file = outFile,append = TRUE, sep = " ",col.names = TRUE,quote = FALSE)
    
    write.table("Trainig Results: ", file = outFile,append = TRUE, sep = " ",col.names = FALSE,quote = FALSE)
    write.table(predictors(lm.train.res), file = outFile,append = TRUE, sep = " ",col.names = TRUE,quote = FALSE)
    write.table("Test Results: ", file = outFile,append = TRUE, sep = " ",col.names = FALSE,quote = FALSE)
    write.table(predictors(lm.test.res), file = outFile,append = TRUE, sep = " ",col.names = TRUE,quote = FALSE)
    
    write.table("Full Statistics: ", file = outFile,append = TRUE, sep = " ",col.names = FALSE,quote = FALSE)
    write.table(my.stats, file = outFile,append = TRUE, sep = " ",col.names = TRUE,quote = FALSE)
  }
  
  # my.stats.full <- c(my.stats.dsInfo,my.stats.10CV,my.stats.LOOCV)   # merge the CV results into one list that contains the names of each field!
  
  return(my.stats)  # return a list with statistics
}

#====================================================================================================
# 8.3. PLS regression (caret)
#====================================================================================================
PLSreg <- function(my.datf.train,my.datf.test,sCV,iSplit=1,fDet=FALSE,outFile="") {
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
  
  return(my.stats)  # return a list with statistics
}

#====================================================================================================
# 8.3W. PLS regression with wrapper feature selection (caret)
#====================================================================================================
PLSregWSel <- function(my.datf.train,my.datf.test,sCV,iSplit=1,fDet=FALSE,outFile="") {
  library(caret)
  
  net.c = my.datf.train[,1] # dependent variable is the first column in Training set
  RegrMethod <- "pls.WSel" # type of regression
  
  # Define the CV conditions
  ctrlw <- rfeControl(method = 'boot', number = 25,saveDetails=T)
  ctrl  <- trainControl(method = sCV, number = 10,repeats = 1,
                        summaryFunction = defaultSummary,savePred=T)

  subsetsx<- seq(5,dim(my.datf.train)[2]-1, by = 10)
  
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
    
    write.table("NNet variable importance: ", file = outFile,append = TRUE, sep = " ",col.names = FALSE,quote = FALSE)
    AppendList2txt(varImp(pls.fit),outFile)
    
    write.table("Full Statistics: ", file = outFile,append = TRUE, sep = " ",col.names = FALSE,quote = FALSE)
    write.table(my.stats, file = outFile,append = TRUE, sep = " ",col.names = TRUE,quote = FALSE)
  }
  
  return(my.stats)  # return a list with statistics
}


#====================================================================================================
# 8.4 Lasso Regression (caret)
#====================================================================================================
LASSOreg <- function(my.datf.train,my.datf.test,sCV,iSplit=1,fDet=FALSE,outFile="") {
  library(caret)
  #attach(my.datf.train)   # make available the names of variables from training dataset
  net.c = my.datf.train[,1] # dependent variable is the first column in Training set
  RegrMethod <- "lasso.RMSE" # type of regression
  
  # Define the CV conditions
  ctrl<- trainControl(method = sCV, number = 10,repeats = 10,
                      summaryFunction = defaultSummary)
  
  # Train the model using only training set
  set.seed(iSplit)
  # Train the model using only training set
  set.seed(iSplit)
  las.fit<- train(net.c~.,data=my.datf.train,
                  method = 'lasso', tuneLength = 10, trControl = ctrl,
                  metric = 'RMSE',tuneGrid=expand.grid(.fraction= seq(0.1,1,by=0.1))) # preProc = c('center', 'scale')
  
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
    write.table(predictors(las.fit), file = outFile,append = TRUE, sep = " ",col.names = TRUE,quote = FALSE)
    
    write.table("Trainig Results: ", file = outFile,append = TRUE, sep = " ",col.names = FALSE,quote = FALSE)
    write.table(predictors(lm.train.res), file = outFile,append = TRUE, sep = " ",col.names = TRUE,quote = FALSE)
    write.table("Test Results: ", file = outFile,append = TRUE, sep = " ",col.names = FALSE,quote = FALSE)
    write.table(predictors(lm.test.res), file = outFile,append = TRUE, sep = " ",col.names = TRUE,quote = FALSE)
    
    write.table("Full Statistics: ", file = outFile,append = TRUE, sep = " ",col.names = FALSE,quote = FALSE)
    write.table(my.stats, file = outFile,append = TRUE, sep = " ",col.names = TRUE,quote = FALSE)
  }
  
  return(my.stats)  # return a list with statistics
}

#====================================================================================================
# 8.5. RBF network with the DDA algorithm regression (caret)
#====================================================================================================
RBF_DDAreg <- function(my.datf.train,my.datf.test,sCV,iSplit=1,fDet=FALSE,outFile="") {
  library(caret)
  #attach(my.datf.train)   # make available the names of variables from training dataset
  net.c = my.datf.train[,1] # dependent variable is the first column in Training set
  RegrMethod <- "rbfDDA" # type of regression
  
  # Define the CV conditions
  ctrl<- trainControl(method = sCV, number = 10,repeats = 10,
                      summaryFunction = defaultSummary)
  
  # Train the model using only training set
  set.seed(iSplit)
  
  rbf.fit<- train(net.c~.,data=my.datf.train,
                  method = 'rbfDDA',trControl = ctrl,
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
    write.table(predictors(rbf.fit), file = outFile,append = TRUE, sep = " ",col.names = TRUE,quote = FALSE)
    
    write.table("Trainig Results: ", file = outFile,append = TRUE, sep = " ",col.names = FALSE,quote = FALSE)
    write.table(predictors(lm.train.res), file = outFile,append = TRUE, sep = " ",col.names = TRUE,quote = FALSE)
    write.table("Test Results: ", file = outFile,append = TRUE, sep = " ",col.names = FALSE,quote = FALSE)
    write.table(predictors(lm.test.res), file = outFile,append = TRUE, sep = " ",col.names = TRUE,quote = FALSE)
    
    write.table("Full Statistics: ", file = outFile,append = TRUE, sep = " ",col.names = FALSE,quote = FALSE)
    write.table(my.stats, file = outFile,append = TRUE, sep = " ",col.names = TRUE,quote = FALSE)
  }
  
  return(my.stats)  # return a list with statistics
}

#====================================================================================================
# 8.6 SVM Radial Regression (caret)
#====================================================================================================
SVLMreg <- function(my.datf.train,my.datf.test,sCV,iSplit=1,fDet=FALSE,outFile="") {
  library(caret)
  #attach(my.datf.train)   # make available the names of variables from training dataset
  net.c = my.datf.train[,1] # dependent variable is the first column in Training set
  RegrMethod <- "svmRadial.RMSE" # type of regression
  
  # Define the CV conditions
  ctrl<- trainControl(method = sCV, number = 10,repeats = 10,
                      summaryFunction = defaultSummary)

  # Train the model using only training set
  set.seed(iSplit)
  
  svmL.fit<- train(net.c~.,data=my.datf.train,
                   method = 'svmRadial', tuneLength = 10, trControl = ctrl,
                   metric = 'RMSE',
                   tuneGrid=expand.grid(.sigma=seq(0,1,0.1),.C= c(1:10))) # preProc = c('center', 'scale')
  
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
    write.table(predictors(svmL.fit), file = outFile,append = TRUE, sep = " ",col.names = TRUE,quote = FALSE)
    
    write.table("Trainig Results: ", file = outFile,append = TRUE, sep = " ",col.names = FALSE,quote = FALSE)
    write.table(predictors(lm.train.res), file = outFile,append = TRUE, sep = " ",col.names = TRUE,quote = FALSE)
    write.table("Test Results: ", file = outFile,append = TRUE, sep = " ",col.names = FALSE,quote = FALSE)
    write.table(predictors(lm.test.res), file = outFile,append = TRUE, sep = " ",col.names = TRUE,quote = FALSE)
    
    write.table("Full Statistics: ", file = outFile,append = TRUE, sep = " ",col.names = FALSE,quote = FALSE)
    write.table(my.stats, file = outFile,append = TRUE, sep = " ",col.names = TRUE,quote = FALSE)
  }
  
  return(my.stats)  # return a list with statistics
}

#====================================================================================================
# 8.8 Neural Network Regression (caret)
#====================================================================================================
NNreg <- function(my.datf.train,my.datf.test,sCV,iSplit=1,fDet=FALSE,outFile="") {
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
                 linout=TRUE, trace = FALSE,MaxNWts=20000,
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
    write.table(predictors(nn.fit), file = outFile,append = TRUE, sep = " ",col.names = TRUE,quote = FALSE)
    
    write.table("Trainig Results: ", file = outFile,append = TRUE, sep = " ",col.names = FALSE,quote = FALSE)
    write.table(predictors(lm.train.res), file = outFile,append = TRUE, sep = " ",col.names = TRUE,quote = FALSE)
    write.table("Test Results: ", file = outFile,append = TRUE, sep = " ",col.names = FALSE,quote = FALSE)
    write.table(predictors(lm.test.res), file = outFile,append = TRUE, sep = " ",col.names = TRUE,quote = FALSE)
    
    #write.table("NNet fitting resuts: ", file = outFile,append = TRUE, sep = " ",col.names = FALSE,quote = FALSE)
    #AppendList2txt(list(nn.fit),outFile)
    
    write.table("NNet variable importance: ", file = outFile,append = TRUE, sep = " ",col.names = FALSE,quote = FALSE)
    AppendList2txt(varImp(nn.fit),outFile)
    
    write.table("Full Statistics: ", file = outFile,append = TRUE, sep = " ",col.names = FALSE,quote = FALSE)
    write.table(my.stats, file = outFile,append = TRUE, sep = " ",col.names = TRUE,quote = FALSE)
  }
  return(my.stats)  # return a list with statistics
}


# ------------------------------------------------------------------------------------------------------
# WRAPPER METHODS
# ------------------------------------------------------------------------------------------------------
GLMregW <- function(my.datf,my.datf.train,my.datf.test,fDet=FALSE,outFile="") { # with wrapper
  # to be implemented
  return("")  # return statistics
}
