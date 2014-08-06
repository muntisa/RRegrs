# ======================================================================
# Regression function based on caret package
# ----------------------------------------------------------------------
# contact: Cristian R Munteanu | BiGCaT - UM    | muntisa@gmail.com
#          Georgia Tsiliki     | ChemEng - NTUA | g_tsiliki@hotmail.com

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

#====================================================================================================
# 8.2- GLM stepwise - based on AIC (caret); Generalized Linear Model with Stepwise Feature Selection
#====================================================================================================
# Inputs:
# - my.datf.train,my.datf.test = training and test dataset frames
# - fDet = flag for detais (TRUE/FALSE)
# - outFile = output file for GLM details
# Output:
# - list of statistics (Regrrs, Step, tr.RMSE.10CV, tr.R2.10CV, tr.RMSESD.10CV,
#   tr.R2SD.10CV, ts.RMSE, ts.R2, both.adjR2, tr.adjR2, ts.adjR2
#   (tr = train, ts = test, both = tr+ts = full dataset, 10Cv = 10-fold cross-validation)
# ---------------------------------------------------------------------------------------------------
GLMreg <- function(my.datf.train,my.datf.test,fDet=FALSE,outFile="") {
  attach(my.datf.train)   # make available the names of variables from training dataset
  
  CVtypes    <- c("repeatedcv","LOOCV")             # CV types: 10-CV and LOOCV
  RegrMethod <- "GLM.AIC"                           # type of regression
  nCases     <- dim(my.datf.train)[1]               # number of cases
  nFeatures  <- dim(my.datf.train)[2] - 1           # number of input features (it could be different with the fitted model features!)
  FeatureList<- paste(names(my.datf.train)[2:nFeatures], collapse="+") # list of input feature names, from second name, first name is the output variable
  OutVar     <- names(my.datf.train)[1]             # name of predicted variable (dependent variable)
  
  # Data set info
  my.stats.dsInfo <- list("RegrMethod"= RegrMethod,           # regression method name
                          "NoCases"= as.numeric(nCases),      # no. of cases
                          "InNoVars"= as.numeric(nFeatures),  # no. of input features
                          "InFeatures"= FeatureList,          # list with feature names
                          "PredVar" = OutVar,                 # name of predicted variable (dependent variable)
                          "SplitNo"= 1)                       # index for the splitting; default is 1 but it will be modified outside the function for each split
  
  # For each type of CV do all the statistics (get 1 list for each CV and merge all)
  for (cv in 1:length(CVtypes)) {
    
    ctrl<- trainControl(method = CVtypes[cv], number = 10,repeats = 10,
                        summaryFunction = defaultSummary)
    
    # Training the model using only training set
    set.seed(cv)
    lm.fit<- train(net.c~.,data=my.datf.train,
                   method = 'glmStepAIC', tuneLength = 10, trControl = ctrl,
                   metric = 'RMSE')
    
    #------------------------------
    # Training RESULTS
    #------------------------------
    RMSE.tr  <- lm.fit$results[,2]
    R2.tr    <- lm.fit$results[,3]
    if (cv == 1){ # if 10-fold CV
      RMSEsd.tr <- lm.fit$results[,4]
      R2sd.tr   <- lm.fit$results[,5]
    }
    if (cv == 2){ # if 10-fold CV
      RMSEsd.tr <- 0 # formulas will be added later!
      R2sd.tr   <- 0 # formulas will be added later!
    }
    
    #------------------------------------------------
    # RMSE & R^2, for train/ test respectively
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
    
    # Generate a list with statistics for each cross-validation type
    # --------------------------------------------------------------------
    # There are 3 lists for Data set info, 10-fold CV and LOOCV
    # These 3 lists will be merged into one in order obtain the function output including the header (statistics names)
    
    if (cv == 1){  # output for 10-fold CV (it shold be modified with unique bloque for each cv, dynamic labels of the fields, etc.)
      my.stats.10CV <- list("NoModelFeats.10CV"  = as.numeric(noFeats.fit),
                            "ModelFeats.10CV"    = Feats.fit,
                            "adjR2.tr.10CV"  = as.numeric(adjR2.tr),
                            "RMSE.tr.10CV"   = as.numeric(RMSE.tr),
                            "R2.tr.10CV"     = as.numeric(R2.tr),
                            "RMSEsd.tr.10CV" = as.numeric(RMSEsd.tr),
                            "R2sd.tr.10CV"   = as.numeric(R2sd.tr),
                            "adjR2.ts.10CV"= as.numeric(adjR2.ts),
                            "RMSE.ts.10CV" = as.numeric((lm.test.res["RMSE"][[1]])),
                            "R2.ts.10CV"   = as.numeric((lm.test.res["Rsquared"][[1]])),
                            "corP.ts.10CV" = as.numeric(corP.ts),
                            "adjR2.both.10CV" = as.numeric(adjR2.both),
                            "RMSE.both.10CV"  = as.numeric(RMSE.both),
                            "R2.both.10CV"    = as.numeric(r2.both))
      
      # default values: "Regrrs"="GLM","Step"=1 (or for one single use of funtion)
    }
    if (cv == 2){  # output for LOOCV
      my.stats.LOOCV <- list("NoModelFeats.10CV"= as.numeric(noFeats.fit),
                             "ModelFeats.10CV"  = Feats.fit,
                             "adjR2.tr.LOOCV"   = as.numeric(adjR2.tr),
                             "RMSE.tr.LOOCV"    = as.numeric(RMSE.tr),
                             "R2.tr.LOOCV"      = as.numeric(R2.tr),
                             "RMSEsd.tr.LOOCV"  = as.numeric(RMSEsd.tr),
                             "R2sd.tr.LOOCV"    = as.numeric(R2sd.tr),
                             "adjR2.ts.LOOCV" = as.numeric(adjR2.ts),
                             "RMSE.ts.LOOCV"  = as.numeric((lm.test.res["RMSE"][[1]])),
                             "R2.ts.LOOCV"    = as.numeric((lm.test.res["Rsquared"][[1]])),
                             "corP.ts.LOOCV"  = as.numeric(corP.ts),
                             "adjR2.both.LOOCV"  = as.numeric(adjR2.both),
                             "RMSE.both.LOOCV"   = as.numeric(RMSE.both),
                             "R2.both.LOOCV"     = as.numeric(r2.both))
    } # default values: "Regrrs"="GLM","Step"=1 (or for one single use of funtion)
    
    # TO ADD !!!!!!!!!! here or in the main script !!!!!!!!!!!!
    # dsFileName, DateTime
    
    #---------------------------------------------------------------------
    # Write to file DETAILS for GLM for each cross-validation method
    #---------------------------------------------------------------------
    if (fDet==TRUE) {   # if flag for details if true, print details about any resut
      write.table("Fitting Summary: ", file = outFile,append = TRUE, sep = " ",col.names = FALSE,quote = FALSE)
      write.table(data.frame(summary(lm.fit)$coefficients), file = outFile,append = TRUE, sep = " ",col.names = TRUE,quote = FALSE)
      
      write.table("Predictors: ", file = outFile,append = TRUE, sep = " ",col.names = FALSE,quote = FALSE)
      write.table(predictors(lm.fit), file = outFile,append = TRUE, sep = " ",col.names = TRUE,quote = FALSE)
      
      write.table("Trainig Results: ", file = outFile,append = TRUE, sep = " ",col.names = FALSE,quote = FALSE)
      write.table(predictors(lm.train.res), file = outFile,append = TRUE, sep = " ",col.names = TRUE,quote = FALSE)
      write.table("Test Results: ", file = outFile,append = TRUE, sep = " ",col.names = FALSE,quote = FALSE)
      write.table(predictors(lm.test.res), file = outFile,append = TRUE, sep = " ",col.names = TRUE,quote = FALSE)
    }
    #--------------------------------------
  } # END cross-validations
  my.stats.full <- c(my.stats.dsInfo,my.stats.10CV,my.stats.LOOCV)   # merge the CV results into one list that contains the names of each field!
  
  #--------------------------------------
  # Write to file DETAILS for GLM
  #--------------------------------------
  if (fDet==TRUE) {   # if flag for details if true, print details about any resut
  
    write.table("Full Statistics: ", file = outFile,append = TRUE, sep = " ",col.names = FALSE,quote = FALSE)
    write.table(my.stats.full, file = outFile,append = TRUE, sep = " ",col.names = TRUE,quote = FALSE)
  }
  
  return(my.stats.full)  # return a list with statistics
}

# ------------------------------------------------------------------------------------------------------
GLMregW <- function(my.datf,my.datf.train,my.datf.test,fDet=FALSE,outFile="") { # with wrapper
  # to be implemented
  return("")  # return statistics
}
