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
#----------------------------------------------------------------------------------------------------------------------
rmse.funct<- function(y,y.new){
  return(sqrt(mean((y.new - y)^2)))
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
  
  # CV types: 10-CV and LOOCV
  CVtypes = c("repeatedcv","LOOCV")
  
  # for each type of CV do all the statistics (get 1 list for each CV and merge all)
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
    RMSE = lm.fit$results[,2]
    Rsquared = lm.fit$results[,3]
    if (cv == 1){ # if 10-fold CV
      RMSE_SD = lm.fit$results[,4]
      Rsquared_SD = lm.fit$results[,5]
    }
    
    #------------------------------------------------
    # RMSE & R^2, for train/ test respectively
    #------------------------------------------------
    lm.train.res<- getTrainPerf(lm.fit)
    lm.test.res <- postResample(predict(lm.fit,my.datf.test),my.datf.test[,1])
    
    #------------------------------------------------
    # Adj R2: y obs, y predicted, No of predictors
    #------------------------------------------------
    ds.full = rbind(my.datf.train,my.datf.test)
    adjR2_both = r2.adj.funct(ds.full[,1], predict(lm.fit,ds.full), length(predictors(lm.fit)))
    adjR2_train= r2.adj.funct(my.datf.train[,1],predict(lm.fit,my.datf.train),length(predictors(lm.fit)))
    adjR2_test = r2.adj.funct(my.datf.test[,1], predict(lm.fit,my.datf.test), length(predictors(lm.fit)))
    
    # generate a list with statistics for each cross-validation type
    if (cv == 1){  # output for 10-CV
      my.stats.10CV = list("RegrMethod"="GLM.AIC","SplitNo"=1,"RMSE.tr.10CV"= RMSE,"R2.tr.10CV" = Rsquared,
                           "RMSEsd.tr.10CV" = RMSE_SD,"R2sd.tr.10CV" = Rsquared_SD,
                           "RMSE.ts.10CV" = (lm.test.res["RMSE"][[1]]),"R2.ts.10CV" = (lm.test.res["Rsquared"][[1]]),
                           "adjR2.both.10CV" = adjR2_both, "adjR2.tr.10CV" = adjR2_train, "adjR2.ts.10CV" = adjR2_test)
      # default values: "Regrrs"="GLM","Step"=1 (or for one single use of funtion)
    }
    if (cv == 2){  # output for LOOCV
      my.stats.LOOCV = list("RMSE.tr.LOOCV"= RMSE,"R2.tr.LOOCV" = Rsquared,
                            "RMSE.ts.LOOCV" = (lm.test.res["RMSE"][[1]]),"R2.ts.LOOCV" = (lm.test.res["Rsquared"][[1]]),
                            "adjR2.both.LOOCV" = adjR2_both, "adjR2.tr.LOOCV" = adjR2_train, "adjR2.ts.LOOCV" = adjR2_test)
    }
    
    # TO ADD !!!!!!!!!! here or in the main script !!!!!!!!!!!!
    # dsFileName, DateTime, CaseNo, FeatureNo, FeatureList, OutVar
    # RMSE.both10CV, RMSE.bothLOOCV, R2.both10CV, R2.bothLOOCV
    # 
    
    
    #--------------------------------------
    # Write to file DETAILS for GLM
    #--------------------------------------
    if (fDet==TRUE) {  
      sink(outFile) # to be modified to append results!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      print("Training Set: ");print(summary(my.datf.train))
      print("Test Set: ");print(summary(my.datf.test))
      print("Fitting: ");print(lm.fit)
      print("Predictions: ");print(predictors(lm.fit))
      print("Trainig Results: ");print(lm.train.res)
      print("Test Results: ");print(lm.test.res)
      if (cv == 1){
        print("10-fold CV statistics: ");print(data.frame(my.stats.10CV))
      }
      if (cv == 2){
        print("LOOCV statistics: ");print(data.frame(my.stats.LOOCV))
      }
      sink()
      # file.show(outFile)
    }
    #--------------------------------------
    
  } # END cross-validations
  my.stats.full = c(my.stats.10CV,my.stats.LOOCV)   # merge the CV results into one list
  return(my.stats.full)  # return a list with statistics
}

# ------------------------------------------------------------------------------------------------------
GLMregW <- function(my.datf,my.datf.train,my.datf.test,fDet=FALSE,outFile="") { # with wrapper
  # to be implemented
  return("")  # return statistics
}
