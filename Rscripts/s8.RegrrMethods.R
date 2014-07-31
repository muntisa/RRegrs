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
  my.stats<- list() # create empty result
  
  # specify CV parameters = 10-fold cross validation
  ctrl<- trainControl(method = 'repeatedcv', number = 10,repeats = 10,
                      summaryFunction = defaultSummary)
  
  # Training the model using only training set
  set.seed(2)
  attach(my.datf.train)
  lm.fit<- train(net.c~.,data=my.datf.train,
                 method = 'glmStepAIC', tuneLength = 10, trControl = ctrl,
                 metric = 'RMSE')
  
  #------------------------------
  # Training RESULTS
  #------------------------------
  RMSE = lm.fit$results[,2]
  Rsquared = lm.fit$results[,3]
  RMSE_SD = lm.fit$results[,4]
  Rsquared_SD = lm.fit$results[,5]
  
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
  
  my.stats = list("Regrrs"="GLM","Step"=1,"tr.RMSE.10CV"= RMSE,"tr.R2.10CV" = Rsquared,
                  "tr.RMSESD.10CV" = RMSE_SD,"tr.R2SD.10CV" = Rsquared_SD,
                  "ts.RMSE.10CV" = (lm.test.res["RMSE"][[1]]),"ts.R2.10CV" = (lm.test.res["Rsquared"][[1]]),
                  "both.adjR2.10CV" = adjR2_both, "tr.adjR2.10CV" = adjR2_train, "ts.adjR2.10CV" = adjR2_test)
  # default values: "Regrrs"="GLM","Step"=1 (or for one single use of funtion)
  
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
    print("Final Statistics: ");print(data.frame(my.stats))
    sink()
    # file.show(outFile)
  }
  return(my.stats)  # return a list with statistics
}

# ------------------------------------------------------------------------------------------------------
GLMregW <- function(my.datf,my.datf.train,my.datf.test,fDet=FALSE,outFile="") { # with wrapper
  # to be implemented
  return("")  # return statistics
}
