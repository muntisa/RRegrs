# Libraries and external custom functions
library(caret)                  # add package caret

# -----------------------------------
# (1.2) Load the ORIGINAL DATASET
# -----------------------------------
cat("-> [1] Loading original dataset ...\n")
# (it can contain errors, correlations, near zero variance columns)
ds.dat0 <- read.csv("ds.csv",header=T)                              # original dataset frame

# resolving the text to number errors for future calculations
ds.indx<- colnames(ds.dat0)[2:dim(ds.dat0)[2]]                    # FEATURE names (no dependent variable)
ds.dat1<- ds.dat0[1:dim(ds.dat0)[1],2:dim(ds.dat0)[2]]            # dataset as columns
ds.dat1<- apply(ds.dat1,1,function(x)as.numeric(as.character(x))) # dataset as row vectors to be used with caret!!!

# dependent variable
net.c<- ds.dat0[,1]
net.c<- as.numeric(as.character(net.c)) # values
# full ds frame with training and test
ds<- as.data.frame(cbind(net.c,t(ds.dat1)))

source("s6.DsSplit.R")
PathDataSet = "D:/GitHubs/RRegrs/Rscripts/_Lab"
iSeed=1                 # to reapeat the ds splitting, different values of seed will be used
dsList  <- DsSplit(ds,3/4,TRUE,PathDataSet,1) # return a list with 2 datasets = dsList$train, dsList$test
# get train and test from the resulted list
ds.train<- dsList$train
ds.test <- dsList$test

# VARIABLE NAME CHANGE !!!!!!
my.datf.train <- ds.train 
my.datf.test <- ds.test


# =============================================================================
# INFO from the main script
# =============================================================================
iSplit=1
sCV="repeatedcv"
outFile="test.csv"

# TRAINING MODEL
# --------------------
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

# LEVERAGE
# --------------------

pred.tr     <- predict(pls.fit,my.datf.train) # predicted Y
pred.ts     <- predict(pls.fit,my.datf.test)  # predicted Y

noFeats.fit <- length(predictors(pls.fit))    # no. of features from the fitted model
Feats.fit   <- paste(predictors(pls.fit),collapse="+") # string with the features included in the fitted model
FeatImp <- varImp(pls.fit, scale = F)


# COPY FROM HERE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
pdf(file=paste(outFile,".",sCV,".","split",iSplit,".pdf"))
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