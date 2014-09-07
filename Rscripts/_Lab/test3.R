# VARIABLE NAME CHANGE !!!!!!
my.datf.train <- ds.train 
my.datf.test <- ds.test

# =============================================================================
# INFO from the main script
# =============================================================================
iSplit=1
sCV="repeatedcv"
outFile="test.csv"

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
# --------------------
pred.tr     <- predict(lm.fit,my.datf.train) # predicted Y for training
pred.ts     <- predict(lm.fit,my.datf.test)  # predicted Y for test
noFeats.fit <- length(predictors(lm.fit))    # no. of features from the fitted model
Feats.fit   <- paste(predictors(lm.fit),collapse="+") # string with the features included in the fitted model
FeatImp <- varImp(lm.fit, scale = F)

# COPY FROM HERE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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

# PDF plots
# --------------------------------------------------------------
pdf(file=paste(outFile,".",sCV,".","split",iSplit,".pdf"))
# par(mfrow = c(3, 4)) # all plots into one page!

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

# Cook's distance
plot(cook.dists,
     main="Cook's Distance for Fitted Model",
     xlab="Index", ylab="Cook Distance")

for (p in 1:6) {
  plot(fitModel, which=p, cook.levels=cutoff.Cook)
}

# plot(FeatImp, top = components,main="Feature Importance") # ERROR !
dev.off()
# --------------------------------------------------------------