# =============================================================================
# INFO from the main script
# =============================================================================
iSplit=1
sCV="repeatedcv"
net.c = ds.train[,1]   # make available the names of variables from training dataset
RegrMethod <- "lm" # type of regression

# Define the CV conditions
ctrl<- trainControl(method = sCV, number = 10,repeats = 10,
                    summaryFunction = defaultSummary)

# Train the model using only training set
set.seed(iSplit)

pls.fit<- train(net.c~.,data=ds.train,
                method = 'pls', tuneLength = 10, trControl = ctrl,
                metric = 'RMSE',
                tuneGrid=expand.grid(.ncomp=c(1:(dim(ds.train)[2]-1))))

# --------------------

fitModel <- pls.fit$finalModel

# =============================================================================
# Assessment of Applicability Domain (plot leverage)
# =============================================================================

# Fitted vs Residuals
plot(fitted(fitModel),residuals(fitModel),
     main="Fitted vs. Residuals for Fitted Model",
     xlab="Fitted", ylab="Residuals")
abline(h = 0, lty = 2)

resids <- residuals(fitModel) # residuals

# ??
#ds.train
# hat.tr <- hat(ds.train[,2:dim(ds.train)[2]])#gives you the hat values!

#IF model is lm, glm, then the following commend dies the same thing
#(again i'm using notation from the original code i've uploaded)

hat.fit <- hatvalues(fitModel)          # hat values
hat.fit.df <- as.data.frame(hat.fit)    # hat data frame
hat.mean <- mean(hat.fit)               # mean hat values
hat.fit.df$warn <- ifelse(hat.fit.df[, 'hat.fit']>3*hat.mean, 'x3',ifelse(hat.fit.df[, 'hat.fit']>2*hat.mean, 'x2', '-' ))

#THRESHOLD values: 3m/n
# where m is the number of parameters, and n number of observations

thresh.lever<- (3*(dim(ds.train)[2]-1))/dim(ds.train)[1] # leverage thresh
hat.problems<- data.frame(hat.fit[hat.fit>thresh.lever]) # points with high leverage

# Leverage plots
plot(hat.fit, type = "h",
     main="Leverage for Fitted Model",
     xlab="Index", ylab="Hat")
abline(h = thresh.lever, lty = 2, col="red") # leverage thresh

#OTHER things that can be calculated are, e.g.
cook.dists<- cooks.distance(fitModel)

# Cook distances plot
plot(cook.dists,
     main="Cook's Distance for Fitted Model",
     xlab="Index", ylab="Cook Distance")

cutoff.Cook <- 4/((nrow(ds.train)-length(fitModel$coefficients)-2)) 
for (p in 1:6) {
  plot(fitModel, which=p, cook.levels=cutoff.Cook)
}

infl <- influence(fitModel)#produces several statistics of the kind

#PACKAGE 'ouliers' gives you options for outliers tests which is pretty good
