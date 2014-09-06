#
#Methods included are: pls, lasso, lm, lmAIC, svmLinear, svmRadial, neural network (single hidden layer),
#som, rbf.
#
#For those which do not include feature selection, rfe() function has been used (stands for recursive feature selection).
#
#REsampling: 10fold CV (no repeats), boostrap for rfe using the default 25 samples
#
#adjsuted R2 and RMSE functions are included at the top.
#

r2.adj.funct<- function(y,y.new,num.pred){#y==y, y.new=predicted, num.pred=number of idependent variables (predictors)
	y.mean<- mean(y)
	x.in<- sum((y-y.new)^2)/sum((y-y.mean)^2)
	x.in<- 1-x.in #r squared

	x.in<- (1-x.in)*((length(y)-1)/(length(y)-num.pred-1))
	x.in<- 1 - x.in 
	return(x.in)
}

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


rmse.funct<- function(y,y.new){
	return(sqrt(mean((y.new - y)^2)))
}


#########


	ctrl<- rfeControl(method = 'boot', number = 25,saveDetails=T)
	ctrl1<- trainControl(method = 'repeatedcv', number = 10,repeats = 1,
	summaryFunction = defaultSummary,savePred=T)

subsets<- seq(5,dim(my.datf.train)[2]-1,by=10)#CHECk

	set.seed(2)
	pls.fit<- rfe(my.y.vec~.,data=my.datf.train, 
	method = 'pls',
	rfeControl = ctrl, trControl=ctrl1, sizes=subsets, importance=T,
	metric = 'RMSE', preProc = c('center', 'scale'),
	tuneGrid=expand.grid(.ncomp=c(1:5)))#15)))
	#,tuneLength = 10, sizes=subsets,
	#ncomp== the number of components to include in the model
	pls.fit ; predictors(pls.fit)
	pls.fit$resample
	pls.fit$fit
	varImp(pls.fit)#?varImp: what is the importance statistic

	set.seed(2)
	las.fit<- train(my.y.vec~.,data=my.datf.train,
	method = 'lasso', tuneLength = 10, trControl = ctrl1,
	metric = 'RMSE',tuneGrid=expand.grid(.fraction= seq(0.1,1,by=0.1)),
	preProc = c('center', 'scale'))
	las.fit ; predictors(las.fit)

	set.seed(2)
	lm.fit<- train(my.y.vec~.,data=my.datf.train,
	method = 'lm', tuneLength = 10, trControl = ctrl1, importance=T,
	metric = 'RMSE', preProc = c('center', 'scale'))
	lm.fit ; predictors(lm.fit)

	set.seed(2)
	lmAIC.fit<- train(my.y.vec~.,data=my.datf.train,
	method = 'glmStepAIC', tuneLength = 10, trControl = ctrl1,importance=T,
	metric = 'RMSE', preProc = c('center', 'scale'))
	lmAIC.fit ; predictors(lmAIC.fit)

	set.seed(2)
	svmL.fit<- rfe(my.y.vec~.,data=my.datf.train,
	method = 'svmLinear', tuneLength = 10,sizes=subsets,importance=T,
	rfeControl = ctrl, trControl=ctrl1,
	metric = 'RMSE',preProc = c('center', 'scale'),
	tuneGrid=expand.grid(.C= c(1:10)))
	svmL.fit ; predictors(svmL.fit)

	set.seed(2)
	svmR.fit<- rfe(my.y.vec~.,data=my.datf.train,
	method = 'svmRadial', tuneLength = 10,sizes=subsets,importance=T,
	rfeControl = ctrl, trControl=ctrl1,
	metric = 'RMSE',preProc = c('center', 'scale'),
	tuneGrid=expand.grid(.sigma=seq(0,1,0.1),.C= c(1:10)))
	svmR.fit ; predictors(svmR.fit)

	set.seed(2)
	nn.fit<- rfe(my.y.vec~.,data=my.datf.train,
	method = 'nnet',sizes=subsets,importance=T,
	rfeControl = ctrl,trControl=ctrl1,
	linout=TRUE, trace = TRUE,MaxNWts=20000,preProc = c('center', 'scale'),
	tuneGrid=expand.grid(.size=c(1,5,10,15),.decay=c(0,0.001,0.1)))
	#Grid parameters are appearing at the print out of the model
	#size==#of units in hidden layer, decay==parameter of weight decay (default:0)
	nn.fit ; predictors(nn.fit)

	set.seed(2)
	som.fit<- rfe(my.y.vec~.,data=my.datf.train,
	method = 'xyf',sizes=subsets,
	rfeControl = ctrl, trControl=ctrl1,preProc = c('center', 'scale'),
	tuneGrid=expand.grid(.xdim=c(1:5),.ydim=c(2:5),
				.xweight=c(0.2,0.5,0.7,0.9),.topo=c(1,2)))
	#Grid parameters are appearing at the print out of the model
	#size==#of units in hidden layer, decay==parameter of weight decay (default:0)
	som.fit ; predictors(som.fit)


	set.seed(2)
	rbf.fit<- train(my.y.vec~., data=my.datf.train,
	method = 'rbfDDA',trControl = ctrl1, preProc = c('center', 'scale'),
	tuneGrid=expand.grid(.negativeThreshold=seq(0,1,0.1)),importance=T)
	#nprune==maximum number of terms in the prune model (i.e. number of descriptors/proteins)
	#degree==maximum degree of interactions, 1 stands for additive model
	rbf.fit ; predictors(rbf.fit)



 pls.fit$results[which.min(pls.fit$results[,2]),]
 bT<- pls.fit$fit$bestTune[1,1]
 predAll<- pls.fit$fit$pred
 predAll.w<- which(predAll[,4]==bT)
 predSel<- predAll[predAll.w,1:2]# 1st column the y_pred, 2nd column the y_obs
 r2.adj.lm.funct(predSel[,2],predSel[,1],pls.fit$bestSubset)
 pls.fit$optVariables # selected variables
 pls.test.res<- postResample(predict(pls.fit,my.datf.test),my.datf.test[,1])


