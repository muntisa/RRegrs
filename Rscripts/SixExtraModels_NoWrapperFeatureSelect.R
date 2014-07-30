#
#Methods included are: lm (no cv), pls, lasso, svmLinear, neural network (single hidden layer), rbf.
#
#REsampling: 10fold CV (no repeats)
#
##
#

	lm.fit<- lm(as.matrix(my.y.vec)~as.matrix(my.descr.dat))#no cv
	lm.fit.sum<- summary(lm.fit)
	lm.fit.sum

	set.seed(2)
	pls.fit<- train(my.y.vec~.,data=my.datf.train,
	method = 'pls', tuneLength = 10, trControl = ctrl,
	metric = 'RMSE',tuneGrid=expand.grid(.ncomp=c(1:15)),preProc = c('center', 'scale'))
	#ncomp== the number of components to include in the model
	pls.fit

	set.seed(2)
	las.fit<- train(my.y.vec~.,data=my.datf.train,
	method = 'lasso', tuneLength = 10, trControl = ctrl,
	metric = 'RMSE',tuneGrid=expand.grid(.fraction= seq(0.1,1,by=0.1)),preProc = c('center', 'scale'))
	#fraction= fraction of full solution
	las.fit

	set.seed(2)
	svmL.fit<- train(my.y.vec~.,data=my.datf.train,
	method = 'svmRadial', tuneLength = 10, trControl = ctrl,
	metric = 'RMSE',preProc = c('center', 'scale'))
	tuneGrid=expand.grid(.sigma=seq(0,1,0.1),.C= c(1:10)))
	#C==cost (a regularization parameter)that controls how much the regression
	#line can adapt to the data smaller values result in more linear, i.e.
	#flat surfaces.
	#sigma==the scale parameter for radial bass function (rbf)
	svmL.fit

	set.seed(2)
	nn.fit<- train(my.y.vec~.,data=my.datf.train,
	method = 'nnet',trControl = ctrl,
	linout=TRUE, trace = FALSE,MaxNWts=20000,preProc = c('center', 'scale'),
	#Grid of tuning parameters to try:
	tuneGrid=expand.grid(.size=c(1,5,10,15),.decay=c(0,0.001,0.1)))
	#Grid parameters are appearing at the print out of the model
	#size==#of units in hidden layer, decay==parameter of weight decay (default:0)
	nn.fit
	varImp(nn.fit)


	set.seed(2)
	rbf.fit<- train(my.y.vec~.,data=my.datf.train,
	method = 'rbfDDA',trControl = ctrl,preProc = c('center', 'scale'),
	tuneGrid=expand.grid(.negativeThreshold=seq(0,1,0.1)))
	#size== number of hidden units
	rbf.fit

