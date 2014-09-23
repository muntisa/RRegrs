## svm regression function
## jseoane
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
  
  
 

