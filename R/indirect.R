
#' @export
design.lasso <- function(train){
  mat = cbind(train$A,train$X)
  colnm=c("dose",paste("X",seq(1,ncol(train$X)),sep=""))
  colnames(mat) =colnm
  sqMat=as.matrix(mat[,1])^2
  colnames(sqMat)='dose:dose'
  newcolnm=paste0(c(paste0(colnm,collapse = '+'),paste0('dose:',colnm[-1])),collapse='+')
  newcolnm=paste0('~0+',newcolnm)
  design_mat=model.matrix(as.formula(newcolnm),data=data.frame(mat))
  design_mat=cbind(design_mat,sqMat)
  return(design_mat)
}
#' @export
cv.lasso<-function(train,test,nfolds,two.sided=FALSE,lower=TRUE,Lambda=c(2^(-(0:10)),0.00000001),...){

  n=length(train$A)
  p=dim(train$X)[2]
  nsize=n/nfolds
  misclass=0

  design.train=design.lasso(train)
  design.test=design.lasso(test)
  pf=rep(1,dim(design.train)[2])
  pf[which(colnames(design.train)%in%c('dose:dose'))]=0
  cvglmfit = cv.glmnet(design.train,train$Y,nfolds=nfolds,penalty.factor=pf)

  pred=predict(cvglmfit,design.test)>train$S
  misclass=sum(test$alpha[,1]*(!pred)*(test$Y>test$S)+test$alpha[,2]*(pred)*(test$Y<test$S))/length(test$Y)
  cost=sum(test$alpha[,1]*(pred)*pmax((test$S-test$Y),0)+test$alpha[,2]*(!pred)*pmax((test$Y-test$S),0))/length(test$Y)

  search=searchInterval(fitted=cvglmfit,test=test,two.sided=two.sided,classification = FALSE,lower=lower)
  measures=measurement(search,test,two.sided,family,lower)
  return(list(fitted=cvglmfit,misclass=misclass,cost=cost,bestPara=cvglmfit$lambda.1se,search=search,pred=pred,measures=measures))
}
#' @export
cv.SVR<-function(train,test,nfolds,two.sided=FALSE,Cost=1.3^((-10):10),lower=TRUE,...){
  n=length(train$A)
  p=dim(train$X)[2]
  nsize=n/nfolds
  misclass<-NA->cost
  prederror=Inf
  trainX = cbind(train$A,train$X)
  testX = cbind(test$A,test$X)

  svmfit=e1071::best.svm(x=trainX,y=train$Y,type='nu-regression',tunecontrol = tune.control(cross=nfolds),cost=Cost)
  bestPara=svmfit$cost
  
  pred=predict(svmfit,newdata=testX)
  pred_test=pred>test$S

  pred=predict(svmfit,newdata=testX,type='response')>test$S
  misclass=(sum(test$alpha[,1]*(pred)*(test$Y<test$S))+sum(test$alpha[,2]*(!pred)*(test$Y>test$S)))/length(test$A)
  cost=(sum(test$alpha[,1]*(pred)*pmax((test$S-test$Y),0))+sum(test$alpha[,2]*(!pred)*pmax((test$Y-test$S),0)))/length(test$A)

  search=searchInterval(fitted=svmfit,test=test,two.sided=two.sided,classification = FALSE,lower=lower)
  measures=measurement(search,test,two.sided,family,lower)
  return(list(bestPara=bestPara,misclass=misclass,cost=cost,search=search,pred=pred,measures=measures))
}
#' @export
cv.foreg<-function(train,test,nfolds,two.sided=FALSE,lower=TRUE,Mtry=c(0.3,0.4,0.5),...){
  n=length(train$A)
  p=dim(train$X)[2]
  nsize=n/nfolds
  Perform=cbind(matrix(0,ncol=1,nrow=length(Mtry)),Mtry)
  Pred_test=matrix(0,ncol=n,nrow=length(Mtry))
  misclass=0
  cost=0
  for(i in 1:nfolds){
    cat(paste("=== Fold ",i,"===",sep=" "))

    index = seq(((i-1)*nsize+1),(i*nsize))
    train0 = list(X=train$X[-index,],A=train$A[-index],Y=train$Y[-index],opt=train$opt[-index],mu=train$mu[-index],S=train$S,propensity=train$propensity[-index])
    test0 = list(X=train$X[index,],A=train$A[index],Y=train$Y[index],opt=train$opt[index],mu=train$mu[index],S=train$S,propensity=train$propensity[index])
    for(k in 1:length(Mtry)){
      Dat=cbind(train0$A,train0$X)
      colnames(Dat)=NULL
      fit.forest=randomForest(Dat,train0$Y,mtry=round(Mtry[k]*p+1))
      Dat=cbind(test0$A,test0$X)
      colnames(Dat)=NULL
      pred=ifelse(predict(fit.forest,Dat)>train$S,TRUE,FALSE)
      Perform[k,1]=  Perform[k,1]+sum(test0$alpha[,2]*(!pred)*(test0$Y>test0$S)+test0$alpha[,1]*(pred)*(test0$Y<test0$S))/length(test0$Y)
      Pred_test[k,index]=pred
    }
  }
  Perform[,1]=Perform[,1]/nfolds
  id=which.min(Perform[,1])
  names(Perform)=c("misclass","Mtry")
  print(Perform[id,])

  Dat=cbind(train$A,train$X)
  colnames(Dat)=NULL
  fit.forest=randomForest(Dat,train$Y,mtry=round(Perform[id,2]*p+1))
  Dat=cbind(test$A,test$X)
  pred=ifelse(predict(fit.forest,Dat)>train$S,TRUE,FALSE)
  misclass=sum(test$alpha[,2]*(!pred)*(test$Y>test$S)+test$alpha[,1]*(pred)*(test$Y<test$S))/length(test$Y)
  cost=sum(test$alpha[,1]*(pred)*pmax((test$S-test$Y),0)+test$alpha[,2]*(!pred)*pmax((test$Y-test$S),0))/length(test$Y)
  percentage=sum(pred*(test$Y>test$S))/sum(pred)

  search=searchInterval(fitted=fit.forest,test=test,two.sided=two.sided,lower=lower,classification = FALSE)
  measures=measurement(search,test,two.sided,family,lower)
  return(list(bestPara=Perform[id,2],misclass=misclass,cost=cost,search=search,pred=pred,percentage=percentage,measures=measures))
}
#' @export
cv.logist<-function(train,test,nfolds,two.sided=FALSE,lower=TRUE,Lambda=c(2^(-(0:10)),0.00000001),...){

  n=length(train$A)
  p=dim(train$X)[2]
  nsize=n/nfolds
  misclass=0

  design.train=design.lasso(train)
  design.test=design.lasso(test)
  pf=rep(1,dim(design.train)[2])
  pf[which(colnames(design.train)%in%c('dose:dose'))]=0
  cvglmfit = cv.glmnet(design.train,train$Y>train$S, family="binomial",nfolds=nfolds,penalty.factor=pf)

  pred=predict(cvglmfit,design.test)>0
  misclass=sum(test$alpha[,1]*(!pred)*(test$Y>test$S)+test$alpha[,2]*(pred)*(test$Y<test$S))/length(test$Y)
  cost=sum(test$alpha[,1]*(pred)*pmax((test$S-test$Y),0)+test$alpha[,2]*(!pred)*pmax((test$Y-test$S),0))/length(test$Y)
  percentage=sum(pred*(test$Y>test$S))/sum(pred)

  search=searchInterval(fitted=cvglmfit,test=test,two.sided=two.sided,lower=lower)
  measures=measurement(search,test,two.sided,family,lower)
  return(list(fitted=cvglmfit,bestPara=cvglmfit$lambda.1se,misclass=misclass,cost=cost,search=search,pred=pred,percentage=percentage,measures=measures))
}

#' @export
cv.SVM<-function(train,test,nfolds,two.sided=FALSE,Cost=1.3^((-10):10),lower=TRUE,...){
  n=length(train$A)
  p=dim(train$X)[2]
  nsize=n/nfolds
  misclass<-NA->cost
  prederror=Inf
  trainX = cbind(train$A,train$X)
  testX = cbind(test$A,test$X)

  svmfit=e1071::best.svm(x=trainX,y=(train$Y>train$S),probability = TRUE,type='C-classification',tunecontrol = tune.control(cross=nfolds),cost=Cost)
  bestPara=svmfit$cost
  
  pred=ifelse(predict(svmfit,newdata=testX,type='response')=='TRUE',1,0)
  misclass=(sum(test$alpha[,1]*(pred)*(test$Y<test$S))+sum(test$alpha[,2]*(!pred)*(test$Y>test$S)))/length(test$A)
  cost=(sum(test$alpha[,1]*(pred)*pmax((test$S-test$Y),0))+sum(test$alpha[,2]*(!pred)*pmax((test$Y-test$S),0)))/length(test$A)
  percentage=sum(pred*(test$Y>test$S))/sum(pred)

  search=searchInterval(fitted=svmfit,test=test,two.sided=two.sided,lower=lower)
  measures=measurement(search,test,two.sided,family,lower)
  return(list(bestPara=bestPara,misclass=misclass,cost=cost,search=search,pred=pred,percentage=percentage,measures=measures))
}


#' @export
cv.forest<-function(train,test,nfolds,two.sided=FALSE,lower=TRUE,Mtry=c(0.2,0.3,0.4),...){

  n=length(train$A)
  p=dim(train$X)[2]
  nsize=n/nfolds
  Perform=cbind(matrix(0,ncol=1,nrow=length(Mtry)),Mtry)
  Pred_test=matrix(0,ncol=n,nrow=length(Mtry))
  misclass=0
  cost=0
  for(i in 1:nfolds){
    cat(paste("=== Fold ",i,"===",sep=" "))

    index = seq(((i-1)*nsize+1),(i*nsize))
    train0 = subsetData(train,-index)
    test0 = subsetData(train,index)
    for(k in 1:length(Mtry)){
      Dat=cbind(train0$A,train0$X)
      colnames(Dat)=NULL
      fit.forest=randomForest(Dat,as.factor(train0$Y>train0$S),mtry=round(Mtry[k]*p+1))
      Dat=cbind(test0$A,test0$X)
      colnames(Dat)=NULL
      pred=ifelse(predict(fit.forest,Dat)==TRUE,TRUE,FALSE)
      Perform[k,1]=  Perform[k,1]+sum(test0$alpha[,2]*(!pred)*(test0$Y>test0$S)+test0$alpha[,1]*(pred)*(test0$Y<test0$S))/length(test0$Y)
      Pred_test[k,index]=pred
    }
  }
  Perform[,1]=Perform[,1]/nfolds
  id=which.min(Perform[,1])
  names(Perform)=c("misclass","Mtry")
  print(Perform[id,])

  Dat=cbind(train$A,train$X)
  colnames(Dat)=NULL
  fit.forest=randomForest(Dat,as.factor(train$Y>train$S),mtry=round(Perform[id,2]*p+1))
  Dat=cbind(test$A,test$X)
  colnames(Dat)=NULL
  pred=ifelse(predict(fit.forest,Dat)==TRUE,TRUE,FALSE)
  misclass=(sum(test$alpha[,1]*(pred)*(test$Y<test$S))+sum(test$alpha[,2]*(!pred)*(test$Y>test$S)))/length(test$A)
  cost=(sum(test$alpha[,1]*(pred)*pmax((test$S-test$Y),0))+sum(test$alpha[,2]*(!pred)*pmax((test$Y-test$S),0)))/length(test$A)
  percentage=sum(pred*(test$Y>test$S))/sum(pred)

  search=searchInterval(fitted=fit.forest,test=test,two.sided=two.sided,lower=lower)
  measures=measurement(search,test,two.sided,family,lower)
  table(search$found)/length(test$A)
  return(list(bestPara=Perform[id,2],misclass=misclass,cost=cost,search=search,pred=pred,percentage=percentage,measures=measures))
}
