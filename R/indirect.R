fit.lasso.model<-function(train,test,alpha=c(0.5,0.5),nfolds,side){
  n=length(train$A)
  length(test$A)=length(test$A)
  design.mat=design.lasso(train)
  pf=rep(1,dim(design.mat)[2])
  pf[which(colnames(design.mat)%in%c('dose:dose'))]=0
  cvglmfit = cv.glmnet(design.mat,train$Y,nfolds=nfolds,penalty.factor=pf)
  coefs = coef(cvglmfit,s=cvglmfit$lambda.min)

  testX = design.lasso(test)
  pred=predict(cvglmfit,newx=testX)

  index_pos=pred>test$S
  index_neg=pred<test$S

  cost1=sum(alpha[1]*index_pos*pmax((test$S-test$Y),0)+alpha[2]*index_neg*pmax((test$Y-test$S),0))/length(test$A)
  misclass=(alpha[1]*sum(index_pos*(test$Y<test$S))+alpha[2]*sum(index_neg*(test$Y>test$S)))/length(test$A)
  return(list(misclass=misclass,cost=cost1,pred_test=index_pos,fitted=cvglmfit))
}


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

fit.SVR.model<-function(train,test,alpha,side,nfolds,lambda=1.3^(-(30:0))){
  n=length(train$A)
  length(test$A)=length(test$A)
  prederror = 100
  trainX = cbind(train$A,train$X)
  expn=c(0:50)

  for(k in 1:length(expn)){
    #svm.model  = svm(x = train$X[nonzero,], y = (train$A[nonzero]), w= weights[nonzero], type="eps-regression",kernel = "radial",epsilon = eps,cost=-log(lambda),scale=FALSE)
    tmpmodel = ksvm(trainX,train$Y,cross = nfolds,C=1.3^expn[k],kernel = "rbfdot",kpar="automatic")
    print(paste(expn[k],tmpmodel@error))
    if(tmpmodel@error < prederror){
      svmfit = tmpmodel
      prederror = svmfit@error
    }
  }
  testX=cbind(test$A,test$X)
  pred=predict(svmfit,newdata=testX)
  pred_test=pred>test$S
  misclass=(alpha[1]*sum((pred_test)*(test$Y<test$S))+alpha[2]*sum((!pred_test)*(test$Y>test$S)))/length(test$A)
  cost=(alpha[1]*sum((pred_test)*pmax((test$S-test$Y),0))+alpha[2]*sum((!pred_test)*pmax((test$Y-test$S),0)))/length(test$A)
  return(list(misclass=misclass,cost=cost,pred_test=pred_test,fitted=svmfit))
}


cv.lasso<-function(train,test,nfolds=5,alpha=c(0.5,0.5),two.sided=FALSE,
                   aL=NULL,aU=NULL,...){

  if (is.null(aL))  aL<-min(train$A)-0.5*sd(train$A)
  if (is.null(aU))  aU<-max(train$A)+0.5*sd(train$A)
  n=length(train$A)
  p=dim(train$X)[2]
  nsize=n/nfolds
  Pred_test=rep(0,n)
  misclass=0
  cost=0

  design.mat=design.lasso(train)
  pf=rep(1,dim(design.mat)[2])
  pf[which(colnames(design.mat)%in%c('dose:dose'))]=0
  cvglmfit = cv.glmnet(design.mat,train$Y,nfolds=nfolds,penalty.factor=pf)
  coefs = coef(cvglmfit,s=cvglmfit$lambda.min)

  testX = design.lasso(test)
  pred=predict(cvglmfit,newx=testX)

  index_pos=pred>test$S
  index_neg=pred<test$S

  cost=sum(alpha[1]*index_pos*pmax((test$S-test$Y),0)+alpha[2]*index_neg*pmax((test$Y-test$S),0))/length(test$A)
  misclass=(alpha[1]*sum(index_pos*(test$Y<test$S))+alpha[2]*sum(index_neg*(test$Y>test$S)))/length(test$A)
  return(list(misclass=misclass,cost=cost1,pred_test=index_pos,fitted=cvglmfit))

  for(i in 1:nfolds){
    cat(paste("=== Fold ",i,"===",sep=" "))
    index = seq(((i-1)*nsize+1),(i*nsize))
    train0 = list(X=train$X[-index,],A=train$A[-index],Y=train$Y[-index],opt=train$opt[-index],mu=train$mu[-index],S=train$S,propensity=train$propensity[-index])
    test0 = list(X=train$X[index,],A=train$A[index],Y=train$Y[index],opt=train$opt[index],mu=train$mu[index],S=train$S,propensity=train$propensity[index])
    fit.lasso=fit.lasso.model(train0,test0,alpha,nfolds=nfolds,model)
    Pred_test[index]=fit.lasso$pred_test
    misclass=misclass+fit.lasso$misclass
    cost=cost+fit.lasso$cost
  }
  misclass=misclass/nfolds
  cost=cost/nfolds

  #plot.FDI(train,Pred_test,main=paste("CV",main))

  fit.lasso=fit.lasso.model(train,test,alpha,nfolds=nfolds,model)
  Pred_test_ind=fit.lasso$pred_test
  misclass_ind=fit.lasso$misclass
  cost_ind=fit.lasso$cost

  #plot.FDI(test,Pred_test_ind,main=paste("independent",main))

  indirectfit=predict_indirect(fitted=fit.lasso$fitted,test=test,two.sided=two.sided,dt=0.1)

  return(list(bestPara=fit.lasso$fitted$lambda.min,misclass=misclass,misclass_ind=misclass_ind,misclass_search=indirectfit$misclass,correlation=indirectfit$correlation,Rsquare=indirectfit$Rsquare,percentage=indirectfit$percentage,Pred_test=Pred_test,Pred_test_ind=Pred_test_ind,cost=cost,cost_ind=cost_ind))
}

cv.SVR<-function(train,test,nfolds,alpha,Cost=c(0.1,1,3,5,8,10),lower=TRUE,two.sided=FALSE,...){
  n=length(train$A)
  p=dim(train$X)[2]
  nsize=n/nfolds
  Pred_test=rep(0,n)
  misclass=0
  cost=0
  prederror=Inf

  trainX = cbind(train$A,train$X)
  testX = cbind(test$A,test$X)

  for(k in 1:length(Cost)){
    tmpmodel = kernlab::ksvm(trainX,train$Y,cross = nfolds,C=Cost[k],kernel = "rbfdot",kpar="automatic")
    print(paste(Cost[k],tmpmodel@error))
    if(tmpmodel@error < prederror){
      svmfit = tmpmodel
      prederror = svmfit@error
    }
  }
  pred=predict(svmfit,newdata=testX)
  pred_test=pred>test$S
  misclass=(alpha[1]*sum((pred_test)*(test$Y<test$S))+alpha[2]*sum((!pred_test)*(test$Y>test$S)))/length(test$A)
  cost=(alpha[1]*sum((pred_test)*pmax((test$S-test$Y),0))+alpha[2]*sum((!pred_test)*pmax((test$Y-test$S),0)))/length(test$A)

  indirectfit=searchInterval(fitted=fit.SVR$fitted,test=test,two.sided=two.sided,dt=0.1)

  return(list(bestPara=fit.SVR$fitted@param$C,misclass=misclass,misclass_ind=misclass_ind,misclass_search=indirectfit$misclass,correlation=indirectfit$correlation,Rsquare=indirectfit$Rsquare,percentage=indirectfit$percentage,Pred_test=Pred_test,Pred_test_ind=Pred_test_ind,cost=cost,cost_ind=cost_ind,indirectfit=indirectfit))
}

cv.foreg<-function(train,test,nfolds,alpha,two.sided,Gamma=c(0.2,0.3,0.4,0.5),main="randomForest Regression",...){
  n=length(train$A)
  length(test$A)=length(test$A)
  p=dim(train$X)[2]
  Pred_test=rep(0,n)
  misclass=0
  cv.forest=e1071::tune.randomForest(x=cbind(train$A,train$X),y=train$Y,mtry=round(Gamma*p+1),tunecontrol=tune.control(sampling = "cross",cross=nfolds))
  fit.forest=randomForest(cbind(train$A,train$X),train$Y,mtry=as.numeric(cv.forest$best.parameters))

  indirectfit=predict_indirect(fitted=fit.forest,test=test,two.sided=two.sided,dt=0.1)

  svrd=predict(fit.forest)>train$S
  cost=(alpha[1]*sum((svrd)*pmax((train$S-train$Y),0))+alpha[2]*sum((!svrd)*pmax((train$Y-train$S),0)))/n
  misclass=(alpha[1]*sum((svrd)*(train$Y<train$S))+alpha[2]*sum((!svrd)*(train$Y>train$S)))/n
  Pred_test=svrd

  svrd=predict(fit.forest,newtrain=cbind(test$A,test$X))>test$S
  cost_ind=(alpha[1]*sum((svrd)*pmax((test$S-test$Y),0))+alpha[2]*sum((!svrd)*pmax((test$Y-test$S),0)))/length(test$A)
  misclass_ind=(alpha[1]*sum((svrd)*(test$Y<test$S))+alpha[2]*sum((!svrd)*(test$Y>test$S)))/length(test$A)
  Pred_test_ind=svrd

  #plot.FDI(train,Pred_test,main=paste("CV",main))
  #plot.FDI(test,Pred_test,main=paste("independent",main))
  return(list(bestPara=as.numeric(cv.forest$best.parameters),misclass=misclass,misclass_ind=misclass_ind,misclass_search=indirectfit$misclass,correlation=indirectfit$correlation,Rsquare=indirectfit$Rsquare,percentage=indirectfit$percentage,Pred_test=Pred_test,Pred_test_ind=Pred_test_ind,cost=cost,cost_ind=cost_ind,indirectfit=indirectfit))

  #return(list(misclass=misclass,misclass_ind=misclass_ind,cost=cost,cost_ind=cost_ind,Pred_test=Pred_test,Pred_test_ind=Pred_test))
}

cv.logist<-function(train,test,nfolds,alpha,two.sided=FALSE,lower=TRUE,Lambda=c(2^(-(0:10)),0.00000001),...){
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
  misclass=sum(alpha[1]*(!pred)*(test$Y>test$S)+alpha[2]*(pred)*(test$Y<test$S))/length(test$Y)
  cost=sum(alpha[1]*(pred)*pmax((test$S-test$Y),0)+alpha[2]*(!pred)*pmax((test$Y-test$S),0))/length(test$Y)

  search=searchInterval(fitted=cvglmfit,test=test,two.sided=two.sided,alpha=alpha)

  return(list(fitted=cvglmfit,misclass=misclass,cost=cost,search=search,pred=pred))
}


cv.SVM<-function(train,test,nfolds,alpha,two.sided=FALSE,Cost=2^(0:15),lower=TRUE,...){
  n=length(train$A)
  p=dim(train$X)[2]
  nsize=n/nfolds
  misclass<-NA->cost

  trainX = cbind(train$A,train$X)
  testX = cbind(test$A,test$X)

  for(k in 1:length(Cost)){
    tmpmodel = kernlab::ksvm(trainX,train$Y>train$S,type='C-svc',cross = nfolds,C=Cost[k],kernel = "rbfdot",kpar="automatic",prob.model=TRUE)
    error=mean((predict(tmpmodel)-1)!=(train$Y>train$S))
    #print(paste(Cost[k],error))
    if(tmpmodel@error < prederror){
      svmfit = tmpmodel
      prederror = error
      bestPara=Cost[k]
    }
  }
  pred=predict(svmfit,newdata=testX,type='response')>test$S
  misclass=(alpha[1]*sum((pred)*(test$Y<test$S))+alpha[2]*sum((!pred)*(test$Y>test$S)))/length(test$A)
  cost=(alpha[1]*sum((pred)*pmax((test$S-test$Y),0))+alpha[2]*sum((!pred)*pmax((test$Y-test$S),0)))/length(test$A)

  search=searchInterval(fitted=svmfit,test=test,two.sided=two.sided,alpha=alpha)
  return(list(bestPara=bestPara,misclass=misclass,cost=cost,search=search,pred=pred))
}




cv.forest<-function(train,test,nfolds,alpha=c(0.5,0.5),two.sided=FALSE,lower=TRUE,Mtry=c(0.2,0.3,0.4),...){
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
      fit.forest=randomForest(Dat,as.factor(train0$Y>train0$S),mtry=round(Mtry[k]*p+1))
      Dat=cbind(test0$A,test0$X)
      colnames(Dat)=NULL
      pred=ifelse(predict(fit.forest,Dat)==TRUE,TRUE,FALSE)
      Perform[k,1]=  Perform[k,1]+sum(alpha[2]*(!pred)*(test0$Y>test0$S)+alpha[1]*(pred)*(test0$Y<test0$S))/length(test0$Y)
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
  pred=ifelse(predict(fit.forest,Dat)==TRUE,TRUE,FALSE)
  misclass=sum(alpha[2]*(!pred)*(test$Y>test$S)+alpha[1]*(pred)*(test$Y<test$S))/length(test$Y)
  cost=sum(alpha[1]*(pred)*pmax((test$S-test$Y),0)+alpha[2]*(!pred)*pmax((test$Y-test$S),0))/length(test$Y)

  search=searchInterval(fitted=fit.forest,test=test,two.sided=two.sided,alpha=alpha)
  return(list(bestPara=Perform[id,2],misclass=misclass,cost=cost,search=search,pred=pred))
}
